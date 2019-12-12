using LineThermalModel

"""
Line parameters that vary with the line (not inherent to conductor material).
"""
struct LineParams
    "Line 'from' node"
    from::Int64
    "Line 'to' node"
    to::Int64
    "pu line resistance"
    rij::Float64
    "pu line reactance"
    xij::Float64
    "Line length in meters"
    length::Float64
end

"""
Line parameters that depend on conductor material.
"""
struct ConductorParams
    "[m] conductor diameter"
    D0::Float64
    "[J/m-C] line heat capacity"
    mCp::Float64
    "[A] max. allowable current"
    Ilim::Float64
    "[Ohms/m] conductor resistance"
    r::Float64
    "[C] highest allowable line temperature"
    Tlim::Float64
    "[W/m-C] conductive heat loss rate coefficient"
    ηc::Float64
    "[W/m-C^4] radiative heat loss rate coefficient"
    ηr::Float64
    "[W/m] solar heat gain rate"
    qs::Float64
end

"""
Returns instance of `ConductorParams` filled with
conductor parameters. Data from Mads's MPC paper.

Currently accepts "waxwing" and "dove" as arguments.
"""
function return_conductor_params(conductor::String)
    if conductor == "waxwing"
        return ConductorParams(
            15.5e-3, 383.0, 439.0, 110e-6, 65.0, 0.955, 2.207e-9, 14.4
        )
    elseif conductor == "dove"
        return ConductorParams(
            23.5e-3, 916.0, 753.0, 60e-6, 69.0, 1.179, 3.346e-9, 21.9
        )
    end
end

function return_conductor_params(acsr_spec::ACSRSpecsMetric, Ilim::Float64, T_a::Float64, T_s::Float64, Tlim::Float64)
#     return ConductorParams(
#         15.5e-3, 383.0, 439.0, 110e-6, 65.0, 0.955, 2.207e-9, 14.4
#     )
# end
    mCp = eq_mCp(acsr_spec)

    emm     = 0.7  # emissivity, [0.23, 0.91]
    phi     = 90.0 # wind/line angle in degrees: wind perpendicular to line.
    H_e     = 61.0 # height above sea level in m. Avg PJM elevation.
    V_w     = 0.61 # wind speed in m/s
    alpha   = 0.9  # solar absorptivity, [0.23, 0.91]
    lat     = 40.0 # latitude in deg
    N       = 161  # day of year: June 10
    Z_l     = 90.0 # line azimuth: West-to-East
    hours_from_noon = 0.0 # time: noon

    D, Al_m, St_m, R, bundle, label = acsr_spec.D, acsr_spec.Al_m, acsr_spec.St_m, acsr_spec.R, acsr_spec.bundle, acsr_spec.label

    A_prime = D # conductor area, m^2 per linear m
    ηr      = eq_eta_r(D, emm)
    T_film  = eq6_T_film(T_s, T_a)
    k_f     = eq15a_k_f(T_film)
    K_angle = eq4a_K_angle(phi)
    p_f     = eq14a_p_f(H_e, T_film)
    mu_f    = eq13a_mu_f(T_film)
    N_Re    = eq2c_N_Re(D, p_f, V_w, mu_f)
    ηc      = eq_eta_c(k_f, K_angle, N_Re)

    omega   = eq_omega(hours_from_noon)
    delta   = eq16b_delta(N)
    H_c     = eq16a_H_c(lat, delta, omega)
    chi     = eq17b_chi(omega, lat, delta)
    C       = eqtable2_C(omega, chi)
    Z_c     = eq17a_Z_c(C, chi)
    theta   = eq9_theta(H_c, Z_c, Z_l)
    Q_s     = eq18_Q_s(H_c)
    K_solar = eq20_K_solar(H_e)
    Q_se    = eq19_Q_se(K_solar, Q_s)
    qs      = eq8_q_s(alpha, Q_se, theta, A_prime)

    return ConductorParams(D, mCp, Ilim, R, Tlim, ηc, ηr, qs)
end

"""
Return (c1, c4, c5, c7) thermal constants used in line temperature IVP.
Arguments:

* `lp` instance of LineParams
* `cp` instance of ConductorParams
* `Tamb` [C] ambient temperature
* `Sb` system base MVA
* `int_length` [s] length of time interval
* `n` number of time intervals
* `T0` [C] initial line temperature
"""
function return_thermal_constants(
    lp::LineParams,
    cp::ConductorParams,
    Tamb::Float64,
    Sb::Float64,
    int_length,
    n,
    T0::Float64
    )
    c1 = compute_c1(cp.mCp, cp.ηc, cp.ηr, Tamb, cp.Tlim)
    c4 = compute_c4(cp.mCp, lp.rij, lp.xij, Sb, lp.length)
    c5 = compute_c5(cp.mCp, cp.ηc, cp.ηr, Tamb, cp.Tlim, cp.qs)
    c7 = compute_c7(int_length, c1, c5, n, T0)
    return c1, c4, c5, c7
end

"""
Returns constant `c1` [1/s]. Arguments:

* `mCp` [J/m-C] is line heat capacity
* `ηc` [W/m-C] is conductive heat loss rate coefficient
* `ηr` [W/m-C^4] is radiative heat loss rate coefficient
* `Tamb` [C] is ambient temperature (of air)
* `Tlim` [C] is highest allowable line temperature
"""
function compute_c1(
    mCp::Float64,
    ηc::Float64,
    ηr::Float64,
    Tamb::Float64,
    Tlim::Float64
    )
    Tmid = (Tamb + Tlim) / 2
    return -(ηc + 4 * ηr * (Tmid + 273)^3) / mCp
end

"""
Return constant `c4` [W/m]. Arguments:

* `mCp` [J/m-C] is line heat capacity
* `r` [pu] is line resistance
* `x` [pu] is line reactance
* `Sb` [VA] is system base
* `L` [m] is line length
"""
function compute_c4(
    mCp::Float64,
    r::Float64,
    x::Float64,
    Sb::Float64,
    L::Float64
    )
    return r * Sb / (3 * mCp * L * x^2)
end

"""
Returns constant `c5` [W/m]. Arguments:

* `mCp` [J/m-C] is line heat capacity
* `ηc` [W/m-C] is conductive heat loss rate coefficient
* `ηr` [W/m-C^4] is radiative heat loss rate coefficient
* `Tamb` [C] is ambient temperature (of air)
* `Tlim` [C] is highest allowable line temperature
* `q_solar` [W/m] is the solar heat gain rate
"""
function compute_c5(
    mCp::Float64,
    ηc::Float64,
    ηr::Float64,
    Tamb::Float64,
    Tlim::Float64,
    q_solar::Float64
    )
    Tmid = (Tamb + Tlim) / 2
    TmidK = Tmid + 273
    return mCp \ (ηc * Tamb - ηr * ((TmidK)^4 - (Tamb + 273)^4) + 4 * ηr * Tmid * (TmidK)^3 + q_solar)
end

"""
Returns constant `c7` [C]. Arguments:

* `il` [s] is time interval length
* `c1` [1/s] is a constant
* `c5` [W/m] is a constant
* `n` [unitless] is the number of time intervals
* `T0` [C] is the initial steady-state line temp
"""
function compute_c7(
    il,
    c1::Float64,
    c5::Float64,
    n::Int64,
    T0::Float64
    )
    sum_coeff = sum([(ℯ^(il * c1))^i - (ℯ^(il * c1))^(i - 1) for i in 1:n])
    return (ℯ^(il * c1))^n * T0 + (c5 / c1) * sum_coeff
end

"""
Return line's final temperature. Arguments:

* `int_length` [s] is length of each interval
* `c1` [1/s] is a constant
* `c4` [W/m] is a constant
* `c7` [C] is a constant
* `n` [unitless] is the number of time intervals
* `θij` [rad] is the vector of angle differences (sorted by time interval)
"""
function compute_T(
    int_length::Float64,
    c1::Float64,
    c4::Float64,
    c7::Float64,
    n::Int64,
    θij::Float64
    )
    angle_coeffs = [(e^(int_length * c1))^i - (e^(int_length * c1))^(i - 1) for i in 1:n]
    return c7 + (c4 / c1) * dot(angle_coeffs, flipud(θij).^2)
end

"""
    temp_values = temperature_trajectory(i, o, idx[, time_steps=10, fixed_wind=[]])

Return temperature trajectory corresponding to instance of
`InstantonInput` and `InstantonOutput`.
`idx` is the event index, with 1 denoting the first set of data in o.
(With default sorting, this will be the instanton.)

This function is useful mostly as a means of validating the temporal deviation scanning algorithm. The first value of temp_values should match the line's initial temperature in i. The final value should match the temperature limit.
"""
function temperature_trajectory(
    i::InstantonInput,
    o::InstantonOutput,
    idx,
    time_steps=10,
    fixed_wind=[]
    )

    i_idx = o.analytic_lines[idx]

    # Line parameters:
    line = i.lines[i_idx]
    from, to = line
    r_ij = i.res[i_idx]
    x_ij = i.reac[i_idx]
    L_ij = i.line_lengths[i_idx]

    Tamb = i.Tamb
    T0 = i.T0[i_idx]
    Tlim = i.Tlim[i_idx]
    cond_params = return_conductor_params(i.line_conductors[i_idx], i.current_limits[i_idx], Tamb, T0, Tlim)
    mCp, ηc, ηr, qs, Tlim = cond_params.mCp, cond_params.ηc, cond_params.ηr, cond_params.qs, cond_params.Tlim
    Tmid = (Tamb + Tlim) / 2

    n = size(i.Y, 1)
    nr = length(i.Ridx)
    numSteps = length(i.time_values) - 1

    Sb = i.Sb

    "Solution to approx. heat balance IVP"
    temp_eq(t, T0, a, b) = (T0 + b / a) * exp(a * t) - b / a

    # Fixed wrt power flow
    therm_a = mCp \ (-ηc - 4 * ηr * (Tmid + 273)^3)

    temp_trajectory = Vector()

    if isempty(fixed_wind)
        fixed_wind = vcat(o.x[idx]...)
    end

    # DC power flow
    fixed_A = fixed_wind_A(numSteps, i.Y, i.ref, i.k)
    fixed_P = expand_renewable_vector(fixed_wind, i.Ridx, n, numSteps)
    fixed_b = fixed_wind_b(n, numSteps, i.G0, i.R0 + fixed_P, i.D0)
    fixed_x = fixed_A \ fixed_b
    angles, alpha = return_angles(fixed_x, n, numSteps)

    fixed_diffs = return_angle_diffs(angles, line)

    temp_values = [T0]
    power_flow = Float64[]

    eval_times = range(0; stop=i.time_values.step, length=time_steps)

    for θij in fixed_diffs
        f_loss_pu = r_ij * (θij / x_ij)^2 # pu
        f_loss_si = f_loss_pu * Sb / (3 * L_ij) # W/m
        push!(power_flow, (Sb / 1e6) * θij / x_ij) # MW
        therm_b = mCp \ (f_loss_si + ηc * Tamb - ηr * ((Tmid + 273)^4 -
            (Tamb + 273)^4) + 4 * ηr * Tmid * (Tmid + 273)^3 + qs)
        temp_values = [temp_values; temp_eq.(eval_times, T0, therm_a, therm_b)[2:end]]
        T0 = temp_values[end]
    end
    return temp_values
end

"""
    ss_temps = steady_state_temps(i, fixed_wind=[])

Return steady state temperatures corresponding to line data, generation, and load encoded in i::InstantonInput. Optionally pass in `fixed_wind` to specify more or less wind generation.
"""
function steady_state_temps(i::InstantonInput, fixed_wind=[])
    T0 = 50.0 # arbitrary, as we will find asymptotic temp
    Tamb = i.Tamb
    Sb = i.Sb

    n = size(i.Y, 1)
    nr = length(i.Ridx)

    ss_temps = Float64[]
    for (line_idx, line) in enumerate(i.lines)
        from, to = line
        r_ij = i.res[line_idx]
        x_ij = i.reac[line_idx]
        L_ij = i.line_lengths[line_idx]
        Tlim = i.Tlim[line_idx]

        cond_params = return_conductor_params(i.line_conductors[line_idx], i.current_limits[line_idx], Tamb, T0, Tlim)
        mCp, ηc, ηr, qs, Tlim = cond_params.mCp, cond_params.ηc, cond_params.ηr, cond_params.qs, cond_params.Tlim
        Tmid = (Tamb + Tlim) / 2

        # Fixed wrt power flow
        therm_a = mCp \ (-ηc - 4 * ηr * (Tmid + 273)^3)

        if isempty(fixed_wind)
            # assume 0 deviation from forecast
            fixed_wind = zeros(nr)
        end

        fixed_A = fixed_wind_A(1, i.Y, i.ref, i.k)
        fixed_P = expand_renewable_vector(fixed_wind, i.Ridx, n, 1)

        # Use initial injections to compute steady-state temps
        fixed_b = fixed_wind_b(n, 1, i.G0[1:n], i.R0[1:n] + fixed_P, i.D0[1:n])
        fixed_x = fixed_A \ fixed_b
        angles, alpha = return_angles(fixed_x, n, 1)
        fixed_diffs = return_angle_diffs(angles, line)

        θij = fixed_diffs[1]
        f_loss_pu = r_ij * (θij / x_ij)^2 # pu
        f_loss_si = f_loss_pu * Sb / (3 * L_ij) # W/m
        # push!(power_flow, (Sb / 1e6) * θij / x_ij)
        therm_b = mCp \ (f_loss_si + ηc * Tamb - ηr * ((Tmid + 273)^4 -
            (Tamb + 273)^4) + 4 * ηr * Tmid * (Tmid + 273)^3 + qs)
        push!(ss_temps, -therm_b / therm_a)
    end
    return ss_temps
end

"""
    temp_trajectories = temperature_trajectories(i, o, event_idx[, time_steps=10])

Inputs:
- `i::InstantonInput`
- `d::Vector{Vector{Float64}}`, a deviation pattern similar to an element of the `x` field of `InstantonOutput`
- `time_steps`: number of evenly-spaced points at which to evaluate the line temperature equation per time step.

Output: dictionary mapping line tuples to temperature trajectories.

*Return temperature trajectories for all lines in a power system, given an instance of InstantonInput and a deviation pattern.
"""
function temperature_trajectories(
    i::InstantonInput,
    d::Vector{Vector{Float64}},
    time_steps::Int64=10
    )

    "Solution to approx. heat balance IVP"
    temp_eq(t, T0, a, b) = (T0 + b / a) * exp(a * t) - b / a

    eval_times = range(0; stop=i.time_values.step, length=time_steps)

    n = size(i.Y, 1)
    nr = length(i.Ridx)
    numSteps = length(i.time_values) - 1
    Tamb = i.Tamb
    Sb = i.Sb

    # i_idx = o.analytic_lines[event_idx]
    fixed_wind = vcat(d...)

    # DC power flow: obtain all angles
    fixed_A = fixed_wind_A(numSteps, i.Y, i.ref, i.k)
    fixed_P = expand_renewable_vector(fixed_wind, i.Ridx, n, numSteps)
    fixed_b = fixed_wind_b(n, numSteps, i.G0, i.R0 + fixed_P, i.D0)
    fixed_x = fixed_A \ fixed_b
    angles, alpha = return_angles(fixed_x, n, numSteps)

    # loop over all lines
    temp_trajectories = Vector{Vector{Float64}}()
    for (idx, line) in enumerate(i.lines)
        # Line parameters
        from, to = line
        r_ij = i.res[idx]
        x_ij = i.reac[idx]
        L_ij = i.line_lengths[idx]
        T0 = i.T0[idx]
        Tlim = i.Tlim[idx]

        # Thermal parameters
        cp = return_conductor_params(i.line_conductors[idx], i.current_limits[idx], Tamb, T0, Tlim)
        mCp, ηc, ηr, qs, Tlim = cp.mCp, cp.ηc, cp.ηr, cp.qs, cp.Tlim
        Tmid = (Tamb + Tlim) / 2
        therm_a = mCp \ (-ηc - 4 * ηr * (Tmid + 273)^3)

        # Angle differences
        angle_diffs = return_angle_diffs(angles, line)

        traj = [T0]
        for θij in angle_diffs
            f_loss_pu = r_ij * (θij / x_ij)^2 # pu
            f_loss_si = f_loss_pu * Sb / (3 * L_ij) # W/m
            # push!(power_flow, (Sb / 1e6) * θij / x_ij) # MW
            therm_b = mCp \ (f_loss_si + ηc * Tamb - ηr * ((Tmid + 273)^4 -
                (Tamb + 273)^4) + 4 * ηr * Tmid * (Tmid + 273)^3 + qs)
            append!(traj, temp_eq.(eval_times, T0, therm_a, therm_b)[2:end])
            T0 = traj[end]
        end
        push!(temp_trajectories, traj)
    end
    return temp_trajectories
end

"""
    temp_trajectories = temperature_trajectories(i, o, event_idx[, time_steps=10])

Inputs:
- `i::InstantonInput`
- `o::InstantonOutput`
- `event_idx`: 1 indicates lowest deviation scanning objective score, and so on
- `time_steps`: number of evenly-spaced points at which to evaluate the line temperature equation per time step.

Output: dictionary mapping line tuples to temperature trajectories.

*Return temperature trajectories for all lines in a power system, given instances of InstantonInput and InstantonOutput, and an event index (which indexes into o). Only one line (the one corresponding to the event index) should reach its limit temperature.*
"""
function temperature_trajectories(
    i::InstantonInput,
    o::InstantonOutput,
    event_idx::Int64;
    time_steps::Int64=10
    )
    d = o.x[event_idx]
    return temperature_trajectories(i, d, time_steps)
end
