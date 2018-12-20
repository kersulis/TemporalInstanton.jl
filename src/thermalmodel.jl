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
    "[pu] line resistance"
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

"""
Return (c1, c4, c5, %()c7) thermal constants used in line temperature IVP.
Arguments:

* `lp` instance of LineParams
* `cp` instance of ConductorParams
* `Tamb` [C] ambient temperature
"""
function return_thermal_constants(
    lp::LineParams,
    cp::ConductorParams,
    Tamb::Float64,
    Sb::Float64,
    int_length::Float64,
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
    return mCp \ (-ηc - 4 * ηr * (Tmid + 273)^3)
end

"""
Return constant `c4` [W/m]. Arguments:

* `mCp` [J/m-C] is line heat capacity
* `r` [pu] is line resistance
* `x` [pu] is line reactance
* `Sb` [W] is system base MVA
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
    il::Float64,
    c1::Float64,
    c5::Float64,
    n::Int64,
    T0::Float64
    )
    sum_coeff = sum([(e^(il * c1))^i - (e^(il * c1))^(i - 1) for i in 1:n])
    return (e^(il * c1))^n * T0 + (c5 / c1) * sum_coeff
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
