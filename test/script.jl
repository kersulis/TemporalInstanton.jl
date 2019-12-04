# using Distributed
# addprocs(2)

# @everywhere using TemporalInstanton, LineThermalModel, JLD, PowerModels, Memento

using TemporalInstanton, JLD, PowerModels, Memento

setlevel!(getlogger(PowerModels), "error")

# fpath = joinpath(@__DIR__, "data", "caseRTS96.jld")
# JLD.load(fpath)

fpath = joinpath(@__DIR__, "data", "nesta_case118_ieee.m")
# fpath = joinpath(@__DIR__, "data", "pglib_opf_case30_ieee.m")
# fpath = joinpath(@__DIR__, "data", "pglib_opf_case5_pjm.m")
# fpath = joinpath("/home/jk/jdocuments/projects/sdp-opf-decomp/pglib-opf", "pglib_opf_case300_ieee.m")

# fpath = joinpath("/home/jk/jdocuments/projects/sdp-opf-decomp/pglib-opf", "pglib_opf_case2383wp_k.m")

nd = PowerModels.parse_file(fpath)

i = build_instanton_input(fpath)
conventional_to_renewable!(i)
# set_temperatures!(i; Tamb=40.0)

# one hour in ten-minute intervals
set_timing!(i, 0:600:3600)

# nominal injections and loads
Gp, Dp, Rp = i.G0, i.D0, i.R0
Gp = Gp ./ sum(Gp)
Dp = Dp ./ sum(Dp)
Rp = Rp ./ sum(Rp)

G0 = 44.0 * [1.0 * Gp; 1.0 * Gp; 1.0 * Gp; 1.0 * Gp; 1.0 * Gp; 1.0 * Gp]
D0 = 65.0 * [1.0 * Dp; 1.0 * Dp; 1.0 * Dp; 1.0 * Dp; 1.0 * Dp; 1.0 * Dp]
R0 = 10.0 * [1.0 * Rp; 1.0 * Rp; 1.0 * Rp; 1.0 * Rp; 1.0 * Rp; 1.0 * Rp]


# Gp, Dp, Rp = i.G0, i.D0, i.R0
# set_timing!(i, 0:600:600)
# G0 = 1.0 * Gp
# D0 = 1.0 * Dp
# R0 = 1.0 * Rp


# set_timing!(i, 0:600:1200)
# # nominal injections and loads
# Gp, Dp, Rp = i.G0, i.D0, i.R0
#
# G0 = 1.0 * [1.0 * Gp; 1.0 * Gp]
# D0 = 1.0 * [1.0 * Dp; 1.0 * Dp]
# R0 = 1.0 * [1.0 * Rp; 1.0 * Rp]

set_injections!(i; G0=G0, D0=D0, R0=R0)

# set initial line temperatures based on initial gen/load, by doing
# instantaneous power flow to get an angle difference, and assuming
# that heat input to the line is constant forever.
# Use analytic asymptotic value of heat balance equation.
set_temperatures!(i; Tamb=40.0, Tlim=fill(100.0, length(i.lines)))

i.auto_prec = 0.03

o = solve_temporal_instanton(i)
o.score
# include("../src/plot.jl")
# ss_temps = steady_state_temps(i)
# i.T0 |> sort

## temperature trajectories

time_steps = 10
fixed_wind = []
idx = 1

# import TemporalInstanton.temperature_trajectory
# include("../src/plot.jl")
# include("../src/powerflow.jl")
traj = temperature_trajectory(i, o, idx, time_steps)
tt = range(0; stop=i.time_values[end], length=length(traj))

using Plots
p = plot(tt, traj; color=:red)

trajs = temperature_trajectories(i, o, idx)
for ti in values(trajs)
    plot!(tt, ti; color=:gray, alpha=0.7)
end
p
## SVD of set of all deviation vectors
using Plots
gr(
    markerstrokewidth=0.3,
    markerstrokecolor=:white,
    markersize=6,
    label=""
)
using LinearAlgebra

normalize(x::Vector) = x / norm(x)

X = hcat((vcat(xi...) |> normalize for xi in o.x)...)

U, s, V = svd(X)

[norm(X[:, i]) for i in 1:size(X, 2)]

scatter(1:length(s), s)

# maximum((U * Diagonal(S) * V') - X)
##

"""
Estimate heat inputs for all lines in a network by simulating
heating of each line using approximate thermal model developed by Kersulis based on Mads line loss approximation.

Inputs:
* `network_data`: output dict of PowerModels.parse_file
* `T_a`: ambient air temp, Celsius

Outputs:
* `convection`: convective heat rate coefficients
* `radiation`: radiation heat loss rate coefficients
* `insolation`: solar heat rate coefficients
"""
function estimate_temperature_limits(network_data::Dict{String, Any}, T_a::Float64=40.0)

    acsr_specs = acsr_interpolation(network_data)

    emm     = 0.7  # emissivity, [0.23, 0.91]
    phi     = 90.0 # wind/line angle in degrees: wind perpendicular to line.
    H_e     = 61.0 # height above sea level in m. Avg PJM elevation.
    V_w     = 0.61 # wind speed in m/s
    alpha   = 0.9  # solar absorptivity, [0.23, 0.91]
    lat     = 40.0 # latitude in deg
    N       = 161  # day of year: June 10
    Z_l     = 90.0 # line azimuth: West-to-East
    hours_from_noon = 0.0 # time: noon

    T_s     = T_a # conductor surface temperature, C

    # compute heat rate coefficients for all lines
    convection = Dict{String, Float64}()
    radiation = Dict{String, Float64}()
    insolation = Dict{String, Float64}()
    # R and bundle are not necessary here
    for (id, a) in acsr_specs
        D, Al_m, St_m, R, bundle, label = a.D, a.Al_m, a.St_m, a.R, a.bundle, a.label

        A_prime = D # conductor area, m^2 per linear m
        eta_r   = eq_eta_r(D, emm)
        T_film  = eq6_T_film(T_s, T_a)
        k_f     = eq15a_k_f(T_film)
        K_angle = eq4a_K_angle(phi)
        p_f     = eq14a_p_f(H_e, T_film)
        mu_f    = eq13a_mu_f(T_film)
        N_Re    = eq2c_N_Re(D, p_f, V_w, mu_f)
        eta_c   = eq_eta_c(k_f, K_angle, N_Re)

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
        eta_s   = eq8_q_s(alpha, Q_se, theta, A_prime)

        convection[id] = eta_c
        radiation[id] = eta_r
        insolation[id] = eta_s
    end

    # now use heat rate coefficients to simulate line temps
    # long enough that they reach steady state
    return convection, radiation, insolation
end

estimate_temperature_limits(network_data)

#
# """
# Use bisection to iteratively compute steady-state temperature
# of a line with provided convective, radiative, and solar heat
# rates, and the given current.
#
# ISSUE: need to use R(Tavg), otherwise you can't compute temp
# limit this way.
# """
# function steady_state_temp(qc, qr, qs, I, R; tol=1e-6)
#     err = Inf
#     while err > tol
