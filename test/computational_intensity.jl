using Distributed
addprocs(2)

@everywhere using TemporalInstanton, LineThermalModel, JLD, PowerModels, Memento

setlevel!(getlogger(PowerModels), "error")

const PGLIB_PATH = "/home/jk/jdocuments/projects/sdp-opf-decomp/pglib-opf"

##

filenames = [
    "pglib_opf_case14_ieee.m";
    "pglib_opf_case24_ieee_rts.m";
    "pglib_opf_case30_ieee.m";
    "pglib_opf_case39_epri.m";
    "pglib_opf_case57_ieee.m";
    "pglib_opf_case73_ieee_rts.m";
    "pglib_opf_case89_pegase.m";
    "pglib_opf_case118_ieee.m";
    "pglib_opf_case162_ieee_dtc.m";
    "pglib_opf_case240_pserc.m";
    "pglib_opf_case300_ieee.m";
    "pglib_opf_case500_tamu.m";
    "pglib_opf_case588_sdet.m";
    "pglib_opf_case1354_pegase.m";
    "pglib_opf_case1888_rte.m";
    "pglib_opf_case1951_rte.m";
    "pglib_opf_case2383wp_k.m";
    "pglib_opf_case3120sp_k.m";
    # "pglib_opf_case4661_sdet.m";
    # "pglib_opf_case6515_rte.m"
]

# precompile
fpath = joinpath(PGLIB_PATH, "pglib_opf_case30_ieee.m")
nd = PowerModels.parse_file(fpath)

i = build_instanton_input(fpath)
conventional_to_renewable!(i)
set_temperatures!(i; Tamb=40.0)

# one hour in ten-minute intervals
set_timing!(i, 0:600:3600)

# nominal injections and loads
Gp, Dp, Rp = i.G0, i.D0, i.R0
Gp = Gp ./ sum(Gp)
Dp = Dp ./ sum(Dp)
Rp = Rp ./ sum(Rp)

Gp, Dp, Rp = i.G0, i.D0, i.R0
set_timing!(i, 0:600:600)
G0 = 1.0 * Gp
D0 = 1.0 * Dp
R0 = 1.0 * Rp

set_injections!(i; G0=G0, D0=D0, R0=R0)

# set initial line temperatures based on initial gen/load, by doing
# instantaneous power flow to get an angle difference, and assuming
# that heat input to the line is constant forever.
# Use analytic asymptotic value of heat balance equation.
set_temperatures!(i; Tamb=40.0, Tlim=fill(100.0, length(i.lines)))

o = solve_temporal_instanton(i)

## loop over all specified PGLIB systems

linetimes = Dict{String, Vector{Float64}}()
dimension = Dict{String, Tuple}()
times = Dict{String, Float64}()

for fname in filenames
    fpath = joinpath(PGLIB_PATH, fname)

    nd = PowerModels.parse_file(fpath)

    i = build_instanton_input(fpath)
    conventional_to_renewable!(i)
    set_temperatures!(i; Tamb=40.0)

    # let final temperature be initial temp plus 10 deg
    set_temperatures!(i; Tamb=40.0, Tlim=i.T0 .+ 10.0)

    # half hour in ten-minute intervals
    set_timing!(i, 0:600:1800)

    # nominal injections and loads
    Gp, Dp, Rp = i.G0, i.D0, i.R0
    # Gp = Gp ./ sum(Gp)
    # Dp = Dp ./ sum(Dp)
    # Rp = Rp ./ sum(Rp)

    G0 = 1.0 * [1.0 * Gp; 1.0 * Gp; 1.0 * Gp] #; 1.0 * Gp; 1.0 * Gp; 1.0 * Gp]
    D0 = 1.0 * [1.0 * Dp; 1.0 * Dp; 1.0 * Dp] #; 1.0 * Dp; 1.0 * Dp; 1.0 * Dp]
    R0 = 1.0 * [1.0 * Rp; 1.0 * Rp; 1.0 * Rp] #; 1.0 * Rp; 1.0 * Rp; 1.0 * Rp]

    set_injections!(i; G0=G0, D0=D0, R0=R0)

    t = time()
    o = solve_temporal_instanton(i; maxlines=100)
    times[fname] = time() - t

    n = size(i.Y, 1)
    nr = length(i.Ridx)
    T = length(i.time_values) - 1

    linetimes[fname] = o.linetimes
    dimension[fname] = (n, nr, T)
    println("finished $(fname)")
end

##
using Plots
using Statistics: mean

gr(;
    markerstrokecolor="white",
    label="",
    markersize=5,
    markercolor="black",
    alpha=0.7
)
n = [dimension[k][1] for k in filenames]
nr = [dimension[k][2] for k in filenames]
T = [dimension[k][3] for k in filenames]

dim = (n .+ nr .+ 2) .* T

mean_linetime = mean.([linetimes[k] for k in filenames])
time_vec = [times[k] for k in filenames]

# scatter(n, mean_linetime; yscale=:log10, xscale=:log10)

##
scatter(
    dim, mean_linetime .* 1e3;
    yscale=:log10, xscale=:log10,
    xlabel="Variables (log)",
    ylabel="Mean QCQP solution time (ms, log)"
)

##
scatter(
    dim, time_vec ./ 60;
    xscale=:log10, yscale=:log10,
    xlabel="Variables (log)", ylabel="Total solution time (min)"
)

##
using JLD

JLD.save("/home/jk/.julia/dev/TemporalInstanton/2019-06-27-timing.jld", "linetimes", linetimes, "dimension", dimension, "times", times)
