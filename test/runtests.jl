using Test, TemporalInstanton, PowerModels, Memento

setlevel!(getlogger(PowerModels), "error")

# Set up and run a small basic test case

fpath = joinpath(@__DIR__, "data", "nesta_case118_ieee.m")
# fpath = joinpath(@__DIR__, "data", "pglib_opf_case30_ieee.m")
# fpath = joinpath(@__DIR__, "data", "pglib_opf_case5_pjm.m")
# fpath = joinpath("/home/jk/jdocuments/projects/sdp-opf-decomp/pglib-opf", "pglib_opf_case300_ieee.m")

# fpath = joinpath("/home/jk/jdocuments/projects/sdp-opf-decomp/pglib-opf", "pglib_opf_case2383wp_k.m")

nd = PowerModels.parse_file(fpath)

i = build_instanton_input(fpath)
conventional_to_renewable!(i)

# half hour in ten-minute intervals
set_timing!(i, 0:600:1800)

# nominal injections and loads
Gp, Dp, Rp = i.G0, i.D0, i.R0
Gp = Gp ./ sum(Gp)
Dp = Dp ./ sum(Dp)
Rp = Rp ./ sum(Rp)

G0 = 18.0 * [1.0 * Gp; 1.0 * Gp; 1.0 * Gp]
D0 = 30.0 * [1.0 * Dp; 1.0 * Dp; 1.0 * Dp]
R0 = 10.0 * [1.0 * Rp; 1.0 * Rp; 1.0 * Rp]

set_injections!(i; G0=G0, D0=D0, R0=R0)

set_temperatures!(i; Tamb=40.0, Tlim=fill(100.0, length(i.lines)))

# auto-precision, inverse of auto-correlation
i.auto_prec = 0.8

o = solve_temporal_instanton(i)

@show o.score[1]
@show o.α[1]
## basic functionality tests
@testset "Basic functionality" begin
    atol = 1e-4
    @test nd isa Dict{String, Any}
    # @test i.D0[310] ≈ 0.06412 atol=atol
    @test size(i.Y) == (118, 118)
    @test i.line_conductors[1].label == "Oriole"
    @test i.line_conductors[3].bundle == 2
    @test i.line_lengths[1] ≈ 114766.5379 atol=atol
end

## Now run tests from other files
include("temperature.jl")
include("outputchecks.jl")
