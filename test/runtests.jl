using TemporalInstanton
# using PowerModels

# include("../src/TemporalInstanton.jl")

fpath = joinpath(@__DIR__, "data", "caseRTS96.mat")

TemporalInstanton.load_power_system(fpath)

using Test
