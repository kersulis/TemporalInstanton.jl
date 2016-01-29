addprocs(2)
@everywhere include("../src/TemporalInstanton.jl")
@everywhere using TemporalInstanton
d = testcase("timing")
solve_temporal_instanton(d,silent=true)
@time solve_temporal_instanton(d,silent=true);
@time solve_temporal_instanton(d,silent=true);
@time solve_temporal_instanton(d,silent=true);
