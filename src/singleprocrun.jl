include("../src/TemporalInstanton.jl")
using TemporalInstanton
d = testcase("timing")
solve_temporal_instanton(d,silent=true)
@time solve_temporal_instanton(d,silent=true);
@time solve_temporal_instanton(d,silent=true);
@time solve_temporal_instanton(d,silent=true);
