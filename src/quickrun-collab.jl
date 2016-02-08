include("TemporalInstanton.jl")
using TemporalInstanton, MAT
d = matread("../collaboration/D20160203-1.mat")["System"];
ps = dict2type(d);
d = ps2instantoninput(ps);
solve_temporal_instanton(d)
