addprocs(2)
@everywhere include("../src/TemporalInstanton.jl")
@everywhere include("../src/mat2tmpinst.jl")
@everywhere using TemporalInstanton

# compile code
timing = testcase("timing")
@time o = solve_temporal_instanton(timing;silent=true);

case = "case96"
maxlines = 20
reps = 30 # used to illustrate solution time variation
line_times_96 = Vector{Vector{Float64}}()

for rep in 1:reps
    num_farms_vec = collect(2:2:72)
    line_times = Vector{Vector{Float64}}()

    for i in 1:length(num_farms_vec)
        num_farms = num_farms_vec[i]
        penetration = 0.7
        d = mat2tmpinst(case,penetration,num_farms,fill_default=true)
        d.Y = timing.Y
        # o = solve_temporal_instanton(d,silent=true, maxlines=maxlines);
        # run with this line overnight:
        o = solve_temporal_instanton(d,silent=true);
        push!(line_times,o.linetimes)
    end
    push!(line_times_96, [mean(line_times[i]) for i in 1:length(line_times)])
end

# store the data that took a half hour to obtain
using JLD
save("../data/2017-05-12-line_times_96.jld", "line_times_96", line_times_96)
