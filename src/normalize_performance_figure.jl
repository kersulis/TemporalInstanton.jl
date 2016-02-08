addprocs(2) # vary number of concurrent processes here
@everywhere include("../src/TemporalInstanton.jl")
@everywhere include("../src/mat2tmpinst.jl")
@everywhere using TemporalInstanton

# compile everything with this run:
d = testcase("timing")
@time o = solve_temporal_instanton(d);
sort(o.score)

cnames = casenames()
gen_count = Vector{Int64}()
maxlines = 40
for case in cnames
    mpc = loadcase(case,describe=false)
    push!(gen_count,length(unique(mpc["gen"][:,1])))
end
# cases = cnames[sortperm(gen_count)[[5;6;7;8;12;13;16;17;18;19]]]
# cases = cnames[sortperm(gen_count)[[5;6;7;8;12;13;16;17;18;20;23;]]]#27;]]]#28;29;30]]]
cases = cnames[sortperm(gen_count)[[6;7;8;12;13;17;18;20;23]]]
score_output = Array(Vector{Tuple{Float64,Int64}},0)
sec_elapsed = Vector{Float64}()
bytes_alloc = Vector{Float64}()
sec_in_gc = Vector{Float64}()
wind_count = Vector{Int64}()
line_times = Vector{Vector{Float64}}()
decision_vars = Vector{Int64}()
for i in 1:length(cases)
    case = cases[i]
    println("starting $case")
#     num_farms = length(unique(loadcase(case,describe=false)["gen"][:,1]))
#     push!(wind_count,num_farms)
    penetration = 0.7 # penetration of 50 %
    d = mat2tmpinst(case,0,penetration,fill_default=true)

    timed_results = @timed solve_temporal_instanton(d,maxlines=maxlines);
    o = timed_results[1]
    push!(score_output,o.score)
    push!(sec_elapsed,timed_results[2])
    push!(bytes_alloc,timed_results[3])
    push!(sec_in_gc,timed_results[4])
    push!(line_times,o.linetimes)
    push!(decision_vars,length(o.x[1][1])*length(o.x[1]))
end

lines_processed = []
for score_vec in score_output
    push!(lines_processed, sum(score_vec.!=Inf))
end

line_times = [mean(l) for l in line_times]

using JLD
save("../data/line_times_polish.jld",
"line_times",line_times,
"decision_vars",decision_vars,
"cases",cases)
