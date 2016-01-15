addprocs(2) # vary number of concurrent processes here
@everywhere include("../src/TemporalInstanton.jl")
@everywhere include("../src/mat2tmpinst.jl")
@everywhere using TemporalInstanton

date = Dates.format(now(),"yyyy-mm-dd")
case = "case96"
note = "rts-96 scaling"
reps = 10

i = testcase("timing")
# compile everything
solve_temporal_instanton(i,silent=true)

open("../data/rts96_scaling.csv","w") do f
    # header
    writecsv(f,["date";"note";"rep";"nl";"nb";
        "nr";"nt";"tt";"mlt";"alloc";"gc"]')
end

# rts-96 scaling
for rep in 1:reps
    num_farms_vec = collect(2:2:72)
    for i in 1:length(num_farms_vec)
        num_farms = num_farms_vec[i]
        penetration = 0.7
        d = mat2tmpinst(case,num_farms,penetration,fill_default=true)
        nl = length(d.lines)
        nb = length(d.k)
        nr = length(d.Ridx)
        nt = round(Int64,length(d.D0)/nb)

        t = @timed solve_temporal_instanton(d,silent=true)
        tt = t[2]
        mlt = mean(t[1].linetimes)
        alloc = t[3]
        gc = t[4]
        open("../data/rts96_scaling.csv","a") do f
            writecsv(f,[date note rep nl nb nr nt tt mlt alloc gc])
        end
    end
end
