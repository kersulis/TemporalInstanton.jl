addprocs(2) # vary number of concurrent processes here
@everywhere include("../src/TemporalInstanton.jl")
@everywhere include("../src/mat2tmpinst.jl")
@everywhere using TemporalInstanton

date = Dates.format(now(),"yyyy-mm-dd")
case = "case2383wp"
note = "polish scaling"
reps = 10
fname = "polish_scaling.csv"

i = testcase("timing")
# compile everything
solve_temporal_instanton(i,silent=true)

open("../data/$fname","w") do f
    writecsv(f,["date";"note";"rep";"nl";"nb";
        "nr";"nt";"tt";"mlt";"alloc";"gc"]')
end

# polish scaling
for rep in 1:reps
    # larger than 1500 gets shaky on my laptop.
    # (tons of gc, each analysis takes over an hour)
    num_farms_vec = collect(50:100:1500)
    for i in 1:length(num_farms_vec)
        num_farms = num_farms_vec[i]
        penetration = 0.7
        d = mat2tmpinst(case,num_farms,penetration,fill_default=true)
        nl = length(d.lines)
        nb = length(d.k)
        nr = length(d.Ridx)
        nt = round(Int64,length(d.D0)/nb)

        t = @timed solve_temporal_instanton(d,maxlines = 20,silent=true)
        tt = t[2]
        mlt = mean(t[1].linetimes)
        alloc = t[3]
        gc = t[4]
        open("../data/$fname","a") do f
            writecsv(f,[date note rep nl nb nr nt tt mlt alloc gc])
        end
    end
end    
