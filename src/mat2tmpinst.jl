using MatpowerCases, HDF5, JLD

""" 20150-09-06. Modified version of mat2tmpinst.
Load (and generate) everything needed to perform temporal
instanton analysis for any network supported by MatpowerCases.

* `num_wind_farms`: how many wind farms do you want?
* `wind_penetration`: desired value of sum(G0)/sum(R0)

Note: each wind farm's output is sum(Gp)*wind_penetration/num_wind_farms
"""
function mat2tmpinst(
    name::String,
    wind_penetration::Float64,
    num_wind_farms::Int64 = 0;
    return_as_type::Bool = true,
    fill_default = false
    )
    mpc = loadcase(name,describe=false)

    bus_orig = mpc["bus"][:,1]
    bus_simple = collect(1:length(bus_orig))

    genBus = mpc["gen"][:,1]
    for i in bus_simple
        genBus[genBus.==bus_orig[i]] = bus_simple[i]
    end

    try
        Sb = mpc["baseMVA"]
    catch
        warn("Base MVA missing from mpc data.")
        Sb = 100.0
    end
    Gp_long = mpc["gen"][:,2]

    f = round(Int64,mpc["branch"][:,1]) # "from bus" ...
    t = round(Int64,mpc["branch"][:,2]) # ... "to bus"
    for i in bus_simple
        f[f.==bus_orig[i]] = bus_simple[i]
        t[t.==bus_orig[i]] = bus_simple[i]
    end
    r = mpc["branch"][:,3]              # resistance, pu
    x = mpc["branch"][:,4]              # reactance, pu
    b = mpc["branch"][:,5]              # susceptance, pu

    Y = createY(f,t,x)

    Gp = zeros(length(bus_simple))
    for i in bus_simple
        Gp[convert(Int64,i)] = sum(Gp_long[find(genBus.==i)])/Sb
    end

    Dp = mpc["bus"][:,3]./Sb

    # convert generators into wind farms:
    if num_wind_farms == 0
        num_wind_farms = length(find(Gp))
    end
    windsize = sum(Gp)*wind_penetration/num_wind_farms
#     Ridx = rand(1:length(Gp),num_wind_farms)
    # Ridx = rand(bus_simple,num_wind_farms)#sortperm(Dp)[1:num_wind_farms]
    Ridx = randperm(length(bus_simple))[1:num_wind_farms]
    Rp = zeros(length(Gp))
    Rp[Ridx] = windsize

    Sb = Sb*1e6 # convert from MW to W
    ref = 1
    lines = [(f[i],t[i]) for i in 1:length(f)]
    lines = convert(Array{Tuple{Int64,Int64},1},lines)

#     res = r
#     reac = x

    # Allow each generator to participate equally in droop response.
    k = Float64[]
    for i = 1:length(Gp)
        if Gp[i] != 0
            push!(k, 1/length(find(Gp)))
        else
            push!(k,0)
        end
    end

    # use RTS-96 line lengths to generate similar line lengths
    line_lengths = load("../data/7k_line_lengths.jld","line_lengths")[1:length(lines)]

    # temporary (re-use rts-96 line conductor parameters)
    line_conductors = fill("waxwing",length(line_lengths))

    if return_as_type
        d = InstantonInputData(Ridx,Y,Gp,Dp,Rp,Sb,ref,lines,r,x,k,
        line_lengths,line_conductors,NaN,NaN,NaN,0.0:0.0,Array{Float64,2}())

        if fill_default
            d.Tamb = 35. # C
            d.T0 = 60. #46. # initial line steady-state temp
            d.time_values = 0:30:300 # five minutes in 30-sec steps
            d.int_length = 300. # seconds = 5 min

            # 6 time steps
            Gp,Dp,Rp = (d.G0, d.D0, d.R0)
            d.G0 = [Gp;Gp;Gp;Gp;Gp;Gp]
            d.D0 = [Dp;Dp;Dp;Dp;Dp;Dp]
            d.R0 = [Rp;Rp;Rp;Rp;Rp;Rp]
        end
        return d
    else
        return Ridx,Y,Gp,Dp,Rp,Sb,ref,lines,r,x,k,line_lengths,line_conductors
    end
end
