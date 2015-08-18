# Data loading and manipulation:
using JLD, MatpowerCases
include("tmp_inst_rts96.jl")
include("readPolishData.jl")

function load_rts96_data(;return_as_type=false)
    function return_line_conductors(bus_names,bus_voltages,from,to)
        numLines = length(from)
        node2voltage(node) = bus_voltages[find(bus_names.==node)][1]
        volt2cond(volt) = volt < 300 ? "waxwing" : "dove" # waxwing for 138, dove for 200+
        line_voltages = Array(Float64,0)

        for i in 1:numLines
            Vfrom = node2voltage(from[i])
            Vto = node2voltage(to[i])
            Vline = max(Vfrom,Vto)
            push!(line_voltages, Vline)
        end
        line_conductors = [volt2cond(volt) for volt in line_voltages]
    end

    ####### LOAD DATA ########
    psData = psDataLoad()

    # unpack psDL (boilerplate):
    (Sb,f,t,r,x,b,Y,bustype,
    Gp,Gq,Dp,Dq,Rp,Rq,
    Pmax,Pmin,Qmax,Qmin,Plim,
    Vg,Vceiling,Vfloor,
    busIdx,N,Nr,Ng,k) = unpack_psDL(psData)

    Sb = 100e6 # overwrite "100.0"

    res = r
    reac = x

    ####### LINK DATA ########
    # Static
    Ridx = find(Rp) # Vector of renewable nodes
    Y = full(Y) # Full admittance matrix (ref not removed)
    ref = 1 # Index of ref node
    k = k # Conventional generator participation factors
    lines = [(f[i],t[i]) for i in 1:length(f)]
    lines = convert(Array{Tuple{Int64,Int64}},lines)
    line_lengths = load("../data/RTS-96\ Data/line_lengths.jld", "line_lengths")

    mpc = loadcase("case96",describe=false)
    from = mpc["branch"][:,1]
    to = mpc["branch"][:,2]
    bus_voltages = mpc["bus"][:,10]
    line_conductors = return_line_conductors(busIdx,bus_voltages,from,to);

    if return_as_type
        return InstantonInputData(Ridx,Y,Gp,Dp,Rp,Sb,ref,lines,res,reac,k,line_lengths,line_conductors,NaN,NaN,NaN,NaN,[])
    else
        return Ridx,Y,Gp,Dp,Rp,Sb,ref,lines,res,reac,k,line_lengths,line_conductors
    end
end
