# Data loading and manipulation:
using JLD, MatpowerCases
include("tmp_inst_rts96.jl")

type InstantonInputData
    Ridx::Vector{Int64}
    Y::SparseMatrixCSC{Float64,Int64}
    G0::Vector{Float64}
    D0::Vector{Float64}
    R0::Vector{Float64}
    Sb::Float64
    ref::Int64
    lines::Array{Tuple{Int64,Int64},1}
    res::Vector{Float64}
    reac::Vector{Float64}
    k::Vector{Float64}
    line_lengths::Vector{Float64}
    line_conductors::Vector{ASCIIString}
    Tamb::Float64
    T0::Float64
    int_length::Float64
    time_values
    corr::Array{Float64,2}
end

type InstantonOutputData
    score::Vector{Float64}
    x::Vector{Vector{Vector{Float64}}}
    θ::Vector{Vector{Vector{Float64}}}
    α::Vector{Vector{Float64}}
    diffs::Vector{Vector{Float64}}
    xopt::Vector{Vector{Float64}}
end

function load_rts96_data(;return_as_type=false)
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
    Y = Y # Full admittance matrix (ref not removed)
    ref = 1 # Index of ref node
    k = k # Conventional generator participation factors
    lines = [(f[i],t[i]) for i in 1:length(f)]
    lines = convert(Array{Tuple{Int64,Int64},1},lines)
    line_lengths = load("../data/RTS-96\ Data/line_lengths.jld", "line_lengths")

    mpc = loadcase("case96",describe=false)
    from = mpc["branch"][:,1]
    to = mpc["branch"][:,2]
    bus_voltages = mpc["bus"][:,10]
    line_conductors = return_line_conductors(busIdx,bus_voltages,from,to)

    if return_as_type
        return InstantonInputData(
            Ridx,
            Y,
            Gp,
            Dp,
            Rp,
            Sb,
            ref,
            lines,
            res,
            reac,
            k,
            line_lengths,
            line_conductors,
            NaN,
            NaN,
            NaN,
            NaN,
            Array{Float64,2}())
    else
        return Ridx,Y,Gp,Dp,Rp,Sb,ref,lines,res,reac,k,line_lengths,line_conductors
    end
end

function mat2tmpinst(name)
    mpc = loadcase(name,describe=false)
    bus_i = mpc["bus"][:,1]
    genBus = mpc["gen"][:,1]
    Sb = mpc["baseMVA"]
    Gp_long = mpc["gen"][:,2]

    f = round(Int64,mpc["branch"][:,1]) # "from bus" ...
    t = round(Int64,mpc["branch"][:,2]) # ... "to bus"
    r = mpc["branch"][:,3]              # resistance, pu
    x = mpc["branch"][:,4]              # reactance, pu
    b = mpc["branch"][:,5]              # susceptance, pu

    Y = createY(f,t,x)

    Gp = zeros(length(bus_i))
    for i in bus_i
        Gp[convert(Int64,i)] = sum(Gp_long[find(genBus.==i)])/Sb
    end

    Dp = mpc["bus"][:,3]./Sb

    # convert half the generators into wind farms:
    Rp = zeros(length(Gp))
    for i in 1:length(Gp)
        if Gp[i] < mean(Gp)
            Rp[i] = Gp[i]
            Gp[i] = 0
        end
    end

    Ridx = find(Rp)

    Sb = Sb*1e6 # convert from MW to W

    ref = 1

    lines = [(f[i],t[i]) for i in 1:length(f)]
    lines = convert(Array{Tuple{Int64,Int64},1},lines)

    res = r
    reac = x

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
    line_lengths = load("../data/polish_line_lengths.jld",
        "line_lengths")[1:length(lines)]

    # temporary (re-use rts-96 line conductor parameters)
    line_conductors = fill("waxwing",length(line_lengths))
    convert(Array{ASCIIString},line_conductors)

    return  Ridx, Y,
            Gp, Dp, Rp,
            Sb, ref, lines,
            res, reac, k,
            line_lengths, line_conductors
end

""" Create an admittance matrix for AC power flow.
"""
function createY(
    f::Vector{Int64},
    t::Vector{Int64},
    r::Vector{Float64},
    x::Vector{Float64},
    b::Vector{Float64}
    )
    G = 1./r
    G[G.==Inf] = 0
    B = 1./x
    y = complex(G,B)
    return sparse([f; t; t; f],[t; f; t; f],[-y; -y; y + b./2; y + b./2])
end

""" Create an admittance matrix for DC power flow.
"""
function createY(
    f::Vector{Int64},
    t::Vector{Int64},
    x::Vector{Float64}
    )
    y = 1./x
    return sparse([f; t; t; f],[t; f; t; f],[-y; -y; y; y])
end

""" Use bus voltage level to determine appropriate conductor type. TODO: replace with Jon's conductor interpolation code.
"""
function return_line_conductors(bus_names,bus_voltages,from,to)
    numLines = length(from)
    node2voltage(node) = bus_voltages[find(bus_names.==node)][1]
    volt2cond(volt) = volt < 300 ? "waxwing" : "dove"
    line_voltages = Array(Float64,0)

    for i in 1:numLines
        Vfrom = node2voltage(from[i])
        Vto = node2voltage(to[i])
        Vline = max(Vfrom,Vto)
        push!(line_voltages, Vline)
    end
    line_conductors = [volt2cond(volt) for volt in line_voltages]
    convert(Array{ASCIIString},line_conductors)
    return line_conductors
end
