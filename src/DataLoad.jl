using JLD, MAT, MatpowerCases

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
    time_values::FloatRange{Float64}
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

function load_rts96_data(; return_as_type::Bool = false)
    path="../src/caseRTS96.mat" # Assumes current dir is nbs
    caseRTS96 = matread(path) # import MATLAB workspace

    # Connect relevant MATLAB variables to Julia
    bus_i = caseRTS96["bus_i"][:,1] # vector of unique bus indices (73)
    bus = caseRTS96["bus"][:,1] # vector of generator bus indices (99)
    Sb = caseRTS96["Sb"]::Float64 # base MVA
    bustype = round(Int64,caseRTS96["type"][:,1])::Vector{Int64}

    Gp_long = caseRTS96["Pg"] # conventional active power output, divide by Sb later
    Gq_long = caseRTS96["Qg"] # conventional reactive power output, divide by Sb later

    Rp = (caseRTS96["Wind"]./Sb)[:,1] # renewable active generation, has zero where no farm exists

    Dp = caseRTS96["Pd"]./Sb + Rp # I previously subtracted wind from active demand
    Dp = Dp[:,1]
    f = round(Int64,caseRTS96["fbus"][:,1])::Vector{Int64} # "from bus" ...
    t = round(Int64,caseRTS96["tbus"][:,1])::Vector{Int64} # ... "to bus"
    r = caseRTS96["r"][:,1]::Vector{Float64} # resistance, pu
    x = caseRTS96["x"][:,1]::Vector{Float64} # reactance, pu

    # map bus numbers to 1:73
    for i = 1:length(f)
        if f[i] < 201
            f[i] -= 100
        elseif f[i] < 301
            f[i] -= 176
        else
            f[i] -= 252
        end
    end
    for i = 1:length(t)
        if t[i] < 201
            t[i] -= 100
        elseif t[i] < 301
            t[i] -= 176
        else
            t[i] -= 252
        end
    end

    busidx = Int64[1]
    for i = 2:length(bus)
        if bus[i] < 201
            push!(busidx, bus[i]-100)
        elseif bus[i] < 301
            push!(busidx, bus[i]-176)
        else
            push!(busidx, bus[i]-252)
        end
    end
    busidx = convert(Vector{Int64},busidx)

    # area 1: buses 1-24
    # area 2: buses 25-48
    # area 3: buses 49-73

    Gp = zeros(length(bus_i))

    for i in unique(busidx)
        Gp[i] = sum(Gp_long[find(busidx .== i)])/Sb
    end
    # Now Gp and Gq reflect active and reactive generation at buses 1:73 consecutively.

    Ridx = find(Rp)
    Y = createY(f,t,x)::SparseMatrixCSC{Float64,Int64}

    # Allow each generator to participate equally in droop response.
    # Note: this only applies to analysis types with droop response!
    k = Vector{Float64}()
    for i = 1:length(Gp)
        if Gp[i] != 0.0
            push!(k, 1/length(find(Gp)))
        else
            push!(k,0.0)
        end
    end

    ref = 1 # index of reference node
    lines = [(f[i],t[i]) for i in 1:length(f)]
    lines = convert(Array{Tuple{Int64,Int64},1},lines)
    line_lengths = load("../data/RTS-96\ Data/line_lengths.jld", "line_lengths")
    Sb = 100e6 # overwrite "100.0"

    mpc = loadcase("case96",describe=false)
    from = convert(Vector{Int64},mpc["branch"][:,1])
    to = convert(Vector{Int64},mpc["branch"][:,2])
    bus_voltages = mpc["bus"][:,10]
    line_conductors = return_line_conductors(round(Int64,bus_i),bus_voltages,from,to)

    if return_as_type
        return  InstantonInputData(
            Ridx,
            Y,
            Gp,
            Dp,
            Rp,
            Sb,
            ref,
            lines,
            r,
            x,
            k,
            line_lengths,
            line_conductors,
            NaN,
            NaN,
            NaN,
            NaN,
            Array{Float64,2}()
            )
        else
            return Ridx,Y,Gp,Dp,Rp,Sb,ref,
                lines,r,x,k,
                line_lengths,line_conductors
    end
end

""" Load (and generate) everything needed to perform
temporal instanton analysis for any network
supported by MatpowerCases
"""
function mat2tmpinst(name::ASCIIString; return_as_type::Bool = false)
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
    Rp = zeros(length(Gp))
    for i in 1:length(Gp)
        if Gp[i] < mean(Gp[find(Gp)])
            Rp[i] = Gp[i]
            Gp[i] = 0
        end
    end

    Ridx = find(Rp)

    Sb = Sb*1e6 # convert from MW to W

    ref = 1

    lines = [(f[i],t[i]) for i in 1:length(f)]

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
    line_lengths = load("../data/polish_line_lengths.jld","line_lengths")[1:length(lines)]

    # temporary (re-use rts-96 line conductor parameters)
    line_conductors = fill("waxwing",length(line_lengths))

    if return_as_type
        return InstantonInputData(Ridx,Y,Gp,Dp,Rp,Sb,ref,
            lines,res,reac,k,
            line_lengths,line_conductors,
            NaN,NaN,NaN,NaN,Array{Float64,2}())
    else
        return Ridx,Y,Gp,Dp,Rp,Sb,ref,
            lines,res,reac,k,
            line_lengths,line_conductors
    end
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
function return_line_conductors(
    bus_names::Vector{Int64},
    bus_voltages::Vector{Float64},
    from::Vector{Int64},
    to::Vector{Int64}
    )
    numLines = length(from)
    function node2voltage(node,bus_names,bus_voltages)
        bus_voltages[find(bus_names.==node)][1]
    end
    volt2cond(volt) = volt < 300 ? "waxwing" : "dove"
    line_voltages = Array(Float64,0)

    for i in 1:numLines
        Vfrom = node2voltage(from[i],bus_names,bus_voltages)
        Vto = node2voltage(to[i],bus_names,bus_voltages)
        Vline = max(Vfrom,Vto)
        push!(line_voltages, Vline)
    end
    line_conductors = [volt2cond(volt) for volt in line_voltages]
    convert(Array{ASCIIString},line_conductors)
    return line_conductors
end
