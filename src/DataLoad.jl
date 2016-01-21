using JLD, MAT, MatpowerCases

"""
Contains fields for all data required to perform
temporal instanton analysis. Rather than passing a
list of arguments to `solve_temporal_instanton`,
one can simply pass an instance of this type.
"""
type InstantonInputData
    "Bus indices having variable (renewable) generation"
    Ridx::Vector{Int64}
    "Admittance matrix, entries in pu"
    Y::SparseMatrixCSC{Float64,Int64}
    "Conventional (fixed) generation injections at all nodes and time steps; pu"
    G0::Vector{Float64}
    "Demands at all nodes and time steps; pu"
    D0::Vector{Float64}
    "Wind injection forecasts for all nodes and time steps; pu"
    R0::Vector{Float64}
    "System base MVA, in VA"
    Sb::Float64
    "Index of system angle reference bus (need not be slack bus)"
    ref::Int64
    "Each line is represented as a tuple of the form (from,to)"
    lines::Array{Tuple{Int64,Int64},1}
    "Resistance of each line in `lines`; pu"
    res::Vector{Float64}
    "Reactance of each line in `lines`; pu"
    reac::Vector{Float64}
    "Conventional generator participation factors; must sum to 1.0"
    k::Vector{Float64}
    "Lengths of all lines in `lines`; meters"
    line_lengths::Vector{Float64}
    "Conductor type for all lines in `lines`"
    line_conductors::Vector{ASCIIString}
    "System ambient temperature in degrees C"
    Tamb::Float64
    "Initial transmission line temperature in degrees C"
    T0::Float64
    "Length of optimization horizon in seconds"
    int_length::Float64
    "Times at which generator dispatch, demand, and wind forecast are updated"
    time_values::FloatRange{Float64}
    "Variable (renewable) generation correlation matrix"
    corr::Array{Float64,2}
end

"""
Contains all data coming out of temporal instanton analysis.
`process_instanton_results` returns an instance of this type
by default to keep the workspace clean.
"""
type InstantonOutputData
    "Score vector; each entry corresponds to entry of `lines`"
    score::Vector{Tuple{Float64,Int64}}
    "Variable (renewable) generation forecast deviation vectors"
    x::Vector{Vector{Vector{Float64}}}
    "Voltage angles at all nodes and time steps"
    θ::Vector{Vector{Vector{Float64}}}
    "Mismatch between generation dispatch + forecast and demand (amount that must be taken up by droop response)"
    α::Vector{Vector{Float64}}
    "Angle difference between endpoints of instanton candidate line across all time steps"
    diffs::Vector{Vector{Float64}}
    "Vector containing all solution information in raw form; used for debugging"
    xopt::Vector{Vector{Float64}}
    "Time taken to solve each line"
    linetimes::Vector{Float64}
end

"""
Loads RTS-96 data from Jenny's ARPA-E data, supplements with
conductor type and line length information, and returns as
an instance of `InstantonInputData`.
"""
function load_rts96_data(; return_as_type::Bool = true)
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
            0.0:0.0,
            Array{Float64,2}()
            )
        else
            return Ridx,Y,Gp,Dp,Rp,Sb,ref,
                lines,r,x,k,
                line_lengths,line_conductors
    end
end

"""
    mat2tmpinst(name) -> d

Loads (and generates) everything needed to perform
temporal instanton analysis for any network supported by
[MatpowerCases.jl](https://github.com/kersulis/MatpowerCases.jl).

Output `d` is an instance of `InstantonOutputData`.
"""
function mat2tmpinst(name::ASCIIString)
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
    line_lengths = load("../data/polish_line_lengths.jld","line_lengths")[1:length(lines)]

    # temporary (re-use rts-96 line conductor parameters)
    line_conductors = fill("waxwing",length(line_lengths))

    return InstantonInputData(Ridx,Y,Gp,Dp,Rp,Sb,ref,
            lines,res,reac,k,
            line_lengths,line_conductors,
            NaN,NaN,NaN,NaN:NaN:NaN,Array{Float64,2}())
end

"""
Use bus voltage level to determine appropriate conductor type.
TODO: replace with Jon's conductor interpolation code.
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

include("mat2tmpinst.jl")
"Test case used for timing analysis"
function testcase(name::ASCIIString)
    if name == "timing"
        d = load_rts96_data(return_as_type=true);
        # Thermal model parameters:
        d.Tamb = 35. # C
        d.T0 = 60. #46. # initial line steady-state temp

        d.time_values = 0:30:300 # five minutes in 30-sec steps
        d.int_length = 300. # seconds = 5 min
        Gp,Dp,Rp = d.G0,d.D0,d.R0
        d.G0 = [0.7*Gp;0.7*Gp;0.7*Gp;0.7*Gp;0.7*Gp;0.7*Gp]
        d.D0 = [0.9*Dp;0.9*Dp;0.9*Dp;0.9*Dp;0.9*Dp;0.9*Dp]
        d.R0 = [Rp;1.1*Rp;1.2*Rp;1.3*Rp;1.4*Rp;1.5*Rp]
        return d
    elseif name == "polish200"
        d = load("../data/polish200.jld","d")
    else
        error("No match.")
    end
end
