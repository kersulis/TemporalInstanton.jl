using JLD
using SparseArrays: SparseMatrixCSC
using PowerModels: parse_file
using Statistics: mean
using LineThermalModel

"""
Contains fields for all data required to perform
temporal instanton analysis. Rather than passing a
list of arguments to `solve_temporal_instanton`,
one can simply pass an instance of this type.
"""
mutable struct InstantonInput
    "Bus indices having variable (renewable) generation"
    Ridx::Vector{Int64}
    "Admittance matrix, entries in pu"
    Y::SparseMatrixCSC{Float64, Int64}
    "Conventional (fixed) generation injections at all nodes and time steps; pu"
    G0::Vector{Float64}
    "Demands at all nodes and time steps; pu"
    D0::Vector{Float64}
    "Wind injection forecasts for all nodes and time steps; pu"
    R0::Vector{Float64}
    "System base in VA"
    Sb::Float64
    "Index of system angle reference bus (need not be slack bus)"
    ref::Int64
    "Each line is represented as a tuple of the form (from,to)"
    lines::Array{Tuple{Int64, Int64}, 1}
    "Resistance of each line in `lines`; pu"
    res::Vector{Float64}
    "Reactance of each line in `lines`; pu"
    reac::Vector{Float64}
    "Conventional generator participation factors; must sum to 1.0"
    k::Vector{Float64}
    "Lengths of all lines in `lines`; meters"
    line_lengths::Vector{Float64}
    "Conductor type for all lines in `lines`"
    line_conductors::Vector{ACSRSpecsMetric}
    "Current limits for all lines, in Amps"
    current_limits::Vector{Float64}
    "System ambient temperature in degrees C"
    Tamb::Float64
    "Initial transmission line temperatures in degrees C"
    T0::Vector{Float64}
    "Limit temperatures for transmission lines in degrees C"
    Tlim::Vector{Float64}
    "Times at which generator dispatch, demand, and wind forecast are updated"
    time_values::StepRangeLen{Int64}
    "Precision coefficient for temporal auto-correlation"
    auto_prec::Float64
end

"""
Contains all data coming out of temporal instanton analysis.
`process_instanton_results` returns an instance of this type
by default to keep the workspace clean.
"""
struct InstantonOutput
    "Score vector; each entry corresponds to entry of `lines`"
    score::Vector{Tuple{Float64, Int64}}
    "Indices of lines for which instanton analysis was performed"
    analytic_lines::Vector{Int64}
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

"Test case used for timing analysis"
function testcase(name::String)
    if name == "timing"
        d = load_rts96_data(return_as_type=true)
        # Thermal model parameters:
        d.Tamb = 35. # C
        d.T0 = 60. #46. # initial line steady-state temp

        d.time_values = 0:30:300 # five minutes in 30-sec steps
        d.int_length = 300. # seconds = 5 min
        Gp, Dp, Rp = d.G0, d.D0, d.R0
        d.G0 = [0.7 * Gp; 0.7 * Gp; 0.7 * Gp; 0.7 * Gp; 0.7 * Gp; 0.7 * Gp]
        d.D0 = [0.9 * Dp; 0.9 * Dp; 0.9 * Dp; 0.9 * Dp; 0.9 * Dp; 0.9 * Dp]
        d.R0 = [Rp; 1.1 * Rp; 1.2 * Rp; 1.3 * Rp; 1.4 * Rp; 1.5 * Rp]
        return d
    elseif name == "polish200"
        d = load("../data/polish200.jld","d")
    else
        error("No match.")
    end
end

"""
    d = build_instanton_input(fpath)

Given a path to a .m file, return an instance of InstantonInput.
"""
function build_instanton_input(fpath::String)
    network_data = parse_file(fpath)
    @assert network_data["per_unit"]

    nb = length(network_data["bus"])

    bus_orig = sort([b["bus_i"] for b in values(network_data["bus"])])
    bus_simple = collect(1:length(bus_orig))
    bus_orig2simple = Dict(zip(bus_orig, bus_simple))

    Sb = network_data["baseMVA"] * 1e6

    # branch key => orig_bus
    f = Dict((k => br["f_bus"] for (k, br) in network_data["branch"]))
    t = Dict((k => br["t_bus"] for (k, br) in network_data["branch"]))

    r = Dict((k => br["br_r"] for (k, br) in network_data["branch"]))
    x = Dict((k => br["br_x"] for (k, br) in network_data["branch"]))
    b = Dict((k => (br["b_fr"] + br["b_to"]) for (k, br) in network_data["branch"]))

    current_limits = return_current_limits(network_data)
    line_lengths = estimate_length(network_data)
    line_conductors = acsr_interpolation(network_data)

    f_vec = Int64[]
    t_vec = Int64[]
    r_vec = Float64[]
    x_vec = Float64[]
    b_vec = Float64[]
    line_lengths_vec = Float64[]
    line_conductors_vec = ACSRSpecsMetric[]
    current_limits_vec = Float64[]

    # f, t, r, x, b all aligned
    for k in keys(f)
        push!(f_vec, bus_orig2simple[f[k]])
        push!(t_vec, bus_orig2simple[t[k]])
        push!(r_vec, r[k])
        push!(x_vec, x[k])
        push!(b_vec, b[k])
        push!(line_lengths_vec, line_lengths[k])
        push!(line_conductors_vec, line_conductors[k])
        push!(current_limits_vec, current_limits[k])
    end

    lines = collect(zip(f_vec, t_vec))

    # Y = createY(f, t, x)
    Y = real(createY(f_vec, t_vec, x_vec, r_vec, b_vec))

    # aggregate generation by bus
    Gp = Dict((b["bus_i"] => 0.0 for b in values(network_data["bus"])))
    for g in values(network_data["gen"])
        Gp[g["gen_bus"]] += g["pg"]
    end

    # convert gen dict to vec
    Gp_vec = zeros(nb)
    for (k, v) in Gp
        Gp_vec[bus_orig2simple[k]] = v
    end

    ng = length(findall(collect(values(Gp)) .> 0))
    participation = zeros(nb)
    participation[findall(Gp_vec .> 0)] .= 1 / ng

    Dp = zeros(nb)
    for l in values(network_data["load"])
        Dp[bus_orig2simple[l["load_bus"]]] = l["pd"]
    end

    ref = 1

    ## TODO: does any of the above code break if these are strings?

    Ridx = Int64[]

    Tamb = NaN
    T0 = Float64[]
    Tlim = Float64[]
    time_values = 0:1:0
    auto_prec = 0.0

    return InstantonInput(
        Ridx, Y, Gp_vec, Dp, Float64[], Sb, ref, lines, r_vec, x_vec,
        participation, line_lengths_vec, line_conductors_vec, current_limits_vec,
        Tamb, T0, Tlim, time_values, auto_prec
    )
end

function conventional_to_renewable!(d::InstantonInput, Ridx::Vector{Int64}=Int64[])
    Gp = d.G0

    # if idx passed in, convert those generators to renewable

    # default behavior: convert generators with below-average output
    # to renewable
    if isempty(Ridx)
        indices = findall(Gp .> 0)
        mean_gen = mean(Gp[indices])

        Ridx = findall(0 .< Gp .< mean_gen)
    end

    Rp = zeros(length(Gp))
    Rp[Ridx] = Gp[Ridx]
    for idx in Ridx
        Gp[idx] = 0.0
    end
    # Gp[Ridx] = 0.0

    d.G0 = Gp
    d.Ridx = Ridx
    d.R0 = Rp
    return d
end

function set_temperatures!(inst_input::InstantonInput; Tamb::Float64=NaN, T0::Vector{Float64}=Float64[], Tlim::Vector{Float64}=Float64[])

    if isnan(Tamb)
        @warn "Tamb empty; filling with 40.0"
        inst_input.Tamb = 40.0
    else
        inst_input.Tamb = Tamb
    end

    if isempty(Tlim)
        @warn "Tlim empty; filling with 100.0"
        inst_input.Tlim = fill(100.0, length(inst_input.lines))
    else
        inst_input.Tlim = Tlim
    end

    if isempty(T0)
        inst_input.T0 = steady_state_temps(inst_input)
    else
        inst_input.T0 = T0
    end
    return inst_input
end

function set_timing!(inst_input::InstantonInput, time_values::StepRange{Int64})
    inst_input.time_values = time_values
    return inst_input
end

function set_injections!(
    inst_input::InstantonInput;
    G0::Vector{Float64}=inst_input.G0,
    D0::Vector{Float64}=inst_input.D0,
    R0::Vector{Float64}=inst_input.R0
    )

    inst_input.G0 = G0
    inst_input.D0 = D0
    inst_input.R0 = R0
    return inst_input
end
