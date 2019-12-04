# Output checks: is solution feasible?

## VULNERABILITY
# show wind injections at instanton scenario; flag any negative values

# issue warning if other line temperatures increase

# Distinguish between lines that are radial and transmission lines

## REDISPATCH
# check generator output limits
# use network data dictionary and run DC pf on deviation scanning output scenario. Compare total generation at each bus with the sum of generator limits at that bus.

export check_decision_variables, check_generator_outputs

"""
    genlim = generator_limits(nd)

Return generator limits as a dict, with bus ids as keys and (pmin, pmax) tuples as values.
"""
function generator_limits(nd::Dict{String, Any})
    nb = nd["bus"] |> length

    genbus = Dict{String, Int64}()
    pminmax = Dict{String, Tuple}()
    for (gid, gd) in nd["gen"]
        genbus[gid] = gd["gen_bus"]
        pminmax[gid] = (gd["pmin"], gd["pmax"])
    end

    genlim = Vector{Tuple}(undef, nb)
    for bus in values(genbus)
        pmin = 0.0
        pmax = 0.0
        for gid in keys(genbus)
            if genbus[gid] == bus
                pmin += pminmax[gid][1]
                pmax += pminmax[gid][2]
            end
        end
        genlim[bus] = (pmin, pmax)
    end
    return genlim
end

"""
    P = return_injections(network_data, o, event_idx, time_interval)

Inputs:
- `network_data`: The dictionary returned by PowerModels.parse_file
- `o`: An instance of InstantonOutput
- `event_idx`: Integer representing which deviation pattern to use
- `time_interval`: Which time interval to run DC power flow for

Given the network data dictionary and an instance of InstantonOutput, along with an event idx (e.g. 1 corresponds to the most likely deviation pattern) and a time index, perform DC power flow and return the injections.
"""
function return_injections(network_data::Dict, o::InstantonOutput, event_idx::Integer=1)
    B = TemporalInstanton.return_B(network_data)
    P = Vector{Vector{Float64}}()
    for theta in o.θ[event_idx]
        push!(P, B * theta)
    end
    return P
end

"""
    D = return_demand(i)

Given an instance of `InstantonInput`, return a vector of T vectors, where each vector provides the active power demand for all buses at the corresponding time interval (T total intervals).
"""
function return_demand(i::InstantonInput)
    D0 = i.D0
    D = Vector{Vector{Float64}}()
    T = length(i.time_values) - 1
    n = length(D0) / T |> Integer
    for t in 1:T
        push!(D, D0[(n * (t - 1) + 1):(n * t)])
    end
    return D
end

"""
    G = return_generation(i)

Given an instance of `InstantonInput`, return a vector of T vectors, where each vector provides the active power generation at all buses for the corresponding time interval (T total intervals).
"""
function return_generation(i::InstantonInput, o::InstantonOutput, event_idx::Integer=1)
    G0, k, A = i.G0, i.k, o.α[event_idx]

    G = Vector{Vector{Float64}}()
    T = length(i.time_values) - 1
    n = length(G0) / T |> Integer

    for (t, α) in zip(1:T, A)
        G0_t = G0[(n * (t - 1) + 1):(n * t)]
        push!(G, G0_t .+ α * k)
    end
    return G
end

"""
    warnings = check_generator_outputs(nd, i, o, event_idx)

Check generation corresponding to deviation scanning event `event_idx`. Prints warning messages if any generator output is below its lower limit or above its upper limit.
"""
function check_generator_outputs(nd::Dict{String, Any}, i::InstantonInput, o::InstantonOutput, event_idx::Integer=1)
    gen_limits = generator_limits(nd)
    generation = return_generation(i, o, event_idx)

    nb = length(gen_limits)
    T = length(generation)

    warnings = Vector{Dict}()
    for t in 1:T
        gen_t = generation[t]

        for i in 1:nb
            ! isassigned(gen_limits, i) && continue

            lim_i, G_i = gen_limits[i], gen_t[i]
            l, u = lim_i

            warn_string = "Time interval $(t), bus $(i): Generator output $(round(G_i; digits=2)) is "
            if G_i < l
                @warn warn_string * "below lower limit of $(l)."
                push!(
                    warnings, Dict(
                        "time_interval" => t,
                        "bus" => i,
                        "output" => G_i,
                        "l" => l
                        )
                    )
            elseif G_i > u
                @warn warn_string * "above upper limit of $(u)."
                push!(
                    warnings, Dict(
                        "time_interval" => t,
                        "bus" => i,
                        "output" => G_i,
                        "u" => u
                        )
                    )
            end
        end
    end
    return warnings
end

"""
    R = return_decision_variables(i, o[, event_idx])

Given instances of `InstantonInput` and `InstantonOutput`, return a vector of T vectors, where each vector provides the decision variable values for all buses at the corresponding time interval (T total intervals) (0 if there is no decision variable at a bus).
"""
function return_decicion_variables(i::InstantonInput, o::InstantonOutput, event_idx::Integer=1)
    R0, Ridx = i.R0, i.Ridx
    deviations = o.x[event_idx]
    R = Vector{Vector{Float64}}()
    T = length(i.time_values) - 1
    n = length(R0) / T |> Integer

    for t in 1:T
        R_t = R0[(n * (t - 1) + 1):(n * t)]
        R_t[Ridx] .+= deviations[t]
        push!(R,  R_t)
    end
    return R
end

"""
    warnings = check_decision_variables(nd, i, o, event_idx)

Check generation corresponding to deviation scanning event `event_idx`. Prints warning messages if any generator output is below its lower limit or above its upper limit.
"""
function check_decision_variables(i::InstantonInput, o::InstantonOutput, event_idx::Integer=1)
    dec_vars = return_decicion_variables(i, o, event_idx)

    Ridx = i.Ridx
    T = length(i.time_values) - 1

    warnings = Vector{Dict}()
    for t in 1:T
        for i in Ridx
            inj = dec_vars[t][i]
            if inj < 0
                @warn "Time interval $(t), bus $(i): Decision variable generator output $(round(inj; digits=2)) is negative."
                push!(
                    warnings, Dict(
                        "time_interval" => t,
                        "bus" => i,
                        "output" => inj
                        )
                    )
            end
        end
    end
    return warnings
end
