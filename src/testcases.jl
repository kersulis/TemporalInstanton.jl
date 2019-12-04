using PowerModels
using Random: MersenneTwister
using Memento

export testcase_118

setlevel!(getlogger(PowerModels), "error")
"""
    i = testcase_118(path_to_pglib_opf_case118_ieee.m)

Returns an instance of `InstantonInput` corresponding to a static
118-bus test case that provides a useful starting point for a few
time-varying temporal deviation scanning scenarios. Wind generation
is concentrated in one area, while demand is concentrated in another.
"""
function testcase_118(
        casepath::String=joinpath("/home/jk/jdocuments/sdp-opf-decomp/pglib-opf", "pglib_opf_case118_ieee.m")
    )
    rng = MersenneTwister(1)

    nd = PowerModels.parse_file(casepath)
    i = build_instanton_input(casepath)

    Ridx = [1; 4; 6; 8; 10; 11; 15; 18; 19; 25; 31; 32; 36; 114]
    conventional_to_renewable!(i, Ridx)

    # ensure all wind farms have non-zero forecast output
    for idx in 1:length(i.R0)
        if idx in Ridx && i.R0[idx] < 0.1
            i.R0[idx] = rand(rng, 0.6:0.05:1.0)
        end
    end

    # reduce output of huge wind farm
    i.R0[findmax(i.R0)[2]] -= 1.0

    # bump up small generators to flatten curve
    # (since deviation scanning may easily drive their output negative)
    for idx in 1:length(i.G0)
        if i.k[idx] > 0
            if i.G0[idx] < 0.3
                i.G0[idx] = rand(rng, 0.0:0.05:1.0)
            end
        end
    end

    # increase generation to match demand (keeps alpha low)
    G0 = copy(i.G0)
    G0[G0 .> 0] .+= 1.85
    G0[findmax(G0)[2]] -= 1.3
    i.G0 = G0

    D0_mod = copy(i.D0) * 2.1 # does not disturb nodes with 0 demand

    neutral_nodes = findall(D0_mod .== 0)

    # ensure there are no super tiny demands
    small_demand_idx = sortperm(D0_mod)[1:35]
    small_demand_idx = setdiff(small_demand_idx, neutral_nodes)
    D0_mod[small_demand_idx] .+= 0.2


    rural_nodes = [
    10, 9, 8, 30, 33, 19, 18, 15, 13, 16, 14, 11,
    12, 4, 5, 6, 1, 3, 117, 17, 113, 31, 115, 28,
    29, 114, 27, 32, 20, 25, 21, 34, 35, 36] |> sort
    rural_nodes = setdiff(rural_nodes, neutral_nodes)

    urban_nodes = [44; 47; 49; 51; 52; 56; 59; 60; 62; 66; 67; 47]
    urban_nodes = setdiff(urban_nodes, neutral_nodes)

    other_nodes = setdiff(1:118, union(urban_nodes, rural_nodes, neutral_nodes))

    D0_mod[urban_nodes] .+= 0.5


    # reduce demand at rural nodes
    rural_demand = D0_mod[rural_nodes]
    large_demand_idx = findall(rural_demand .> 1.0)
    rural_demand[large_demand_idx] .-= 0.6
    D0_mod[rural_nodes] = rural_demand

    large_demand_idx = findall(D0_mod .> 2.0)
    D0_mod[large_demand_idx] .*= 0.55

    # reduce one crazy high demand
    D0_mod[findmax(D0_mod)[2]] -= 2.0

    D0_mod[[80; 90; 116]] .-= 0.8

    i.D0 = D0_mod

    # re-do temperature calc using updated demand
    set_temperatures!(i; Tamb=40.0, Tlim=fill(100.0, length(i.lines)))

    # just one time step by default
    set_timing!(i, 0:600:600)

    return i
end

"""
    D0 = increase_urban_demand(D0, urban_nodes, T[, seed])

Inputs:
- `D0`: Nominal demand vector containing demand at all nodes and time steps
- `urban_nodes`: Vector of nodes at which demand should increase
- `T`: Number of time steps
- `seed`: Used to seed the random number generator before randomly generating demand increases

Increase demand at urban nodes over time, while leaving other loads relatively constant.
"""
function increase_urban_demand(D0::Vector{Float64}, urban_nodes::Vector{Int64}, T::Int64; seed::Int64=1)
    rng = MersenneTwister(seed)

    D0_final = Vector{Vector{Float64}}()
    push!(D0_final, D0)
    for i in 1:(T - 1)
        x = copy(D0_final[end])
        for idx in urban_nodes
#             x[idx] += 0.08
            x[idx] += rand(rng, 0.0:0.01:0.2)
        end
        push!(D0_final, x)
    end
    return vcat(D0_final...)
end

"""
    R0 = generate_wind_forecast(R0, wind_nodes, T[, seed])

Inputs:
- `R0`: Nominal renewable generation forecast
- `wind_nodes`: Vector of nodes where wind farms are located
- `T`: Number of time steps
- `seed`: Used to seed random number generator before generating wind forecast
"""
function generate_wind_forecast(R0::Vector{Float64}, wind_nodes::Vector{Int64}, T::Int64; seed::Int64=1)
    rng = MersenneTwister(seed)

    R0_final = Vector{Vector{Float64}}()
    push!(R0_final, R0)
    for i in 1:(T - 1)
        x = copy(R0_final[end])
        for idx in wind_nodes
            x[idx] += randn(rng) / 6.0
        end
        push!(R0_final, x)
    end
    return vcat(R0_final...)
end
