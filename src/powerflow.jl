# The following five functions are used to run power flow
# with fixed decision variable values.

using SparseArrays: sparse, blockdiag

"""
    expand_renewable_vector(x,Ridx,N,T) -> P
Starting with `x`, a `length(Ridx)*T` vector containing
generation forecast for all nodes and time steps,
insert zeros to produce an `N*T`-long vector.
"""
function expand_renewable_vector(x, Ridx, N, T)
    idx = Int64[]
    for i in 0:(T - 1)
        append!(idx, Ridx .+ i * N)
    end
    P = zeros(N * T)
    P[idx] = x
    return P
end
"""
    fixed_wind_A(T, Y, ref, k) -> A
Generate the power balance constraint A matrix from
problem dimensions, admittance matrix, and generator
participation factors. Use this function when variable
(renewable) injections are known, to run power flow
and determine angle differences throughout the network.

Assumes the admittance matrix is n-by-n.

Returns A, which is (n + 1) * T -by- (nr + n + 1) * T

* `nr` is the number of wind farms in the network
* `n` is the number of nodes in the network
* `Ridx` is a vector indicating wind farm locations
* `T` is the number of time steps
* `Y` is the admittance matrix (n-by-n)
* `ref` is the index of the angle reference bus
* `k` is the vector of generator participation factors
"""
function fixed_wind_A(T, Y, ref, k)
    function ei(n, i)
        v = zeros(n)
        v[i] = 1.0
        return v
    end

    n = size(Y, 1)

    # A has a block diagonal pattern where each
    # block is Atemp:
    Atemp = sparse([[Y       -k];
                    ei(n, ref)'  0])

    # Now we can tile the Atemp matrix to generate A:
    A = Atemp
    for t in 2:T
        A = blockdiag(A, Atemp)
    end
    return Array(A)
end

"""
    fixed_wind_b(n, T, G0, Pnet, D) -> b
Generate the vector `b` of power balance constraints.
Assumes `G0` and `D` are nT -by- 1 vectors.
"""
function fixed_wind_b(n, T, G0, Pnet, D)
    b = Float64[]
    netGen = G0 + Pnet - D

    for t in 1:T
        startidx = (t - 1) * n + 1
        stopidx = startidx + n - 1
        append!(b, netGen[startidx:stopidx])
        push!(b, 0.0)
    end
    return b
end

"""
    return_angles(fixed_x, N, T) -> angles, alpha
Given fixed injections, the number of nodes in the network,
and the number of time steps, return the set of voltage
angles and mismatches.
"""
function return_angles(fixed_x, N, T)
    angles = Vector{Float64}[]
    alpha = Float64[]
    for i in 1:T
        push!(angles, fixed_x[((N + 1) * (i - 1) + 1):((N + 1) * (i - 1) + N)])
        push!(alpha, fixed_x[(N + 1) * (i - 1) + N + 1])
    end
    return angles, alpha
end

"""
    return_angle_diffs(angles, line)
Compute the angle differences across a particular line `(from, to)` across all time steps.
"""
function return_angle_diffs(angles, line)
    angle_diffs = Float64[]
    f = line[1]
    t = line[2]
    for v in angles
        push!(angle_diffs, v[f] - v[t])
    end
    return angle_diffs
end

"""
Calculate injection shift factor matrix.
Each row corresponds to a line in the network.
Each column corresponds to a node.
Credit to Jonathon Martin for derivation.

Inputs:
* `Y`: full admittance matrix
* `lines`: vector of tuples; each tuple encodes a line as (i,j)
* `ref`: index of angle reference bus
* `k`: vector of generator participation factors
"""
function isf(
    Y::AbstractArray,
    lines::Vector{Tuple{Int64,Int64}},
    ref::Int64,
    k=[NaN]::Vector{Float64}
    )

    Y = Array(Y)
    n, l = (size(Y, 1), length(lines))

    Bflow = zeros(l, n)
    for idx in 1:l
        i, j = lines[idx]
        Bflow[idx, i] =  Y[i, j]
        Bflow[idx, j] = -Y[i, j]
    end

    if length(k) != 1
        Y[:, ref] = k
        B = Y
        Bflow[:, ref] = zeros(l)
    else
        nonref = setdiff(1:n, ref)
        B = Y[nonref, nonref]
        Bflow = Bflow[:, nonref]
    end
    return Bflow / B
end
isf(i::InstantonInput) = isf(i.Y, i.lines, i.ref, i.k)

"""
    createY(f,t,x [,r,b]) -> Y
Create an admittance matrix for AC power flow.
All inputs are real.

* `f`,`t`: vectors encoding all lines (fi,ti)
* `x`: per-unit reactance xi for all lines
* `r`: per-unit resistance ri for all lines
* `b`: per-unit susceptance bi for all lines
"""
function createY(
    f::Vector{Int64},
    t::Vector{Int64},
    x::Vector{Float64},
    r=0.0::Union{Vector{Float64},Float64},
    b=0.0::Union{Vector{Float64},Float64}
    )
    z = r .+ x * 1im
    y = 1 ./ z
    b = b * 1im
    Y = sparse(
        [f; t; t; f],
        [t; f; t; f],
        [-y; -y; (y .+ b ./ 2); (y .+ b ./ 2)])

    # for DC power flow, we typically want a matrix with real entries:
    if r == 0
        return imag(Y)
    else
        return Y
    end
end

"""
    d = return_B(network_data)

Return B matrix for use in DC power flow (P = B x theta)
"""
function return_B(network_data::Dict)
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

    return createY(f_vec, t_vec, x_vec)
end
