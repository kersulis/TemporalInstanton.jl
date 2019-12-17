using LinearAlgebra: eigen, I, Diagonal
using Statistics: mean

"""
    D = impedance_distance(Y)

Given an admittance matrix, return a dense distance matrix `D` such that `D[i, j]` is the impedance distance between nodes i and j.

Impedance distance is an excellent measure of electrical distance, and provides a physically meaningful planar embedding. Unfortunately, small admittance matrix values can lead to large impedances, which results in stretched graphs with tight clusters of nodes and far-flung outliers.
"""
function impedance_distance(Y::AbstractArray)
    nb = size(Y, 1)
    D = zeros(nb, nb)

    Z = inv(Array(Y))
    for i in 1:nb
        for j in 1:nb
            D[i, j] = Z[i, i] + Z[j, j] - Z[i, j] - Z[j, i]
        end
    end
    return abs.(D)
end
impedance_distance(i::InstantonInput) = impedance_distance(i.Y)

"""
Calculate pairwise power transfer distance.
`ref` is the index of the reference bus.
"""
function power_transfer_distance(ISF::Matrix, ref::Integer)
    l, n = size(ISF)
    D = zeros(n, n)
    for i in 1:n
        for j in 1:n
            if i == j
                D[i, j] = 0
            else
                P = zeros(n) # ref element will be deleted
                P[i], P[j] = 1, -1
                P[ref] = 0
                D[i, j] = sum(abs, ISF * P)
            end
        end
    end
    return D
end
function power_transfer_distance(i::InstantonInput)
    Y = copy(i.Y)
    ISF = TemporalInstanton.isf(Y, i.lines, i.ref, i.k)
    return power_transfer_distance(ISF, i.ref)
end

"""
    Xr, λ = mds(D, d)

Inputs:
* `D` is an n x n matrix such that D[i, j] is the distance from object i to object j
* `d` is the desired embedding dimension. The value used should be
capped at the number of positive eigenvalues of K

Outputs:
* `Xr` is an n x d matrix whose rows contains the relative coordinates
 of the n objects and e the n eigenvalues of K
* `λ` are the eigenvalues of K

Note: MDS is only unique up to rotation and translation, so we
enforce the following conventions on Xr:

* [ORDER] Xr[:, i] corresponds to ith largest eigenpair of K
* [EIG] Actual d used is min(d, # positive eigenvalues of K)
* [CENTER] The rows of Xr have zero centroid
* [SIGN] The largest magnitude element of Xr[:, i] is positive
"""
function mds(D::Matrix, d::Integer)
    n = size(D, 1)

    # Compute correlation matrix
    S = D.^2
    P = I - ones(n, n) / n
    K = -0.5 * P * S * P
    K = 0.5 * (K + K') # Force symmetry

    # Compute relative coordinates
    λ, U = eigen(K)
    d = min(d, count(λ .> 0)) # Apply [EIG]
    idx = sortperm(λ; rev=true) ## TODO: Fill in ??
    Xr = U[:, idx[1:d]] * Diagonal(sqrt.(λ[idx[1:d]])) # Apply [ORDER]

    # Apply [CENTER]
    Xr .-= mean(Xr; dims=1) ## TODO: Fill in ??

    # Apply [SIGN]
    Xr .*= sign.(Xr[findmax(abs.(Xr); dims=1)[2]])

    return Xr, λ
end

# """
# Use GraphLayout.jl to compute a spring layout
# for the graph corresponding to a provided
# adjacency matrix. Returns Void; useful only
# for quickly obtaining network graphs in an
# interactive environment.
#
# Optional `C` controls spring layout tightness.
# """
# function simple_graph_adj(adj,C=2.0)
#     loc_x, loc_y = layout_spring_adj(adj,C=C)
#     layout_tree
#     labels = collect(1:size(adj,1))
#
#     draw_layout_adj(adj,loc_x,loc_y,
#     labels=labels,labelc="#000000",
#     arrowlengthfrac=0.0,
#     nodefillc="#FFFFFF",
#     labelsize=2.0)
# end
