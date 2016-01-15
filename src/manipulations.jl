"""
Partition A into A1, A2, A3 where:

* A1 corresponds to wind
* A2 corresponds to angles + mismatch
* A3 corresponds to angle difference vars

Used to find x_star, the min-norm solution to
Ax=b such that `x_star[idx3] = 0`.
"""
function partition_A(
    A::SparseMatrixCSC{Float64,Int64},
    Qobj::SparseMatrixCSC{Float64,Int64},
    T::Int64
    )
    n = size(A,2)
    idx1 = find(diag(Qobj))
    idx2 = setdiff(1:n-T,idx1)
    idx3 = n-T+1:n

    (A1,A2) = (A[:,idx1],A[:,idx2])
    return A1,A2,idx1,idx2
end

"""
x_star is the n-vector by which the problem must
be translated in the first step of the temporal
instanton QCQP solution.

x_star is chosen to be the point in the set Ax=b
nearest to the origin such that x_star[idx3] = 0.
This condition ensures no linear term is introduced
into the quadratic constraint.
"""
function translation_point(
    A1::SparseMatrixCSC{Float64,Int64},
    A2::SparseMatrixCSC{Float64,Int64},
    idx1::Vector{Int64},
    idx2::Vector{Int64},
    n::Int64,
    b::Vector{Float64}
    )
    x_star = zeros(n)
    x_star[[idx1;idx2]] = [A1 A2]\b
    return x_star
end

"""
This function performs the change of variables from x to z,
where z = x - x_star. (For translating a quadratic problem.)
Returns triple `H_of_x` consisting of matrix H, vector h, constant kh.

Arguments:

* `G_of_x` consists of matrix G, vector g, constant kg.
* `x_star` is translation.

Used to perform second step of temporal instanton solution method,
assuming `x_star` is min-norm solution of Ax=b.

To save time, this method does not check for dimension mismatches.
"""
function translate_quadratic(
    G_of_x::Tuple{SparseMatrixCSC{Float64,Int64},Vector{Float64},Float64},
    x_star::Vector{Float64}
    )
    G,g,kg = G_of_x
    H = G
    h = g + 2*G*x_star
    kh = kg + x_star'*(G*x_star) + g'*x_star
    return (H,h,kh[1])
end

"""
Find an orthonormal basis for the nullspace of A.
This matrix may be used to rotate a temporal instanton
problem instance to eliminate all but nullity(A) elements.
"""
function kernel_basis(A::SparseMatrixCSC{Float64,Int64})
    m,n = size(A)
    # Assume A has full row rank
    dim_N = n - m
    F = qrfact(A')
    # B selects last dim_N cols of Q:
    B = [zeros(n-dim_N,dim_N); eye(dim_N)]
    N = sparse(SparseMatrix.SPQR.qmult(SparseMatrix.SPQR.QX, F, SparseMatrix.CHOLMOD.Dense(B)))
    return N
end

"""
Return a basis for the null space of a dense rectangular m-by-n
matrix with rank m. Based on "back substitution algorithm"
as described in Berry 1985 (**DOI: 10.1007/BF01389453**). Input
`q` represents a column order such that the first m columns of
`A[:,q]` are linearly independent.

*Note: if `q` does not cause `A[:,1:m]` to have full rank, U1 will be
singular. The method will fail with a `LAPACKException(k)`, where
`k` is the index of the first zero on the diagonal of U1.*
"""
function kernel_backsubs(
    A::SparseMatrixCSC{Float64,Int64},
    q=1:size(A,2)
    )
    m,n = size(A)

    F = lufact(A[:,q])
    qlu = F[:q]
    U = F[:U][:,invperm(qlu)]
    U1 = U[:,1:m]
    U2 = U[:,m+1:end]

    W = U1\full(U2)
    B = [-W;eye(n-m)][invperm(q),:] # equiv. to
    # pre-multiplying by pvec2mat(invperm(q))
end

"""
Return a column permutation q for a temporal instanton
A matrix such that the first m columns of A form a square,
nonsingular matrix. Useful for finding a null space basis
matrix via the back substitution method.
"""
function permutecols(nr,nb,nt)
    m = (nb+2)*nt
    n = (nr+nb+2)*nt
    q = Vector{Int64}()
    # take non-decision-var cols across all time steps
    for t in 1:nt
        append!(q,collect(nr+1:nr+nb+1)+(t-1)*(nr+nb+1))
    end
    # append aux. angle variable cols
    append!(q,collect(n-nt+1:n))
    q = [q;setdiff(1:n,q)]
    return q
end

"""
Rotate quadratic `G_of_x` by
rotation matrix `R`. Sparse or dense.
"""
function rotate_quadratic{T<:AbstractArray}(
    G_of_x::Tuple{T,Vector{Float64},Float64},
    R::Union{Array{Float64,2},SparseMatrixCSC{Float64,Int64}}
    )
    G,g,kg = G_of_x
    return (R*(G*R'),R*g,kg)
end

"""
Return `K`, the diagonal matrix whose elements are
square roots of elements in the given vector `d`.
"""
function return_K(d::Vector{Float64})
    K = ones(length(d))
    K[find(d)] = sqrt(d[find(d)])
    return spdiagm(K)
end

"""
Use diagonal elements of `Q_of_w[1]` to divide `G_of_w[1]`
into four blocks and `G_of_w[2]` into two blocks.
"""
function partition_B(G_of_w::Tuple,Q::Vector{Float64})
    B,b = G_of_w[1],G_of_w[2]
    i2 = find(Q)
    i1 = setdiff(1:length(Q),i2)
    # i1 = find(diag(Q))
    # i2 = setdiff(1:size(Q,1),i1)
    B11,B12,B21,B22 = B[i1,i1],B[i1,i2],B[i2,i1],B[i2,i2]
    b1 = b[i1]
    b2 = b[i2]
    return B11,B12,B22,b1,b2,i1,i2
end

"""
Obtain Bhat used to eliminate w1.
"""
function return_Bhat(
    B11::Array{Float64,2},
    B12::Array{Float64,2},
    B22::Array{Float64,2},
    b1::Vector{Float64},
    b2::Vector{Float64}
    )
    temp = B12'/B11
    Bhat = B22 - temp*B12
    bhat = b2 - temp*b1
    return Bhat,bhat
end

# """
# Eliminate w2
# """
# function return_Bhat(
#     B11::Array{Float64,2},
#     B12::Array{Float64,2},
#     B22::Array{Float64,2},
#     b1::Vector{Float64},
#     b2::Vector{Float64}
#     )
#     temp = B12/B22
#     Bhat = B11 - temp*(B12')
#     bhat = -(b1 - temp*b2)/2
#     return Bhat,bhat
# end

# """
# Reverse rotations and translations to map
# secular equation solution back to original problem
# space.
# """
# function return_xopt(
#     w2opt::Vector{Float64},
#     B11::Array{Float64,2},
#     B12::Array{Float64,2},
#     b1::Vector{Float64},
#     N::SparseMatrixCSC{Float64,Int64},#Array{Float64,2},
#     U::Array{Float64,2},
#     K::SparseMatrixCSC{Float64,Int64},#Array{Float64,2},
#     x_star::Vector{Float64}
#     )
#     w1opt = -B11\(B12*w2opt + b1/2)
#     wopt = [w1opt;w2opt]
#     xopt = N*sparse(U)*spdiagm(1./diag(K))*wopt + x_star
#     return xopt
# end

"""
    v,w = solvesecular(d,D,s2)
Solve the explicit secular equation
by returning the smallest `v` satisfying

f(v) = w'\*w - s2 = 0,

where `w = d./(D - v)`.

Inputs:

* `d` is a vector with numerator values
* `D` is a vector of poles (no need to sort them)
* `s2` is the target constant; we return solutions s(v) = s2.

Output:

* `v` is the optimal `v` satisfying s(v) = c
* `w` is the vector corresponding to `v`
"""
function solvesecular(
    d::Vector{Float64},
    D::Vector{Float64},
    s2::Float64
    )
    # sort eigenvalues by magnitude:
    p = sortperm(D)
    D = D[p]
    d = d[p]
    # explicit secular eq
    f(v::Float64) = sumabs2(d./(D-v)) - s2
    fprime(v::Float64) = sum((2*d.^2)./((D-v).^3))
    # vector in terms of v
    w(v::Float64) = float(d./(v - D))[invperm(p)]

    if isempty(find(d))
        println("d = $d, no secular equation solution!")
        return 0., w(0.)
    end
    # index of first nonzero of d
    k = findfirst(d)

    iterates = Vector{Float64}()
    v0 = D[k] - abs(d[k])/sqrt(s2)
    push!(iterates,v0)

    stop = false
    while !stop
        v = iterates[end]
        vnew = v - 2*((f(v) + s2)/fprime(v))*(sqrt((f(v)+s2)/s2) - 1)
        push!(iterates,vnew)
        stop = vnew >= v
    end
    v = iterates[end]
    return v, w(v)
end

"""
Return A[i2,i2] - A[i2,i1] x inv(A[i1,i1]) x A[i1,i2].

Uses block LU decomposition.
"""
function schurcomp(A,i1,i2)
    A11 = lufact(A[i1,i1])
    if issparse(A)
        L11,U11,p,q,Rs = A11[:(:)]
        P = pvec2mat(p)
        Q = pvec2mat(q)
        Rs = diagm(Rs)

        L21 = (A[i2,i1]*Q)/UpperTriangular(U11)
        U12 = LowerTriangular(L11)\(P*Rs*A[i1,i2])

        S = A[i2,i2] - L21*U12
        return S,(P,Q,Rs)
    else
        P = pvec2mat(A11[:p])
        println(A11[:p])
        L11 = A11[:L]
        U11 = A11[:U]

        L21 = A[i2,i1]/UpperTriangular(U11)
        U12 = LowerTriangular(L11)\(P*A[i1,i2])

        S = A[i2,i2] - L21*U12
        return S
    end
end

"""
Given a permutation vector, return the corresponding permutation matrix such that
`A[p,:] = P*A`.
See [Wikipedia](https://en.wikipedia.org/wiki/Permutation_matrix#Definition).
"""
function pvec2mat(p)
    n = length(p)
    P = zeros(n,n)
    for i in 1:n
        P[i,:] = ej(n,p[i])
    end
    return sparse(P)
end

"""
Return the j-th standard basis vector with length n.
"""
function ej(n,j)
    v = zeros(n)
    v[j] = 1.0
    return v
end
