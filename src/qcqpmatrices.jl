using SparseArrays: sparse, blockdiag, spzeros, spdiagm, sparsevec, dropzeros

# """
# Generate the objective function matrix given problem dimensions.
# Output `Qobj` is square with dimension (nr + n + 2)T (equal to number of
# variables in QCQP).
# """
# function tmp_inst_Qobj(
#     n::Int64,
#     nr::Int64,
#     T::Int64,
#     corr::Array{Float64,2}=Matrix{Float64}(undef, 0, 0);
#     pad::Bool = true
#     )
#     if isempty(corr)
#         Qobj = spdiagm(0 => repeat([ones(nr); zeros(n + 1)]; outer=T))
#     else
#         prec = sparse(inv(corr))
#         Qobj = blockdiag(prec,spzeros(n + 1, n + 1))
#         for t in 2:T
#             Qobj = blockdiag(Qobj, blockdiag(prec, spzeros(n + 1, n + 1)))
#         end
#         # ensure the objective matrix has norm 1:
#         scale!(Qobj, 1 / norm(full(Qobj)))
#     end
#
#     if !pad
#         return Qobj
#     else
#         # Add T rows and columns of zeros to Q:
#         r, c = size(Qobj)
#         return [[Qobj spzeros(r, T)]; spzeros(T, c + T)]
#     end
# end

"""
Generate the objective function matrix given problem dimensions.
Output `Qobj` is square with dimension (nr + n + 2)T (equal to number of
variables in QCQP).
"""
function tmp_inst_Qobj(
    n::Int64,
    nr::Int64,
    T::Int64,
    auto_prec::Float64=0.0
    )

    Qobj = spdiagm(
        0 => repeat([ones(nr); zeros(n + 1)]; outer=T),
        (n + nr + 1) => repeat([auto_prec * ones(nr); zeros(n + 1)]; outer=(T - 1)),
        -(n + nr + 1) => repeat([auto_prec * ones(nr); zeros(n + 1)]; outer=(T - 1))
    )

    # Add T rows and columns of zeros to Q:
    r, c = size(Qobj)
    return dropzeros([[Qobj spzeros(r, T)]; spzeros(T, c + T)])
end

"""
Generate the power balance constraint A matrix
from problem dimensions, admittance matrix,
and generator participation factors.
Assumes the admittance matrix is n-by-n.

Returns A, which is (n+1)T -by- (nr+n+2)T
(or (n+1)T -by- (nr+n+1)T if pad=false)

* nr is the number of wind farms in the network
* n is the number of nodes in the network
* Ridx is a vector indicating wind farm locations
* T is the number of time steps
* Y is the admittance matrix (n-by-n)
* ref is the index of the angle reference bus
* k is the vector of generator participation factors
"""
function tmp_inst_A1(
    Ridx::Vector{Int64},
    T::Int64,
    Y::SparseMatrixCSC{Float64,Int64},
    ref::Int64,
    k::Vector{Float64};
    pad::Bool = true
    )
    function ei(n::Int64, i::Int64)
        v = spzeros(n, 1)
        v[i] = 1.0
        return v
    end

    n = size(Y, 1)

    # A has a block diagonal pattern where each
    # block is Atemp:
    Atemp = [[-spdiagm(0 => ones(n))    sparse(Y)       -sparsevec(k)];
            spzeros(1, n)   ei(n + 1, ref)'                 ]

    # Remove columns corresponding to non-wind nodes:
    Atemp = Atemp[:, [Ridx; (n + 1):(2 * n + 1)]]

    # Now we can tile the Atemp matrix to generate A:
    A = Atemp
    for t in 2:T
        A = blockdiag(A, Atemp)
    end

    if !pad
        return A
    else
        # pad A with T columns of zeros
        # (rows having coefficients for these new variables
        # depend on the line and are added during the big for loop)
        return [A spzeros((n + 1) * T, T)]
    end
end

"""
Generate the vector b of power balance constraints.
Assumes G0 and D are nT-by-1 vectors.
"""
function tmp_inst_b(
    n::Int64,
    T::Int64,
    G0::Vector{Float64},
    R0::Vector{Float64},
    D0::Vector{Float64};
    pad::Bool = true
    )
    b = Vector{Float64}()
    netGen = G0 + R0 - D0

    for t in 1:T
        startidx = (t - 1) * n + 1
        stopidx = startidx + n - 1
        append!(b, netGen[startidx:stopidx])
        push!(b, 0.0)
    end
    if pad
        # Not sure when I would ever not want to do this...
        append!(b, zeros(T))
    end
    return b
end

"""
Generate `Qconstr` in the temperature constraint of a
temporal instanton problem instance.
"""
function tmp_inst_Qconstr(
    n::Int64,
    nr::Int64,
    T::Int64
    )

    blockdiag(spzeros((nr + n + 1) * T, (nr + n + 1) * T), spdiagm(0 => ones(T)))
end

"""
Augment `A` with `T` additional rows relating
angle difference variables to angle variables.

Returns a T -by- (n+nr+2)T matrix that may be
concatenated with the output of `temp_inst_A`

Arguments:
* n is the number of nodes in the network
* Ridx is a vector indicating wind nodes
* T is the number of time steps
* tau is the thermal coefficient from IEEE 738
* slack is the index of the slack bus
* line is the pair (i, k) indicating the chosen
line
* therm_a is a constant defined by heating parameters
* int_length is interval length in seconds (e.g. 600
    for 10 minutes)
"""
function tmp_inst_A2(
    n::Int64,
    Ridx::Array{Int64,1},
    T::Int64,
    line::Tuple{Int64,Int64},
    therm_a::Float64,
    int_length
    )

    (i, k) = line
    nr = length(Ridx)

    pos = [(nr + n + 1) * (t - 1) + nr for t in 1:T]
    one_pos = collect((n + nr + 1) * T + 1:(n + nr + 2) * T)
    coefs = [sqrt(-exp(therm_a * int_length)^(T - t + 1) +
        exp(therm_a * int_length)^(T - t)) for t in 1:T]

    sparse(
    repeat(1:T; outer=3),
    [pos .+ i; pos .+ k; one_pos],
    [-coefs; coefs; ones(T)]
    )
end
