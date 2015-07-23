# The following five functions are used to run power flow
# with fixed decision variable values.

""" Expand renewable generation vector with zeros.
"""
function expand_renewable_vector(x,Ridx,N,T)
    idx = Array(Integer,0)
    for i = 0:T-1
        append!(idx,Ridx + i*N)
    end
    P = zeros(N*T)
    P[idx] = x
    return P
end

""" Generate the power balance constraint A matrix
from problem dimensions, admittance matrix,
and generator participation factors.
Assumes the admittance matrix is n-by-n.

Returns A, which is (n+1)*T-by-(nr+n+1)*T

* nr is the number of wind farms in the network
* n is the number of nodes in the network
* Ridx is a vector indicating wind farm locations
* T is the number of time steps
* Y is the admittance matrix (n-by-n)
* ref is the index of the angle reference bus
* k is the vector of generator participation factors
"""
function fixed_wind_A(T,Y,ref,k)
    function ei(n,i)
        e = zeros(n)
        e[i] = 1.
        return e
    end

    n = size(Y,1)

    # A has a block diagonal pattern where each
    # block is Atemp:
    Atemp = sparse([[  Y       -k];
                    ei(n,ref)'   0])

    # Now we can tile the Atemp matrix to generate A:
    A = Atemp
    for t = 2:T
        A = blkdiag(A, Atemp)
    end

    return full(A)
end

""" Generate the vector b of power balance constraints.
Assumes G0 and D are nT-by-1 vectors.
"""
function fixed_wind_b(n,T,G0,Pnet,D)
    b = FloatingPoint[]
    netGen = G0 + Pnet - D

    for t = 1:T
        start = (t-1)*n + 1
        stop = start + n - 1
        append!(b,netGen[start:stop])
        push!(b,0.)
    end
    return b
end

function return_angles(fixed_x,N)
    angles = Array(Vector,0)
    alpha = FloatingPoint[]
    for i = 1:T
        push!(angles,fixed_x[(N+1)*(i-1)+1:(N+1)*(i-1)+N])
        push!(alpha,fixed_x[(N+1)*(i-1)+N+1])
    end
    return angles,alpha
end

function return_angle_diffs(angles,line)
    angle_diffs = FloatingPoint[]
    f = line[1]
    t = line[2]
    for v in angles
        push!(angle_diffs,v[f] - v[t])
    end
    return angle_diffs
end
