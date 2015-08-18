""" Return A1, A2, A3 where:
* A1 corresponds to wind
* A2 corresponds to angles + mismatch
* A3 corresponds to angle difference vars

Used to find x_star, the min-norm solution to
Ax=b such that x_star[idx3] = 0.
"""
function partition_A(A,Qobj,T)
    n = size(A,2)
    idx1 = find(diag(Qobj))
    idx2 = setdiff(1:n-T,idx1)
    idx3 = n-T+1:n

    (A1,A2) = (A[:,idx1],A[:,idx2])
    return A1,A2,idx1,idx2,idx3
end

""" x_star is the n-vector by which the problem must
be translated in the first step of the temporal
instanton QCQP solution.

x_star is chosen to be the point in the set Ax=b
nearest to the origin such that x_star[idx3] = 0.
This condition ensures no linear term is introduced
into the quadratic constraint.
"""
function find_x_star(A1,A2,idx1,idx2,n,b)
    x_star = zeros(n)

    try
        x_star[[idx1;idx2]] = sparse([A1 A2])\b
    catch
        # if [A1 A2] has LD rows, use A2 only
        A2 = sparse(A2)
        x_star[idx2] = (b'*A2/(A2'*A2))'
    end
    # x_star[idx2] = sparse(A2[1:end-T,:])\(b[1:end-T])

    return x_star
end

""" This function performs the change of variables from x to z,
where z = x - x_star. (For translating a quadratic problem.)
Returns triple H_of_x consisting of matrix H, vector h, constant kh.

Arguments
G_of_x consists of matrix G, vector g, constant kg.
x_star is translation.

Used to perform second step of temporal instanton solution method,
assuming x_star is min-norm solution of Ax=b.
"""
function translate_quadratic(G_of_x,x_star)
    G,g,kg = G_of_x
    if g == 0
        g = zeros(size(G,1),1)
    end
    H = G
    h = g + 2*G*x_star
    kh = kg + x_star'*G*x_star + g'*x_star
    return (H,h,kh[1])
end

""" Find an orthonormal basis for the nullspace of A.
This matrix may be used to rotate a temporal instanton
problem instance to eliminate all but nullity(A) elements.
"""
function kernel_rotation(A; spqr=true)
    m,n = size(A)

    # Assume A always has full row rank of m.
    # It may be possible for this assumption to fail
    # due to numerics, but a rank() check is expensive.
    dim_N = n - m # dimension of nullspace of A

    if spqr
        F = qrfact(sparse(A'))
        # B selects last dim_N cols of Q:
        B = [zeros(size(A,2)-dim_N,dim_N); eye(dim_N)]
        N = sparse(SparseMatrix.SPQR.qmult(SparseMatrix.SPQR.QX, F, SparseMatrix.CHOLMOD.Dense(B)))
        return N
    else
        q = qr(A'; thin=false)[1]
        return q[:,end-dim_N+1:end]
    end
end

""" Rotate quadratic G_of_x by
rotation matrix R.
"""
function rotate_quadratic(G_of_x,R)
    G,g,kg = G_of_x
    return (R*G*R',R*g,kg)
end

function rotate_matrix(G,R)
    R*G*R'
end

function rotate_vector(G,R)
    R*g
end

""" Return K, the diagonal matrix whose elements are
square roots of eigenvalues of the given matrix D.
"""
function return_K(D)
    K = ones(length(D))
    K[find(D)] = sqrt(D[find(D)])
    return spdiagm(K)
end

function partition_B(G_of_w,Q_of_w)
    B,b = G_of_w[1],G_of_w[2]
    Q = round(Q_of_w[1])
    i2 = find(diag(Q))
    i1 = setdiff(1:size(Q,1),i2)
    B11,B12,B21,B22 = B[i1,i1],B[i1,i2],B[i2,i1],B[i2,i2]
    b1 = b[i1]
    b2 = b[i2]
    return B11,B12,B21,B22,b1,b2
end

function return_Bhat(B11,B12,B22,b1,b2)
    Bhat = B22 - (B12'/B11)*B12
    bhat = b2 - (B12'/B11)*b1
    return round(Bhat,10),bhat
end

function find_w(v,D,d)
    if v == 0
        w = float([-d[i]/(D[i,i]) for i in 1:length(d)])
    else
        w = float([d[i]/(v - D[i,i]) for i in 1:length(d)])
    end
    return w
end

""" Reverse rotations and translations to map
secular equation solution back to original problem
space.
"""
function return_xopt(w2opt,B11,B12,b1,N,U,K,x_star)
    w1opt = -B11\(B12*w2opt + b1/2)
    wopt = [w1opt;w2opt]
    xopt = N*U*diagm(1./diag(K))*wopt + x_star
    return xopt
end
