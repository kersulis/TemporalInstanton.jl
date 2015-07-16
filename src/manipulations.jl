@iprofile begin
function partition_A(A,Qobj,T)
    """ Return A1, A2, A3 where:
    * A1 corresponds to wind
    * A2 corresponds to angles + mismatch
    * A3 corresponds to angle difference vars

    Used to find x_star, the min-norm solution to
    Ax=b such that x_star[idx3] = 0.
    """
    m,n = size(A)
    idx1 = find(diag(Qobj))
    idx2 = setdiff(1:n-T,idx1)
    idx3 = n-T+1:n

    (A1,A2) = (A[:,idx1],A[:,idx2])
    return A1,A2,idx1,idx2,idx3
end

function find_x_star(A1,A2,idx1,idx2,n,b)
    """ x_star is the n-vector by which the problem must
    be translated in the first step of the temporal
    instanton QCQP solution.

    x_star is chosen to be the point in the set Ax=b
    nearest to the origin such that x_star[idx3] = 0.
    This condition ensures no linear term is introduced
    into the quadratic constraint.
    """
    x_star = zeros(n)
    Z = sparse([A1 A2]')
    x_star[[idx1;idx2]] = (Z/(Z'*Z))*b
    return x_star
end

function translate_quadratic(G_of_x,x_star)
    """ This function performs the change of variables from x to z,
    where z = x - x_star. (For translating a quadratic problem.)
    Returns triple H_of_x consisting of matrix H, vector h, constant kh.

    Arguments
    G_of_x consists of matrix G, vector g, constant kg.
    x_star is translation.

    Used to perform second step of temporal instanton solution method,
    assuming x_star is min-norm solution of Ax=b.
    """
    G,g,kg = G_of_x
    if g == 0
        g = zeros(size(G,1),1)
    end
    H = G
    h = g + 2*G*x_star
    kh = kg + x_star'*G*x_star + g'*x_star
    return (H,h,kh[1])
end

function kernel_rotation(A)
    """ Find an orthonormal basis for the nullspace of A.
    This matrix may be used to rotate a temporal instanton
    problem instance to eliminate all but nullity(A) elements.
    """
    m,n = size(A)

    # Assume A always has full row rank.
    #if isposdef(A*A')
    dim_N = n - m
    # else
    #     dim_N = n - rank(A)
    #     warn("A does not have full row rank.")
    # end
    q = qr(A'; thin=false)[1]
    R = circshift(q,(0,dim_N))
    return R
end

function rotate_quadratic(G_of_x,R)
    """ Rotate quadratic G_of_x by
    rotation matrix R.
    """
    G,g,kg = G_of_x
    return (R*G*R',R*g,kg)
end

function return_K(D)
    """ Return K, the diagonal matrix whose elements are
    square roots of eigenvalues of the given matrix D.
    """
    K = ones(length(D))
    K[find(D)] = sqrt(D[find(D)])
    return diagm(K)
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

function return_xopt(w2opt,B11,B12,b1,N,U,K,x_star)
    """ Reverse rotations and translations to map
    secular equation solution back to original problem
    space.
    """
    w1opt = -B11\(B12*w2opt + b1/2)
    wopt = [w1opt;w2opt]
    #xopt = (N*U/K)*wopt + x_star
    xopt = N*U*diagm(1./diag(K))*wopt + x_star
    return xopt
end

end # @iprofile begin
