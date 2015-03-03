module TrustRegionSubproblem

export min_norm, translate_quadratic, kernel_rotation

function min_norm(A::Array{Float64,2},b::Array{Float64,1})
    """ Find the point x* in the set Ax=b which is
    closest to the origin in the 2-norm sense.
    """
    if size(A,1) < size(A,2)
        x_star = \(A,b)
    else
        println("need m < n")
    end
end

function translate_quadratic(
    G::Array{Float64,2},
    g::Array{Float64,1},
    kg::Float64,
    x_star::Array{Float64,1},
    )
    """ This function performs the change of variables from x to z,
    where z = x - x_star. (For translating a quadratic problem.)
    """
    H = G
    h = g - 2*G*x_star
    kh = kg - x_star'*G*x_star - g'*x_star
    return (H,h,kh[1])
end

function kernel_rotation(A::Array{Float64})
    """ Find an orthonormal basis for the 
    nullspace of A. This matrix may be used
    to rotate a temporal instanton problem
    instance.
    """
    m,n = size(A)
    dim_N = n - rank(A)
    # if A has full row rank:
    # dim_N = n - m
    q = qr(A'; thin=false)[1]
    R = circshift(q,(0,dim_N))
end