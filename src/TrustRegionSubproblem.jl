module TrustRegionSubproblem

export min_norm, translate_quadratic, kernel_rotation, tr_translate, tr_kernel_rotate, tr_diag_rotate

# function min_norm(A::Array{Float64,2},b::Array{Float64,1})
function min_norm(A,b)
    """ Find the point x* in the set Ax=b which is
    closest to the origin in the 2-norm sense.
    """
    m,n = size(A)
    if m < n
        x_star = \(A,b)
    else
        println("need m < n")
    end
end

# function translate_quadratic(
#     G::Array{Float64,2},
#     g::Array{Float64,1},
#     kg::Float64,
#     x_star::Array{Float64,1},
#     )
function translate_quadratic(
    G_of_x,x_star)
    """ This function performs the change of variables from x to z,
    where z = x - x_star. (For translating a quadratic problem.)
    Returns triple H_of_x consisting of matrix H, vector h, constant kh.

    Arguments
    G_of_x consists of matrix G, vector g, constant kg.
    x_star is translation.
    """
    G,g,kg = G_of_x
    if g == 0
        g = zeros(size(G,1),1)
    end
    H = G
    h = g - 2*G*x_star
    kh = kg - x_star'*G*x_star - g'*x_star
    return (H,h,kh[1])
end

# function kernel_rotation(A::Array{Float64})
function kernel_rotation(A)
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

# Functions added 2015-03-07
function tr_translate(G_of_x,Q_of_x,A,b)
    """ Starting with min G(x) s.t. ||x||=c, Ax=b,
    translate problem by x_star, the point in 
    {x:Ax=b} closest to the origin. The problem
    becomes min H(x) s.t. ||x||=c', Ax=0.
    (Note: Q(x) == ||x||.)
    
    Returns H_of_x (translated G) and R_of_x
    (translated Q).
    
    Arguments:
    G_of_x has form (G,g,kg), with G matrix, g vector,
    kg constant.
    Q_of_x is entered similarly.
    A and b define the Ax=b relationship.
    """
    x_star = min_norm(A,b)
    
    H_of_x = translate_quadratic(G_of_x,x_star)
    R_of_x = translate_quadratic(Q_of_x,x_star)
    
    return H_of_x,R_of_x
end

function tr_kernel_rotate(G_of_x,A)
    """ Rotate quadratic G_of_x to eliminate all
    but first k variables. Returns quadratic H_of_z
    """
    G,g,kg = G_of_x
    rotation = kernel_rotation(A)
    return (rotation*G,rotation*g,kg)
end

function tr_diag_rotate(G_of_x)
    """ Change variables to obtain quadratic polynomial
    whose quadratic part is diagonal.
    """
    G,g,kg = G_of_x
    D,V = eig(G)
    return (diagm(D),V'*g,kg)
end

end