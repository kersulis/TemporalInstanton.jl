module TrustRegionSubproblem

export
    min_norm, translate_quadratic, kernel_rotation, 
    tr_translate, tr_kernel_rotate, tr_diag_rotate, 
    tr_map_back, tr_trans_rotate, tr_check_diag,
    tr_solve_secular, find_w

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
    h = g + 2*G*x_star
    kh = kg + x_star'   *G*x_star + g'*x_star
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
    
    return H_of_x,R_of_x,x_star
end

function tr_kernel_rotate(G_of_x,A)
    """ Rotate quadratic G_of_x to eliminate all
    but first k variables. Returns quadratic H_of_z.

    Note: equivalent to change of variables from x 
    to z=R*x
    """
    G,g,kg = G_of_x
    R = kernel_rotation(A)
    return (R*G*R',R*g,kg),R
end

function tr_diag_rotate(G_of_x)
    """ Change variables to obtain quadratic polynomial
    whose quadratic part is diagonal.
    """
    G,g,kg = G_of_x
    D,U = eig(G)
    return (diagm(D),U'*g,kg),U'
end

function tr_map_back(w,Rkernel,Reigvec)
    return Rkernel'*Reigvec*w
end

function tr_trans_rotate(G_of_x,Q_of_x,A,b)
    H_of_x,R_of_x,x_star = tr_translate(G_of_x,Q_of_x,A,b)
    # R_of_x is rotation invariant.
    J_of_z,Rkernel = tr_kernel_rotate(H_of_x,A)
    K_of_w,Reigvec = tr_diag_rotate(J_of_z)
    
    return K_of_w,R_of_x,x_star,Rkernel,Reigvec
end

function find_w(mu,D,d,Qtheta)
    # Evaluate secular equation
    w = float([d[i]/(mu*Qtheta[i,i] - D[i,i]) for i in 1:length(d)])
end

function tr_check_diag(D,d,Qtheta,c)
    eps = 1e-8
    check = Bool[]
    w_vals = Array(Vector{Float64},0)
    mu_vals = FloatingPoint[]

    # Check zero:
    mu = 0
    w = find_w(mu,D,d,Qtheta)
    append!(check, [abs((w'*w)[1] - c^2) < eps])
    push!(w_vals,w)
    append!(mu_vals,[mu])

    for i = 1:size(D,1)
        mu = D[i,i]
        w = find_w(mu,D,d,Qtheta)
        append!(check, [abs((w'*w)[1] - c^2) < eps])
        push!(w_vals,w)
        append!(mu_vals,[mu])
    end
    return mu_vals,w_vals,check
end

function tr_solve_secular(D,d,Qtheta,c)
    """ Solve the secular equation via binary search.
    """
    eps = 1e-8
    solutions = Float64[]
    poles = sort(diag(D*Qtheta))
    # Each diagonal element is a pole.
    for i in 1:length(poles)
        
        # Head left first:
        high = poles[i]
        if i == 1
            low = high - abs(poles[i] - poles[i+1])
        else
            low = high - abs(poles[i] - poles[i-1])/2
        end
        
        # Initialize mu:
        mu = (high + low)/2
        w = find_w(mu,D,d,Qtheta)
        diff = (w'*w)[1] - c^2
        diff_old = 0
        stall = false
        while abs(diff) > eps
            if diff == diff_old
                stall = true
                break
            end
            if diff > 0
                high = mu
            else
                low = mu
            end
            mu = (high + low)/2
            w = find_w(mu,D,d,Qtheta)
            diff_old = diff
            diff = (w'*w)[1] - c^2
        end
        if !stall
            push!(solutions,mu)
        end
        
        # Now head right:
        high = poles[i]
        if i == length(poles)
            low = high + abs(poles[i] - poles[i-1])
        else
            low = high + abs(poles[i] - poles[i+1])/2
        end
        
        mu = (high + low)/2
        w = find_w(mu,D,d,Qtheta)
        diff = (w'*w)[1] - c^2
        diff_old = 0
        stall = false
        while abs((w'*w)[1] - c^2) > eps
            if diff == diff_old
                stall = true
                break
            end
            if diff > 0
                high = mu
            else
                low = mu
            end
            mu = (high + low)/2
            w = find_w(mu,D,d,Qtheta)
            diff_old = diff
            diff = (w'*w)[1] - c^2
        end
        if !stall
            push!(solutions,mu)
        end
    end
    return solutions
end

end