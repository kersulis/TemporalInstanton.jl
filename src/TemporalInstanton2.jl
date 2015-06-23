module TemporalInstanton2

export
    partition_A,find_x_star,translate_quadratic,rotate_quadratic,
    kernel_rotation,return_K,partition_B,return_Bhat,find_w,
    solve_secular,return_xopt,solve_temporal_instanton,loop_through_lines,
    loop_through_lines_new

push!(LOAD_PATH, dirname(@__FILE__))
using TemporalInstanton
include("LineThermalModel.jl")

function partition_A(A,Qobj,T)
    """ Return A1, A2, A3 where
    A1 corresponds to wind
    A2 corresponds to angles + mismatch
    A3 corresponds to angle difference vars
    """
    m,n = size(A)
    idx1 = find(diag(Qobj))
    idx2 = setdiff(1:n-T,idx1)
    idx3 = n-T+1:n
    idx = [idx1,idx2,idx3]
    
    A1 = A[:,idx1]
    A2 = A[:,idx2]
    A3 = A[:,idx3]
    return A1,A2,A3,idx1,idx2,idx3
end

function find_x_star(Qobj,A,b,T)
    A1,A2,A3,idx1,idx2,idx3 = partition_A(A,Qobj,T)
    x_star = zeros(size(A,2))
    x_star[[idx1,idx2]] = [A1 A2]\b
    return x_star
end

function translate_quadratic(G_of_x,x_star)
    """ This function performs the change of variables from x to z,
    where z = x - x_star. (For translating a quadratic problem.)
    Returns triple H_of_x consisting of matrix H, vector h, constant kh.

    Arguments
    G_of_x consists of matrix G, vector g, constant kg.
    x_star is translation.

    Assume x_star is min-norm solution of Ax=b.
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

function rotate_quadratic(G_of_x,R)
    """ Rotate quadratic G_of_x by
    rotation matrix R.
    """
    G,g,kg = G_of_x
    return (R*G*R',R*g,kg)
end

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
    return R
end

function return_K(D)
    K = ones(length(D))
    K[find(D)] = sqrt(D[find(D)])
    K = diagm(K)
    return K
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

function solve_secular(D,d,c)
    """ Solve the secular equation via binary search.
    """
    eps = 1e-8
    solutions = Float64[]
    vectors = Array(Vector{Float64},0)
    poles = sort(unique(round(diag(D),10)))
    
    # Each diagonal element is a pole.
    for i in 1:length(poles)
        
        # Head left first:
        high = poles[i]
        if length(poles) == 1
            low = high - high
        elseif i == 1
            low = high - abs(poles[i] - poles[i+1])
        else
            low = high - abs(poles[i] - poles[i-1])/2
        end
        
        # Initialize v:
        v = (high + low)/2
        w = find_w(v,D,d)
        diff = (w'*w)[1] - c
        diff_old = 0
        stall = false
        while abs(diff) > eps
            if diff == diff_old
                stall = true
                break
            end
            if diff > 0
                high = v
            else
                low = v
            end
            v = (high + low)/2
            w = find_w(v,D,d)
            diff_old = diff
            diff = (w'*w)[1] - c
        end
        if !stall
            push!(solutions,v)
            push!(vectors,w)
        end
        
        # Now head right:
        high = poles[i]
        if length(poles) == 1
            low = high + high
        elseif i == length(poles)
            low = high + abs(poles[i] - poles[i-1])
        else
            low = high + abs(poles[i] - poles[i+1])/2
        end
        
        v = (high + low)/2
        w = find_w(v,D,d)
        diff = (w'*w)[1] - c
        diff_old = 0
        stall = false
        while abs(diff) > eps
            if diff == diff_old
                stall = true
                break
            end
            if diff > 0
                high = v
            else
                low = v
            end
            v = (high + low)/2
            w = find_w(v,D,d)
            diff_old = diff
            diff = (w'*w)[1] - c
        end
        if !stall
            push!(solutions,v)
            push!(vectors,w)
        end
    end
    return solutions,vectors
end

function return_xopt(w2opt,B11,B12,b1,N,U,K,x_star)
    w1opt = -B11\(B12*w2opt + b1/2)
    wopt = [w1opt,w2opt]
    xopt = (N*U/K)*wopt + x_star
    return xopt
end

function solve_temporal_instanton(G_of_x,Q_of_x,A,b,T)
    
    Qobj = G_of_x[1]
    c = - Q_of_x[3]
    
    opt = Array(Vector{Float64},0)
    
    # Partition A:
    A1,A2,A3,idx1,idx2,idx3 = partition_A(A,Qobj,T)

    # Find translation point:
    x_star = find_x_star(Qobj,A,b,T)

    # Translate quadratics:
    G_of_y = translate_quadratic(G_of_x,x_star)
    Q_of_y = translate_quadratic(Q_of_x,x_star)

    N = kernel_rotation(A)[:,1:size(A,2) - rank(A)] # take only first k cols

    N1,N2,N3 = N[idx1,:],N[idx2,:],N[idx3,:] # partition N

    G_of_z = rotate_quadratic(G_of_y,N')
    Q_of_z = rotate_quadratic(Q_of_y,N')

    D,U = eig(Q_of_z[1])
    D = round(D,10)

    K = return_K(D)

    G_of_w = rotate_quadratic(G_of_z,(U/K)')
    Q_of_w = rotate_quadratic(Q_of_z,(U/K)')

    B11,B12,B21,B22,b1,b2 = partition_B(G_of_w,Q_of_w)

    Bhat,bhat = return_Bhat(B11,B12,B22,b1,b2)

    eps = 1e-8
    w0 = find_w(0,Bhat,bhat/2)

    if abs((w0'*w0) - c)[1] < eps
        println("v=0 works!")
    end

    solutions, vectors = solve_secular(Bhat,bhat/2,-Q_of_w[3])
    if isempty(solutions)
        return [],Inf
    else
        sol = zeros(length(vectors))
        for i in 1:length(vectors)
            w2 = vectors[i]
            xvec = return_xopt(w2,B11,B12,b1,N,U,K,x_star)
            sol[i] = (xvec'*Qobj*xvec)[1]
            push!(opt,xvec)
        end
    end

    return opt[indmin(sol)],minimum(sol)
end

function loop_through_lines(
    G_of_x,
    Q_of_x,
    A1,
    b,
    n,
    T,
    tau,
    Ridx,
    ref,
    lines)

    nr = length(Ridx)
    
    numLines = length(lines)

    # Initialize vars used to store results:
    score = Float64[]
    α = Array(Vector{Float64},0)

    θ = Array(Array,0)
    x = Array(Array,0)
    diffs = Array(Array,0)
    xopt = Array(Array,0)

    # Loop through all lines:
    for line in lines
        # array of vec. with Float values:
        deviations = Array(Vector{Float64},0)
        angles = Array(Vector{Float64},0)
        alpha = Float64[]

        # Create A2 based on chosen line:
        A2 = tmp_inst_A_scale(n,Ridx,T,tau,ref,line)
        # Stack A1 and A2:
        A = [A1; A2]

        xvec,sol = solve_temporal_instanton(G_of_x,Q_of_x,A,b,T)
        push!(score,sol)
        if isinf(sol)
            push!(deviations,[])
            push!(angles,[])
            push!(alpha,NaN)
            push!(diffs,[])
        else

            # Variable breakdown:
            # (nr+n+1) per time step
            #   First nr are deviations
            #   Next n are angles
            #   Last is mismatch
            # T variables at the end: anglediffs
            for t = 1:T
                push!(deviations,xvec[(nr+n+1)*(t-1)+1:(nr+n+1)*(t-1)+nr])
                push!(angles,xvec[(nr+n+1)*(t-1)+nr+1:(nr+n+1)*(t-1)+nr+n])
                push!(alpha,xvec[(nr+n+1)*(t)])
            end
            push!(diffs,xvec[end-T+1:end])
        end

        push!(x,deviations)
        push!(θ,angles)
        push!(α,alpha)
        push!(xopt,xvec)
    end
    return score,x,θ,α,diffs,xopt
end

function loop_through_lines_new(
    G_of_x,
    Qtheta,
    A1,
    b,
    n,
    T,
    Ridx,
    Sb,
    ref,
    lines,
    res,
    reac,
    line_lengths,
    Tamb,
    T0,
    int_length)

    nr = length(Ridx)
    
    numLines = length(lines)

    # Initialize vars used to store results:
    score = Float64[]
    α = Array(Vector{Float64},0)

    θ = Array(Array,0)
    x = Array(Array,0)
    diffs = Array(Array,0)
    xopt = Array(Array,0)

    # Loop through all lines:
    for idx = 1:numLines
        # thermal model cannot handle zero-length lines:
        if line_lengths[idx] == 0
            continue
        end
        line = lines[idx]
        line_model = LineModel(line[1],
                    line[2],
                    res[idx],
                    reac[idx],
                    line_lengths[idx],
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN,
                    NaN)

        add_thermal_parameters(line_model, "waxwing")

        therm_a = compute_a(line_model.mCp,
                            line_model.ηc,
                            line_model.ηr,
                            Tamb,
                            line_model.Tlim)
        therm_c = compute_c(line_model.mCp,
                            line_model.rij,
                            line_model.xij,
                            Sb,
                            line_model.length)
        therm_d = compute_d(line_model.mCp,
                            line_model.ηc,
                            line_model.ηr,
                            Tamb,
                            line_model.Tlim,
                            line_model.qs)
        therm_f = compute_f(int_length,
                            therm_a,
                            therm_d,
                            T,
                            T0)
        
        # thermal constraint, Q(z) = 0:
        kQtheta = (therm_a/therm_c)*(line_model.Tlim - therm_f)
        Q_of_x = (Qtheta,0,kQtheta)

        # array of vec. with Float values:
        deviations = Array(Vector{Float64},0)
        angles = Array(Vector{Float64},0)
        alpha = Float64[]

        # Create A2 based on chosen line:
        A2 = tmp_inst_A_scale_new(n,Ridx,T,line,therm_a,int_length)
        # Stack A1 and A2:
        A = [A1; A2]

        xvec,sol = solve_temporal_instanton(G_of_x,Q_of_x,A,b,T)
        push!(score,sol)
        if isinf(sol)
            push!(deviations,[])
            push!(angles,[])
            push!(alpha,NaN)
            push!(diffs,[])
        else

            # Variable breakdown:
            # (nr+n+1) per time step
            #   First nr are deviations
            #   Next n are angles
            #   Last is mismatch
            # T variables at the end: anglediffs
            for t = 1:T
                push!(deviations,xvec[(nr+n+1)*(t-1)+1:(nr+n+1)*(t-1)+nr])
                push!(angles,xvec[(nr+n+1)*(t-1)+nr+1:(nr+n+1)*(t-1)+nr+n])
                push!(alpha,xvec[(nr+n+1)*(t)])
            end
            push!(diffs,xvec[end-T+1:end])
        end

        push!(x,deviations)
        push!(θ,angles)
        push!(α,alpha)
        push!(xopt,xvec)
    end
    return score,x,θ,α,diffs,xopt
end

end