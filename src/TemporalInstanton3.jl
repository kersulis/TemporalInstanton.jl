module TemporalInstanton3

using HDF5, JLD, ProgressMeter

export
    solve_instanton_qcqp, solve_temporal_instanton, LineModel


# THERMAL MODELING
type LineModel
    from
    to
    rij
    xij
    length
    D0
    mCp
    Ilim
    r
    Tlim
    ηc
    ηr
    qs
end

function add_thermal_parameters(line_model,conductor_name)
    """ Assign values to fields of LineModel type instance.
    Uses data from Mads's MPC paper.
    """
    if conductor_name == "waxwing"
        line_model.D0   = 15.5e-3
        line_model.mCp  = 383.
        line_model.Ilim = 439.
        line_model.r    = 110e-6
        line_model.Tlim = 65.
        line_model.ηc   = 0.955
        line_model.ηr   = 2.207e-9
        line_model.qs   = 14.4
    elseif conductor_name == "dove"
        line_model.D0   = 23.5e-3
        line_model.mCp  = 916.
        line_model.Ilim = 753.
        line_model.r    = 60e-6
        line_model.Tlim = 69.
        line_model.ηc   = 1.179
        line_model.ηr   = 3.346e-9
        line_model.qs   = 21.9
    end
    return line_model
end

function compute_a(mCp,ηc,ηr,Tamb,Tlim)
    """ Returns constant a [1/s]
    mCp [J/m-C] is line heat capacity
    ηc [W/m-C] is conductive heat loss rate coefficient
    ηr [W/m-C^4] is radiative heat loss rate coefficient
    Tamb [C] is ambient temperature (of air)
    Tlim [C] is highest allowable line temperature
    """
    Tmid = (Tamb + Tlim)/2
    return mCp\(-ηc - 4*ηr*(Tmid + 273)^3)
end

function compute_c(mCp,r,x,Sb,L)
    """ Return constant c [W/m]
    mCp [J/m-C] is line heat capacity
    r [pu] is line resistance
    x [pu] is line reactance
    Sb [W] is system base MVA
    L [m] is line length
    """
    return r*Sb/(3*mCp*L*(x^2))
end

function compute_d(mCp,ηc,ηr,Tamb,Tlim,q_solar)
    """ Returns constant d [W/m]
    mCp [J/m-C] is line heat capacity
    ηc [W/m-C] is conductive heat loss rate coefficient
    ηr [W/m-C^4] is radiative heat loss rate coefficient
    Tamb [C] is ambient temperature (of air)
    Tlim [C] is highest allowable line temperature
    q_solar [W/m] is the solar heat gain rate
    """
    Tmid = (Tamb + Tlim)/2
    return mCp\(ηc*Tamb - ηr*((Tmid + 273)^4 - (Tamb + 273)^4) + 4*ηr*Tmid*(Tmid+273)^3 + q_solar)
end

function compute_f(int_length,a,d,n,T0)
    """ Returns constant f
    int_length [s] is length of each interval
    a [1/s] is a constant
    d [W/m] is a constant
    n [-] is the number of time intervals
    T0 [C] is the initial steady-state line temp
    """
    sum_coeff = sum([(e^(int_length*a))^i - (e^(int_length*a))^(i-1) for i in 1:n])
    return (e^(int_length*a))^n*T0 + (d/a)*sum_coeff
end

# TEMPORAL INSTANTON
function tmp_inst_Qobj(n,nr,T)
    """ Generate the objective function matrix
    Qobj from the problem dimensions.
    Assume no correlation between wind sites.
    """
    Qobj = sparse(diagm(repeat([ones(nr),zeros(n+1)],outer=[T])))
    #Qobj = tmp_inst_pad_Q(full(Qobj),T)
    return Qobj
end

function tmp_inst_pad_Q(Q,T)
    """ Add T rows and T columns of zeros
    to matrix Q
    """
    m,n = size(Q)
    return [[Q zeros(m,T)]; zeros(T,n+T)]
end

function tmp_inst_A(Ridx,T,Y,ref,k)#,tau,line)
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
    
    function ei(n,i)
        e = zeros(n)
        e[i] = 1.
        return e
    end
    
    n = size(Y,1)
    
    # A has a block diagonal pattern where each
    # block is Atemp:
    Atemp = [[-eye(n)    Y       -k];
            zeros(1,n) ei(n,ref)' 0]
    
    # Remove columns corresponding to non-wind nodes:
    Atemp = sparse(Atemp[:,[Ridx,n+1:2*n+1]])
    
    # Now we can tile the Atemp matrix to generate A:
    A = Atemp
    for t = 2:T
        A = blkdiag(A, Atemp)
    end

    # These rows relate auxiliary angle variables to
    # original angle variables:
    #     A2 = tmp_inst_A_scale(n,Ridx,T,tau,ref,line)

    #     # Augment A with new rows:
    #     A = [[full(A) zeros((n+1)*T,T)]; A2]
    
    return A
end

function tmp_inst_b(n,T,G0,P0,D)
    """ Generate the vector b of power balance constraints.
    Assumes G0 and D are nT-by-1 vectors.
    """
    b = FloatingPoint[]
    netGen = G0 + P0 - D

    for t = 1:T
        start = (t-1)*n + 1
        stop = start + n - 1
        append!(b,netGen[start:stop])
        push!(b,0.)
    end

    # Extend b with T additional zeroes:
    # tmp_inst_pad_b(b,T)
    return b
end

function tmp_inst_pad_b(b,T)
    """ Append T zeros to vector b
    """
    append!(b,zeros(T))
end

function tmp_inst_Qtheta(n,nr,T)#,tau)
    """ Generate Q_theta in the temperature constraint
    of a temporal instanton problem instance.
    "line" has the form (i,k), where i and k refer to
    the endpoints of the chosen line.
    """
    Qtheta = zeros((nr+n+2)*T,(nr+n+2)*T)
    
    Qtheta[end-T+1:end,end-T+1:end] = eye(T)
    return Qtheta
end

function tmp_inst_A_scale_new(n,Ridx,T,line,therm_a,int_length)
    """ Augment A with T additional rows relating
    angle difference variables to angle variables.
    
    Returns a T-by-(n+nr+2)*T matrix that may be
    concatenated with the output of temp_inst_A
    
    Arguments:    
    * n is the number of nodes in the network
    * Ridx is a vector indicating wind nodes
    * T is the number of time steps
    * tau is the thermal coefficient from IEEE 738
    * slack is the index of the slack bus
    * line is the pair (i,k) indicating the chosen
    line
    * therm_a is a constant defined by heating parameters
    * int_length is interval length in seconds (e.g. 600 
        for 10 minutes)
    """
    (i,k) = line
    nr = length(Ridx)
    
    A = zeros(T,(nr+n+2)*T)
    
    for t = 1:T
        i_pos = (nr+n+1)*(t-1) + nr + i
        k_pos = (nr+n+1)*(t-1) + nr + k
        coef = sqrt(-exp(therm_a*int_length)^(T-t+1) + exp(therm_a*int_length)^(T-t))
        A[t,i_pos] = -coef
        A[t,k_pos] = coef
        A[t,(n+nr+1)*T + t] = 1
    end
    return A
end

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
    idx = [idx1,idx2,idx3]
    
    A1 = A[:,idx1]
    A2 = A[:,idx2]
    A3 = A[:,idx3]
    return A1,A2,A3,idx1,idx2,idx3
end

function find_x_star(Qobj,A,b,T)
    """ x_star is the vector by which the problem must
    be translated in the first step of the temporal
    instanton QCQP solution.

    x_star is chosen to be the point in the set Ax=b
    nearest to the origin such that x_star[idx3] = 0.
    This condition ensures no linear term is introduced
    into the quadratic constraint.
    """
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
    dim_N = n - rank(A)
    # if A has full row rank:
    # dim_N = n - m
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

function solve_instanton_qcqp(G_of_x,Q_of_x,A,b,T)
    """ This function solves the following quadratically-
    constrained quadratic program:

        min  G_of_x
        s.t. A*x = b
             Q_of_x = 0

    where   G_of_x = x'*Qobj*x,
            Q_of_x = x'*Qtheta*x - c

    Thus, an equivalent problem expression is:

        min  z'*Qobj*z
        s.t. A*z = b
             z'*Qtheta*z = c

    The solution method is due in part to work by
    Dr. Dan Bienstock of Columbia University. It
    involves translating and rotating the problem,
    using partial KKT conditions, and solving the
    resulting secular equation.
    """

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

function solve_temporal_instanton(
    Ridx,
    Y,
    G0,
    P0,
    D0,
    Sb,
    ref,
    lines,
    res,
    reac,
    k,
    line_lengths,
    Tamb,
    T0,
    int_length)
    """ Convenience function used to perform temporal
    instanton analysis on many lines in a system at once.
    """

    n = length(k)
    nr = length(Ridx)
    T = int64(length(find(P0))/nr)
    
    numLines = length(lines)

    # Initialize progress meter:
    prog = Progress(length(find(line_lengths)),1)

    # Initialize vars used to store results:
    score = Float64[]
    α = Array(Vector{Float64},0)

    θ = Array(Array,0)
    x = Array(Array,0)
    diffs = Array(Array,0)
    xopt = Array(Array,0)

    # Create Qobj:
    Qobj = tmp_inst_Qobj(n,nr,T)
    # Augment Qobj with additional rows and columns of zeros:
    Qobj = tmp_inst_pad_Q(full(Qobj),T)

    # Create A1 (only A2 changes during opt.):
    A1 = full(tmp_inst_A(Ridx,T,Y,ref,k))
    A1 = [A1 zeros((n+1)*T,T)]

    # Create b:
    b = tmp_inst_b(n,T,G0,P0,D0)
    # Augment b with new elements:
    tmp_inst_pad_b(b,T)

    # Create Qtheta:
    Qtheta = tmp_inst_Qtheta(n,nr,T)

    # Form objective quadratic:
    G_of_x = (Qobj,0,0)    

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

        xvec,sol = solve_instanton_qcqp(G_of_x,Q_of_x,A,b,T)
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

        next!(prog)
    end
    return score,x,θ,α,diffs,xopt
end

end