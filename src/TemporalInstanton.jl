module TemporalInstanton

export tmp_inst_Qobj, tmp_inst_A, tmp_inst_b, tmp_inst_Qtheta, tmp_inst_A_scale, tmp_inst_pad_b, tmp_inst_pad_Q, temporalInstanton

function tmp_inst_Qobj(n,nr,T)
    """ Generate the objective function matrix
    Qobj from the problem dimensions.
    Assume no correlation between wind sites.
    """
    Qobj = sparse(diagm(repeat([ones(nr),zeros(n)],outer=[T])))
    Qobj = tmp_inst_pad_Q(Qobj,T)
    return Qobj
end

function tmp_inst_A(n,Ridx,T,Y,slack,k,tau,line)
    """ Generate the power balance constraint A matrix
    from problem dimensions, admittance matrix,
    and generator participation factors.
    Assumes the admittance matrix is n-by-n.
    
    Returns A, which is (n+1)*T-by-(n+nr)*T
    
    * n is the number of nodes in the network
    * Ridx is a vector indicating wind nodes
    * T is the number of time steps
    * Y is the admittance matrix (n-by-n)
    * slack is the index of the slack bus
    * k is the vector of generator participation factors
    """
    
    # zero out slack row of Y:
    adm = deepcopy(Y)
    adm[slack,:] = zeros(n)
    
    # remove slack column of Y:
    adm = adm[:,setdiff(1:n,slack)]
    
    # A has a block diagonal pattern where each
    # block is Atemp:
    Atemp = zeros(n+1,2*n)
    Atemp = [[-eye(n) adm -k]; ones(1,n) zeros(1,n-1) 1]
    
    # Remove columns corresponding to non-wind nodes:
    Atemp = sparse(Atemp[:,[Ridx,n+1:2*n]])
    
    # Now we can tile the Atemp matrix to generate A:
    A = Atemp
    for t = 2:T
        A = blkdiag(A, Atemp)
    end

    # These rows relate auxiliary angle variables to
    # original angle variables:
    A2 = tmp_inst_A_scale(n,Ridx,T,tau,slack,line)

    # Augment A with new rows:
    A = [[full(A) zeros((n+1)*T,T)]; A2]
    
    return A
end

function tmp_inst_b(n,T,G0,P0,D)
    """ Generate the vector b of power balance constraints.
    Assumes G0 and D are nT-by-1 vectors.
    """
    b = FloatingPoint[]
    
    # b = zeros((n+1)*T,1)
    netGen = G0 + P0 - D
    for t = 1:T
        start = (t-1)*n + 1
        stop = start + n - 1
        mismatch = -sum(netGen[start:stop])
        append!(b,netGen[start:stop])
        push!(b,mismatch)
    end

    # Extend b with T additional zeroes:
    tmp_inst_pad_b(b,T)

    return b
end

function tmp_inst_Qtheta(n,nr,T,tau)
    """ Generate Q_theta in the temperature constraint
    of a temporal instanton problem instance.
    "line" has the form (i,k), where i and k refer to
    the endpoints of the chosen line.
    """
    Qtheta = zeros((nr+n+1)*T,(nr+n+1)*T)
    
    Qtheta[end-T+1:end,end-T+1:end] = eye(T)
    return sparse(Qtheta)
end

function tmp_inst_A_scale(n,Ridx,T,tau,slack,line)
    """ Augment A with T additional rows relating
    angle difference variables to angle variables.
    
    Returns a T-by-(n+nr+1)*T matrix that may be
    concatenated with the output of temp_inst_A
    
    Arguments:    
    * n is the number of nodes in the network
    * Ridx is a vector indicating wind nodes
    * T is the number of time steps
    * tau is the thermal coefficient from IEEE 738
    * slack is the index of the slack bus
    * line is the pair (i,k) indicating the chosen
    line
    """
    (i,k) = line
    nr = length(Ridx)
    
    A = zeros(T,(nr+n+2)*T)
    
    for t = 1:T
        i_pos = (nr+n+1)*(t-1) + nr + i
        k_pos = (nr+n+1)*(t-1) + nr + k
        coef = tau^((1/2)*(T-t))
        A[t,i_pos] = -coef
        A[t,k_pos] = coef
        A[t,(n+nr+1)*T + t] = 1
    end
    slack_cols = [(nr+n+1)*(t-1) + nr + slack for t in 1:T]
    
    # remove slack columns:
    return sparse(A[:,setdiff(1:(n+nr+2)*T,slack_cols)])
end

function tmp_inst_pad_b(b,T)
    """ Append T zeros to vector b
    """
    append!(b,zeros(T))
end

function tmp_inst_pad_Q(Q,T)
    """ Add T rows and T columns of zeros
    to matrix Q
    """
    m,n = size(Q)
    return [[Q zeros(m,T)]; zeros(T,n+T)]
end

function temporalInstanton(Ridx,Y,slack,k,tau,line,G0,P0,D)
	""" Return an instance of the temporal instanton
	problem.
	"""
    n = size(Y,1)
    T = int(length(G0)/n)
    nr = length(Ridx)
    Qobj = tmp_inst_Qobj(n,nr,T)
    A = tmp_inst_A(n,Ridx,T,Y,slack,k,tau,line)
    b = tmp_inst_b(n,T,G0,P0,D)
    Qtheta = tmp_inst_Qtheta(n,nr,T,tau)
    return Qobj, A, b, Qtheta
end

end