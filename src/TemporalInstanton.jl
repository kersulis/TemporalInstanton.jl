module TemporalInstanton

export tmp_inst_Qobj, tmp_inst_A, tmp_inst_b, tmp_inst_Qtheta

function tmp_inst_Qobj(n,nr,T)
    """ Generate the objective function matrix
    Qobj from the problem dimensions.
    Assume no correlation between wind sites.
    """
    return sparse(diagm(repeat([ones(nr),zeros(n)],outer=[T])))
end

function tmp_inst_A(n,Ridx,T,Y,slack,k)
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
    
    return A
end

function tmp_inst_b(n,T,G0,D)
    """ Generate the vector b of power balance constraints.
    Assumes G0 and D are nT-by-1 vectors.
    """
    b = FloatingPoint[]
    
    # b = zeros((n+1)*T,1)
    netGen = G0 - D
    for t = 1:T
        start = (t-1)*n + 1
        stop = start + n - 1
        mismatch = -sum(netGen[start:stop])
        append!(b,netGen[start:stop])
        push!(b,mismatch)
    end
    return sparse(b)
end

function tmp_inst_Qtheta(n,nr,T,tau,line,slack)
    """ Generate Q_theta in the temperature constraint
    of a temporal instanton problem instance.
    "line" has the form (i,k), where i and k refer to
    the endpoints of the chosen line.
    """
    Qtheta = zeros((n+nr)*T,(n+nr)*T)
    i,k = line
    ei = zeros(n,1)
    ei[i] = 1
    ek = zeros(n,1)
    ek[k] = 1
    
    Q0 = (ei - ek)*(ei - ek)'
    
    Q0 = Q0[:,setdiff(1:n,slack)]
    Q0 = Q0[setdiff(1:n,slack),:]
    
    for t = 1:T
        start = nr + 1 + (nr+n)*(t-1)
        stop = start + n - 2
        Qtheta[start:stop,start:stop] = tau^(T-t)*Q0
    end
    return sparse(Qtheta)
end

function temporalInstanton(n,t,tau,Y,k,G0,D,line)
	""" Return an instance of the temporal instanton
	problem.
	"""
    Qobj = tmp_inst_Qobj(n,T)
    A = tmp_inst_A(n,T,Y,k)
    b = tmp_inst_b(n,T,G0,D)
    Qtheta = tmp_inst_Qtheta(n,T,tau,line)
    return Qobj, A, b, Qtheta
end