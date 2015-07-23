""" Generate the objective function matrix
Qobj from the problem dimensions.
Assume no correlation between wind sites.
"""
function tmp_inst_Qobj(n,nr,T; pad=true)
    Qobj = sparse(diagm(repeat([ones(nr);zeros(n+1)],outer=[T])))
    if !pad
        return Qobj
    else
        # Add T rows and columns of zeros to Q:
        r,c = size(Qobj)
        return [[full(Qobj) zeros(r,T)]; zeros(T,c+T)]
    end
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
function tmp_inst_A1(Ridx,T,Y,ref,k; pad=true)
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
    Atemp = sparse(Atemp[:,[Ridx;n+1:2*n+1]])

    # Now we can tile the Atemp matrix to generate A:
    A = Atemp
    for t = 2:T
        A = blkdiag(A, Atemp)
    end

    if !pad
        return full(A)
    else
        # pad A with columns of zeros (rows will be added during line loop)
        return [full(A) zeros((n+1)*T,T)]
    end
end

""" Generate the vector b of power balance constraints.
Assumes G0 and D are nT-by-1 vectors.
"""
function tmp_inst_b(n,T,G0,P0,D; pad=true)
    b = FloatingPoint[]
    netGen = G0 + P0 - D

    for t = 1:T
        start = (t-1)*n + 1
        stop = start + n - 1
        append!(b,netGen[start:stop])
        push!(b,0.)
    end
    if pad
        # Not sure when I would ever not want to do this...
        append!(b,zeros(T))
    end
    return b
end

""" Generate Q_theta in the temperature constraint
of a temporal instanton problem instance.
"line" has the form (i,k), where i and k refer to
the endpoints of the chosen line.
"""
function tmp_inst_Qtheta(n,nr,T)#,tau)
    Qtheta = zeros((nr+n+2)*T,(nr+n+2)*T)
    Qtheta[end-T+1:end,end-T+1:end] = eye(T)
    return Qtheta
end

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
function tmp_inst_A_scale_new(n,Ridx,T,line,therm_a,int_length)
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
