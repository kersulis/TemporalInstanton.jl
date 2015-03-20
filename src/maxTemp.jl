module maxTemp

export
	maxTempModel,maxTempInstanton

####### IMPORTS: JULIA PACKAGES ########
using JuMP
using Ipopt

push!(LOAD_PATH, dirname(@__FILE__))
using TemporalInstanton

function maxTempInstanton(
    Ridx::Vector,
    Y::Array,
    ref::Int,
    G0::Vector,
    P0::Vector,
    D0::Vector,
    k::Vector,
    tau::Float64,
    lines::Array,
    c::Float64)
    
    n = size(Y,1)
    nr = length(Ridx)
    T = int(length(G0)/n) # infer number of time steps

    Qobj = tmp_inst_Qobj(n,nr,T)
    # Augment Qobj with additional rows and columns of zeros:
    Qobj = tmp_inst_pad_Q(full(Qobj),T)
    A1 = full(tmp_inst_A(Ridx,T,Y,ref,k))
    A1 = [A1 zeros((n+1)*T,T)]
    b = tmp_inst_b(n,T,G0,P0,D0)
    # Augment b with new elements:
    tmp_inst_pad_b(b,T)
    Qtheta = tmp_inst_Qtheta(n,nr,T)#,tau)

    score = Float64[]
    result = Symbol[]
    α = Array(Vector{Float64},0)

    θ = Array(Array,0)
    x = Array(Array,0)
    diffs = Array(Array,0)

    TT = STDOUT # save original STDOUT
    redirect_stdout()

    # Loop through all lines:
    for line in lines
        # array of vec. with Float values:
        deviations = Array(Vector{Float64},0)
        angles = Array(Vector{Float64},0)
        alpha = Float64[]

        # Create instance of instanton problem
        A2 = tmp_inst_A_scale(n,Ridx,T,tau,ref,line)
        # Augment A with new rows:
        A = [A1; A2]

        # Create solver model
        m,z,LinearConstrs = maxTempModel(Qobj,A,b,Qtheta,P0,c)

        # Solve
        status = solve(m)

        # Store results
        xvec = getValue(z)[:]
        push!(score, getObjectiveValue(m))
        push!(result,status)

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
        push!(x,deviations)
        push!(θ,angles)
        push!(α,alpha)
        push!(diffs,xvec[end-T+1:end])
    end
    redirect_stdout(TT)
    return score,result,x,θ,α,diffs
end

function maxTempModel(Qobj,A,b,Qtheta,P0,c)
    # Build a model:
    m = Model(solver = IpoptSolver()) # Use MOSEK to solve model
    numRows,numVars = size(A)
    Nr = length(find(P0))
    @defVar(m,z[1:numVars])
    setObjective(m,:Max,sum([(Qtheta[i,i]*z[i]^2) for i in 1:size(Qtheta,1)]))

    @defConstrRef LinearConstrs[1:numRows]
    for i in 1:numRows
        LinearConstrs[i] = @addConstraint(m, sum([A[i,k]*z[k] for k in 1:numVars]) == b[i])
    end

    windPos = find(diag(Qobj))
    windVal = P0[find(P0)]
    @defConstrRef windLims[1:length(windPos)]
    for i = 1:length(windPos)
        windLims[i] = @addConstraint(m, -c*windVal[i] <= z[windPos[i]] <= c*windVal[i])
    end
    return m,z,LinearConstrs
end

end