using JLD

"""
    solve_instanton_qcqp(Qobj,Qconstr,A,b,T) -> (optVec,minObj)
Solve the following quadratically-
constrained quadratic program:

    min  G_of_x
    s.t. A*x = b
    Q_of_x = 0

where
    G_of_x = x'\*Qobj\*x,
    Q_of_x = x'\*Qtheta\*x - c

Thus, an equivalent problem expression is:

    min  z'\*Qobj\*z
    s.t. A\*z = b
         z'\*Qconstr\*z = c

The solution method is due in part to Dr. Dan
Bienstock of Columbia University. It involves
translating and rotating the problem, applying
partial KKT conditions, and solving the
resulting secular equation.

*Return `NaN` if there is no intersection
between the secular equation and horizontal line.*
"""
function solve_instanton_qcqp(
    Qobj::Tuple{SparseMatrixCSC{Float64,Int64},Vector{Float64},Float64},
    Qconstr::Tuple{SparseMatrixCSC{Float64,Int64},Vector{Float64},Float64},
    A::SparseMatrixCSC{Float64,Int64},
    b::Vector{Float64},
    T::Int64,
    x_star::Vector{Float64}
    )
    m,n = size(A)
    Qobj_orig = Qobj[1]

    ###### October 13, 2015
    # Timing analysis
    # keep track of times:
    # tTrans    = [ti[1] for ti in t]
    # tKern     = [ti[2] for ti in t]
    # tEig      = [ti[3] for ti in t]
    # tRot      = [ti[4] for ti in t]
    # tSec      = [ti[5] for ti in t]
    # tMap      = [ti[6] for ti in t]
    ######

    opt = Array(Vector{Float64},0)

    # Partition A:
    A1,A2,idx1,idx2,idx3 = partition_A(A,Qobj[1],T)

    tic()
    # Find translation point:
    x_star = find_x_star(A1,A2,idx1,idx2,n,b)

    # Translate quadratics:
    Qobj = translate_quadratic(Qobj,x_star)
    Qconstr = translate_quadratic(Qconstr,x_star)
    tTrans = toq()

    # store for Eig:
    Qc = Qconstr[1]

    tic()
    N = kernel_rotation(A) # take only cols spanning N(A)
    tKern = toq()

    tic()
    Qobj = rotate_quadratic(Qobj,N')
    Qconstr = rotate_quadratic(Qconstr,N')
    tRot = toq()

    tic()
    D,U = eig(full(Qconstr[1]))
    # eigs won't work because nev cannot be size(Q,1):
    # D,U = eigs(Q_of_z[1],nev=size(Q_of_z[1],1))
    # save("Qeig.jld","D",D,"U",U,"Qc",Qc,"N",N,"A",A,"Qconstr",Qconstr[1])
    D = round(D,10)

    K = return_K(D)
    Kinv = diagm(1./diag(K))

    Qobj = rotate_quadratic(Qobj,(U*Kinv)')
    Qconstr = rotate_quadratic(Qconstr,(U*Kinv)')

    tEig = toq()

    ####### November 3, 2015
    # Alternative to computing Eig:
    # re-use existing decomposition.
    # replaces section commented out above.
    #
    # tic()
    # # D,U = eig(full(Qconstr[1]))
    # # D = round(D,10)
    # D = diag(Qc)[m+1:n]
    # U = (N')[:,m+1:n]
    #
    # K = return_K(D)
    # Kinv = diagm(1./diag(K))
    #
    # Qobj = rotate_quadratic(Qobj,(U*Kinv)')
    # Qconstr = rotate_quadratic(Qconstr,(U*Kinv)')
    # save("Qconstr.jld","Qconstr",Qconstr)
    # tEig = toq()
    ####### End alternative Eig

    B11,B12,B21,B22,b1,b2 = partition_B(Qobj,Qconstr)
    Bhat,bhat = return_Bhat(B11,B12,B22,b1,b2)

    # tinynumber = 1e-8
    # w0 = find_w(0.0,Bhat,bhat/2)
    # c = -Qconstr[3]
    # if abs((w0'*w0) - c)[1] < tinynumber
    #     println("v=0 works!")
    # end

    num = bhat/2
    poles = unique(round(diag(Bhat),10))
    c = -Qconstr[3]
    # if maxabs(num) > 1/tinynumber
    #     warn("No solution for secular equation; num elements too large")
    #     return [],NaN
    # end

    # solve secular equation
    tic()
    v, w2 = solvesecular(num,poles,c)
    tSec = toq()

    # map back to original problem space
    tic()
    # used to map w back to x space:
    map_mat = N*sparse(U)*spdiagm(1./diag(K))

    w1opt = -B11\(B12*w2 + b1/2)
    wopt = [w1opt;w2]
    xvec = map_mat*wopt + x_star

    objVal = (xvec'*Qobj_orig*xvec)[1]
    tMap = toq()

    # concatenate times
    times = [tTrans;tKern;tEig;tRot;tSec;tMap]

    # return xvec, sol, times
    return xvec,objVal,times
end

"""
Perform temporal instanton analysis on many lines in a system at once. If `@addprocs(n)` has been run, computation is parallelized by splitting the vector of lines into n equal-length parts.

Return `results`, a vector of tuples that may be processed into a more accessible form by `process_instanton_results`.
"""
function solve_temporal_instanton(
    Ridx::Vector{Int64},
    Y::SparseMatrixCSC{Float64,Int64},
    G0::Vector{Float64},
    R0::Vector{Float64},
    D0::Vector{Float64},
    Sb::Float64,
    ref::Int64,
    lines::Array{Tuple{Int64,Int64},1},
    res::Vector{Float64},
    reac::Vector{Float64},
    k::Vector{Float64},
    line_lengths::Vector{Float64},
    line_conductors::Vector{ASCIIString},
    Tamb::Float64,
    T0::Float64,
    int_length::Float64,
    corr::Array{Float64,2} = Array{Float64,2}();
    maxlines::Int64 = 0
    )

    # aggregate timing results:
    timeVecs = Vector{Vector{Float64}}()

    tic()

    n = length(k)
    nr = length(Ridx)
    T = round(Int64,length(find(R0))/nr)
    numLines = length(lines)

    ##########################################
    ##           Matrix Formation           ##
    ##########################################
    # Form objective quadratic:
    Qobj = tmp_inst_Qobj(n,nr,T,corr; pad=true)
    G_of_x = (Qobj,zeros(size(Qobj,1)),0.0)

    # Create A1 (only A2, the bottom part,
    # changes during line loop):
    A1 = tmp_inst_A1(Ridx,T,Y,ref,k; pad=true)

    b = tmp_inst_b(n,T,G0,R0,D0; pad=true)
    Qtheta = tmp_inst_Qconstr(n,nr,T)

    ##########################################
    ##           r = 0 Check                ##
    ##########################################
    # Exclude lines with zero length or zero resistance:
    analytic_lines = intersect(find(line_lengths),find(res))
    println("r=0 check: \tremoving $(length(find(res.==0))) lines")

    ##########################################
    ##        Shift Factor Check            ##
    ##########################################
    # rule out lines where all renewable shift factors are 0:
    ISF = isf(Y,lines,ref,k)[:,Ridx]
    small = 1e-8
    singular_line_idx = find(1 - [maxabs(ISF[i,:])>small for i in 1:size(ISF,1)])
    println("ISF pre-check: \tremoving lines $singular_line_idx")
    analytic_lines = setdiff(analytic_lines,singular_line_idx)

    # truncate to go through subset of remaining lines (testing only):
    if maxlines > 0 && maxlines < length(analytic_lines)
        analytic_lines = analytic_lines[randperm(length(analytic_lines))[1:maxlines]]
    end
    tBuild = toq()
    push!(timeVecs,[tBuild])

    ###########################################
    # 2015-11-03
    # We can find a translation point x_star
    # that works for each QCQP.

    # compute x_star that will work for all lines:
    # x_star = A1\(b[1:end-T])
    # save("xstar.jld","x_star",x_star)
    ###########################################

    ##########################################
    ##        Begin Line Loop               ##
    ##########################################
    # loop through lines (having non-zero length)
    results = @parallel (vcat) for idx in analytic_lines
        tic()
        line = lines[idx]
        conductor_name = line_conductors[idx]
        conductor_params = return_conductor_params(conductor_name)
        # if current line uses a different conductor than previous,
        # re-compute conductor_params:
        # if line_conductors[idx] != conductor_name
        #     conductor_name = line_conductors[idx]
        #     conductor_params = return_conductor_params(conductor_name)
        # end
        # compute line_params based on current line:
        line_params = LineParams(line[1],line[2],res[idx],reac[idx],line_lengths[idx])

        (therm_a,therm_c,therm_d,therm_f) = return_thermal_constants(line_params,conductor_params,Tamb,Sb,int_length,T,T0)

        # thermal constraint, Q(z) = 0:
        kQtheta = (therm_a/therm_c)*(conductor_params.Tlim - therm_f)
        Q_of_x = (Qtheta,zeros(size(Qtheta,1)),kQtheta)

        # Create A2 based on chosen line:
        A2 = tmp_inst_A2(n,Ridx,T,line,therm_a,int_length)
        # Stack A1 and A2:
        A = [A1; A2]::SparseMatrixCSC{Float64,Int64}

        # Computationally expensive part: solving QCQP
        xvec,sol,times = solve_instanton_qcqp(G_of_x,Q_of_x,A,b,T,x_star)
        # push!(timeVecs,times)
        if isempty(xvec)
            xvec,sol = zeros(size(Qobj,1)),sol
        end

        # this is what will be concatenated into `results`:
        xvec,(sol,idx),toq(),times
    end
    return results
end

"""
    solve_temporal_instanton(d, [maxlines=0]) -> results
Convenience method for performing temporal instanton analysis on an instance of `InstantonInputData`.
"""
solve_temporal_instanton(d::InstantonInputData,maxlines::Int64 = 0) = solve_temporal_instanton(
    d.Ridx,
    d.Y,
    d.G0,
    d.R0,
    d.D0,
    d.Sb,
    d.ref,
    d.lines,
    d.res,
    d.reac,
    d.k,
    d.line_lengths,
    d.line_conductors,
    d.Tamb,
    d.T0,
    d.int_length,
    d.corr,
    maxlines = maxlines
)

"""
    process_instanton_results(results,n,nr,T,[return_as_bool=true]) -> o
Accepts `results` output from `solve_temporal_instanton`, returns output data in more human-readable form as instance of InstantonOutputData.

* `n` is the number of nodes in the network
* `nr` is the number of variable (renewable) generation nodes in the network
* `T` is the number of time steps in the analysis
"""
function process_instanton_results(
    results::Vector{Tuple{
        Vector{Float64},
        Tuple{Float64,Int64},
        Float64,
        Vector{Float64}
        }},
    n::Int64,
    nr::Int64,
    T::Int64;
    return_as_type::Bool = true
    )
    # Store results in human-readable form:
    score   = Vector{Tuple{Float64,Int64}}()
    x       = Vector{Vector{Vector{Float64}}}()
    θ       = Vector{Vector{Vector{Float64}}}()
    α       = Vector{Vector{Float64}}()
    diffs   = Vector{Vector{Float64}}()
    xopt    = Vector{Vector{Float64}}()
    linetimes = Vector{Float64}()
    timeVecs = [ri[4] for ri in results]
    save("../data/timing.jld","timeVecs",timeVecs)

    for i in 1:size(results,1)
        xvec,sol,linetime = results[i][1:3]
        # array of vectors with Float64 values:
        deviations = Array(Vector{Float64},0)
        angles = Array(Vector{Float64},0)
        alpha = Float64[]

        push!(score,sol)
        if isinf(sol[1])
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
        push!(linetimes,linetime)
    end
    return InstantonOutputData(score,x,θ,α,diffs,xopt,linetimes)
end
