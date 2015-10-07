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
    T::Int64
    )
    m,n = size(A)
    Qobj_orig = Qobj[1]

    opt = Array(Vector{Float64},0)

    # Partition A:
    A1,A2,idx1,idx2,idx3 = partition_A(A,Qobj[1],T)

    # Find translation point:
    x_star = find_x_star(A1,A2,idx1,idx2,n,b)

    # Translate quadratics:
    Qobj = translate_quadratic(Qobj,x_star)
    Qconstr = translate_quadratic(Qconstr,x_star)

    N = kernel_rotation(A) # take only cols spanning N(A)

    Qobj = rotate_quadratic(Qobj,N')
    Qconstr = rotate_quadratic(Qconstr,N')

    # testing to see if I can make the eig() argument smaller
    #
    # # partition N in the same way I partitioned A:
    # N3 = N[idx3,:]
    # d,v = eigs(N3'*N3; nev=T)
    # # There will be NR*n EVAs, but at most n will be nonzero.
    #
    #
    # Qconstr[1] ==
    ########

    D,U = eig(full(Qconstr[1]))
    # eigs won't work because nev cannot be size(Q,1):
    # D,U = eigs(Q_of_z[1],nev=size(Q_of_z[1],1))
    D = round(D,10)

    K = return_K(D)
    Kinv = diagm(1./diag(K))

    Qobj = rotate_quadratic(Qobj,(U*Kinv)')
    Qconstr = rotate_quadratic(Qconstr,(U*Kinv)')

    B11,B12,B21,B22,b1,b2 = partition_B(Qobj,Qconstr)

    Bhat,bhat = return_Bhat(B11,B12,B22,b1,b2)

    tinynumber = 1e-8
    w0 = find_w(0.0,Bhat,bhat/2)
    c = -Qconstr[3]
    if abs((w0'*w0) - c)[1] < tinynumber
        println("v=0 works!")
    end

    num = bhat/2
    poles = unique(round(diag(Bhat),10))
    if maxabs(num) > 1/tinynumber
        warn("No solution for secular equation; num elements too large")
        return [],NaN
    end

    solutions, vectors = solvesecular(num,poles,c)
    sol = zeros(length(vectors))
    for i in 1:length(vectors)
        w2 = vectors[i]
        xvec = return_xopt(w2,B11,B12,b1,N,U,K,x_star)
        sol[i] = (xvec'*Qobj_orig*xvec)[1]
        push!(opt,xvec)
    end
    if isempty(sol)
        warn("No solution for secular equation; solver returned []")
        return [], NaN
    else
        return opt[indmin(sol)],minimum(sol)
    end
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

    # why does all allocation happen here?
    # (parallel question)
    n = length(k)
    nr = length(Ridx)
    T = round(Int64,length(find(R0))/nr)
    numLines = length(lines)

    # Form objective quadratic:
    Qobj = tmp_inst_Qobj(n,nr,T,corr; pad=true)
    G_of_x = (Qobj,zeros(size(Qobj,1)),0.0)

    # Create A1 (only A2, the bottom part,
    # changes during line loop):
    A1 = tmp_inst_A1(Ridx,T,Y,ref,k; pad=true)

    b = tmp_inst_b(n,T,G0,R0,D0; pad=true)
    Qtheta = tmp_inst_Qconstr(n,nr,T)

    # Exclude lines with zero length or zero resistance:
    nz_line_idx = intersect(find(line_lengths),find(res))

    # truncate to go through subset of lines (testing only):
    if maxlines > 0 && maxlines < length(nz_line_idx)
        # nz_line_idx = nz_line_idx[1:maxlines]
        nz_line_idx = nz_line_idx[randperm(length(nz_line_idx))[1:maxlines]]
    end

    # loop through lines (having non-zero length)
    results = @parallel (vcat) for idx in nz_line_idx
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
        xvec,sol = solve_instanton_qcqp(G_of_x,Q_of_x,A,b,T)
        if isempty(xvec)
            xvec,sol = zeros(size(Qobj,1)),sol
        end

        # this is what will be concatenated into `results`:
        xvec,sol,toq()
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
    results::Array{Tuple{Array{Float64,1},Float64,Float64},1},
    n::Int64,
    nr::Int64,
    T::Int64;
    return_as_type::Bool = true
    )
    # Store results in human-readable form:
    score   = Vector{Float64}()
    x       = Vector{Vector{Vector{Float64}}}()
    θ       = Vector{Vector{Vector{Float64}}}()
    α       = Vector{Vector{Float64}}()
    diffs   = Vector{Vector{Float64}}()
    xopt    = Vector{Vector{Float64}}()
    linetimes = Vector{Float64}()

    for i in 1:size(results,1)
        xvec,sol,linetime = results[i]
        # array of vectors with Float64 values:
        deviations = Array(Vector{Float64},0)
        angles = Array(Vector{Float64},0)
        alpha = Float64[]

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
        push!(linetimes,linetime)
    end
    return InstantonOutputData(score,x,θ,α,diffs,xopt,linetimes)
end
