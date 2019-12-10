using JLD
using SparseArrays: spzeros, sparse
using Distributed: @distributed
using LinearAlgebra: diag, svd, qr
using Random: randperm

"""
    solve_instanton_qcqp(Qobj, Qconstr, A, b, T) -> (optVec, minObj)
Solve the following quadratically-
constrained quadratic program:

    min  G_of_x
    s.t. A*x = b
    Q_of_x = 0

where `G_of_x = x' * Qobj * x` and `Q_of_x = x' * Qtheta * x - c`.

Thus, an equivalent problem expression is:

    min  z' * Qobj * z
    s.t. A * z = b
         z' * Qconstr * z = c

The solution method is due in part to Dr. Dan
Bienstock of Columbia University. It involves
translating and rotating the problem, applying
partial KKT conditions, and solving the
resulting secular equation.

*Return `NaN` if there is no intersection
between the secular equation and horizontal line.*
"""
function solve_instanton_qcqp(
    Qobj::Tuple{SparseMatrixCSC{Float64, Int64},Vector{Float64}, Float64},
    Qconstr::Tuple{SparseMatrixCSC{Float64, Int64}, Vector{Float64}, Float64},
    A::SparseMatrixCSC{Float64, Int64},
    b::Vector{Float64},
    T::Int64,
    N1::SparseMatrixCSC{Float64, Int64}
    )

    m, n = size(A)
    objidx = findall(diag(Qobj[1]) .> 0)
    Qobj_orig = Qobj[1][objidx, objidx]

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

    tTrans = time()
    A1, A2, idx1, idx2 = partition_A(A, Qobj[1], T)
    x_star = translation_point(A1, A2, idx1, idx2, n, b)

    # nb = Int64(m/T - 2)
    # nd = Int64(size(Qobj_orig,1)/T)
    # nt = T

    ## testing: add a column of the null basis of A. Shouldn't change solution.
    # nullbasis = full(kernel_basis(A[:,1:end-6]))
    # x_star[1:end-6] += nullbasis*rand(size(nullbasis,2))

    Qobj = translate_quadratic(Qobj, x_star)

    # (Qconstr is invariant under translation by x_star.)
    tTrans = time() - tTrans

    tKern = time()
    # most of null basis N has been computed and was passed
    # in as N1. The bottom part, N2, is:
    N2 = -A[(end - T + 1):end, 1:(end - T)] * N1

    # stack N1 and N2 to obtain null basis N
    N = [N1; N2]::SparseMatrixCSC{Float64, Int64}
    tKern = time() - tKern

    tEig = time()
    # svds method: faster than eig, svd, and eigs
    # U,Sn,Vn = svds(N2',nsv=T)
    # USnVn = svds(N2'; nsv=T)[1]
    # U, Sn = USnVn.U, USnVn.S

    # Julia 1.0: won't use svds on small matrices
    USnVn = svd(Array(N2'))
    U = USnVn.U[:, 1:T]
    Sn = USnVn.S[1:T]

    # augment with zeros, take qr:
    U = qr([sparse(round.(U; digits=10)) spzeros(n - m, n - m - T)])

    # Julia 1.0: it's easier to extract Q now, but beware permutation.
    # Note: MUST PERMUTE ROWS HERE to cover any case where T > 1
    U = sparse(U.Q[invperm(U.prow), :])
    # U = sparse(U.Q) # ROWS ARE OUT OF ORDER IF YOU ONLY DO THIS

    # U = sparse(round(U,13)) # round sparse is expensive!
    K = [Sn; ones(n - m - T)]
    Kinv = 1 ./ K

    D = [ones(T); zeros(size(U, 1) - T)]
    tEig = time() - tEig

    tMap = time()
    map_mat = (N * (U .* Kinv')) # ::SparseMatrixCSC{Float64, Int64}
    tMap = time() - tMap

    tRot = time()
    Qobj = rotate_quadratic(Qobj, map_mat')
    tRot = time() - tRot

    B11, B12, B22, b1, b2, i1, i2 = partition_B(Qobj, D)
    # temp = full(B12')/Symmetric(full(B11))
    # B11 = lufact(full(B11))
    # temp = B12'/B11
    # save("../data/B.jld","B11",B11)
    # temp = full(B12')/Symmetric(full(B11))

    temp = Array(B12') / Array(B11)
    Bhat = B22 - temp * B12
    bhat = b2 - temp * b1

    # works, but still need to find bhat somehow:
    # Bhat = schurcomp(full(Qobj[1]), i1, i2)

    tSec = time()
    # solve secular equation
    num = bhat / 2
    poles = unique(round.(diag(Bhat); digits=10))
    c = -Qconstr[3]

    v, w2 = solvesecular(num, poles, c)

    if any(w2 .== Inf)
        xvec = fill(Inf, length(x_star))
        objVal = Inf
        times = [tTrans; tKern; tEig; tRot; Inf; Inf]
        return xvec, objVal, times
    end

    tSec = time() - tSec

    # map back to original problem space
    tTemp = time()
    # w1opt = -B11\(B12*w2 + b1/2)
    # re-use temp once more here (saves time)
    w1opt = -temp' * w2 - B11 \ (b1 / 2)
    wopt = zeros(n - m)
    wopt[i1] = w1opt
    wopt[i2] = w2

    xvec = map_mat * wopt + x_star
    xobjonly = xvec[objidx]
    objVal = (xobjonly' * (Qobj_orig * xobjonly))[1]
    tMap += time() - tTemp

    # concatenate times
    times = [tTrans; tKern; tEig; tRot; tSec; tMap]

    return xvec, objVal, times
end

"""
Perform temporal instanton analysis on many lines in a system at once.
If `@addprocs(n)` has been run, computation is parallelized by
splitting the vector of lines into n equal-length parts.

Return an instance of `InstantonOutput`.
"""
function solve_temporal_instanton(
    Ridx::Vector{Int64},
    Y::SparseMatrixCSC{Float64, Int64},
    G0::Vector{Float64},
    R0::Vector{Float64},
    D0::Vector{Float64},
    Sb::Float64,
    ref::Int64,
    lines::Vector{Tuple{Int64, Int64}},
    res::Vector{Float64},
    reac::Vector{Float64},
    k::Vector{Float64},
    line_lengths::Vector{Float64},
    line_conductors::Vector{ACSRSpecsMetric},
    current_limits::Vector{Float64},
    Tamb::Float64,
    T0::Vector{Float64},
    Tlim::Vector{Float64},
    time_values::StepRangeLen{Int64,Int64,Int64},
    auto_prec::Float64=0.0,
    analytic_lines::Vector{Int64}=collect(1:length(lines));
    maxlines::Int64=0,
    silent::Bool=false,
    analyze_overheated_only::Bool=false
    )

    # aggregate timing results:
    timeVecs = Vector{Vector{Float64}}()

    tBuild = time()

    n = length(k)
    nr = length(Ridx)
    T = length(time_values) - 1
    int_length = time_values.step
    numLines = length(lines)

    ##########################################
    ##           Matrix Formation           ##
    ##########################################
    # Form objective quadratic:
    # Qobj = tmp_inst_Qobj(n, nr, T, corr; pad=true)
    Qobj = tmp_inst_Qobj(n, nr, T, auto_prec)
    G_of_x = (Qobj, zeros(size(Qobj, 1)), 0.0)

    # Create A1 (only A2, the bottom part,
    # changes during line loop):
    A1 = tmp_inst_A1(Ridx, T, Y, ref, k; pad=true)

    b = tmp_inst_b(n, T, G0, R0, D0; pad=true)
    Qtheta = tmp_inst_Qconstr(n, nr, T)

    ##########################################
    ##           r = 0 Check                ##
    ##########################################
    # Exclude lines with zero length or zero resistance:
    analytic_lines = intersect(analytic_lines, findall(line_lengths .> 0), findall(res .> 0))
    if length(findall(res .== 0)) != 0 && !silent
        @info "r=0 check: \tremoving $(length(findall(res .== 0))) lines"
    end

    ##########################################
    ##        Shift Factor Check            ##
    ##########################################
    # rule out lines where all renewable shift factors are 0:
    ISF = isf(Y, lines, ref, k)[:, Ridx]
    small = 1e-8
    singular_line_idx = findall(maximum(abs.(ISF); dims=2)[:] .< small)
    if !isempty(singular_line_idx) && !silent
        @info "ISF pre-check: \tremoving $(length(singular_line_idx)) lines"
    end
    analytic_lines = setdiff(analytic_lines, singular_line_idx)

    ##########################################
    ##    Analyze overheated lines only     ##
    ##########################################
    if analyze_overheated_only
        underheated_lines = findall(T0 .< Tlim)
        analytic_lines = setdiff(analytic_lines, underheated_lines)
    end

    if isempty(analytic_lines)
        @warn "No lines to analyze."
        return nothing
    end

    # truncate to go through subset of remaining lines (testing only):
    if maxlines > 0 && maxlines < length(analytic_lines)
        analytic_lines = analytic_lines[randperm(length(analytic_lines))[1:maxlines]]
    end

    tBuild = time() - tBuild
    # @show tBuild

    tBasis = time()
    # Find top part of null basis using A1.
    # reuse for remaining lines.

    # JLD.save("basis.jld", "A1", A1[:, 1:(end - T)])
    N1 = kernel_basis(A1[:, 1:(end - T)])
    tBasis = time() - tBasis
    # @show tBasis

    push!(timeVecs, [tBuild])
    ##########################################
    ##        Begin Line Loop               ##
    ##########################################
    # loop through lines (having non-zero length)
    results = @distributed (vcat) for idx in analytic_lines
        line_time = time()
        line = lines[idx]
        # if isa(line_conductors, Vector{String})
            # conductor_name = line_conductors[idx]
            # conductor_params = return_conductor_params(conductor_name)
            # compute line_params based on current line:
            # line_params = LineParams(line[1], line[2], res[idx], reac[idx], line_lengths[idx])
            # Tinit = T0
        # else
            # line_params = line_conductors[1][idx]
            Ilim = current_limits[idx]
            conductor_params = return_conductor_params(line_conductors[idx], Ilim, Tamb, T0[idx], Tlim[idx])
            line_params = LineParams(line[1], line[2], res[idx], reac[idx], line_lengths[idx])
            Tinit = T0[idx]
        # end

        t_a, t_c, t_d, t_f = return_thermal_constants(line_params, conductor_params, Tamb, Sb, int_length, T, Tinit)

        # thermal constraint, Q(z) = 0:
        kQtheta = (t_a / t_c) * (conductor_params.Tlim - t_f)
        Q_of_x = (Qtheta, zeros(size(Qtheta, 1)), kQtheta)

        # Create A2 based on chosen line:
        A2 = tmp_inst_A2(n, Ridx, T, line, t_a, int_length)

        # Stack A1 and A2:
        A = [A1; A2] # ::SparseMatrixCSC{Float64, Int64}

        # Computationally expensive part: solving QCQP
        xvec, sol, times = solve_instanton_qcqp(G_of_x, Q_of_x, A, b, T, N1)
        if isempty(xvec)
            xvec, sol = zeros(size(Qobj, 1)), sol
        end

        line_time = time() - line_time
        # this is what will be concatenated into `results`:
        # xvec, (sol, idx), toq(), times
        xvec, (sol, idx), line_time, times
    end
    if length(analytic_lines) == 1
        results = [results]
    end
    process_instanton_results(results, analytic_lines, n, nr, T)
end

"""
    solve_temporal_instanton(d, [maxlines=0]) -> results
Convenience method for performing temporal instanton analysis
on an instance of `InstantonInput`.
"""
solve_temporal_instanton(
    d::InstantonInput,
    analytic_lines::Vector{Int64}=collect(1:length(d.lines));
    maxlines::Int64=0,
    silent::Bool=false,
    analyze_overheated_only::Bool=false) = solve_temporal_instanton(
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
    d.current_limits,
    d.Tamb,
    d.T0,
    d.Tlim,
    d.time_values,
    d.auto_prec,
    analytic_lines;
    maxlines=maxlines,
    silent=silent,
    analyze_overheated_only=analyze_overheated_only
)

"""
    process_instanton_results(results, n, nr, T, [return_as_bool=true]) -> o
Accepts `results` output from `solve_temporal_instanton`,
returns output data in more human-readable form as instance of
InstantonOutput.

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
    analytic_lines,
    n::Int64,
    nr::Int64,
    T::Int64;
    return_as_type::Bool=true
    )
    # Store results in human-readable form:
    score   = Tuple{Float64, Int64}[]
    x       = Vector{Vector{Float64}}[]
    θ       = Vector{Vector{Float64}}[]
    α       = Vector{Float64}[]
    diffs   = Vector{Float64}[]
    xopt    = Vector{Float64}[]
    linetimes = Float64[]
    timeVecs = [ri[4] for ri in results]

    # JLD.save("timing.jld", "timeVecs", timeVecs)

    for i in 1:size(results, 1)
        xvec, sol, linetime = results[i][1:3]
        # array of vectors with Float64 values:
        deviations = Vector{Float64}[]
        angles = Vector{Float64}[]
        alpha = Float64[]

        push!(score, sol)
        if isinf(sol[1])
            push!(deviations, [])
            push!(angles, [])
            push!(alpha, NaN)
            push!(diffs, [])
        else
            # Variable breakdown:
            # (nr + n + 1) per time step
            #   First nr are deviations
            #   Next n are angles
            #   Last is mismatch
            # T variables at the end: angle diffs
            for t in 1:T
                push!(deviations, xvec[((nr + n + 1) * (t - 1) + 1):((nr + n + 1) * (t - 1) + nr)])
                push!(angles, xvec[((nr + n + 1) * (t - 1) + nr + 1):((nr + n + 1) * (t - 1) + nr + n)])
                push!(alpha, xvec[(nr + n + 1) * t])
            end
            push!(diffs, xvec[(end - T + 1):end])
        end

        push!(x, deviations)
        push!(θ, angles)
        push!(α, alpha)
        push!(xopt, xvec)
        push!(linetimes, linetime)
    end
    p = sortperm([s[1] for s in score])
    return InstantonOutput(score[p], analytic_lines[p], x[p], θ[p], α[p], diffs[p], xopt[p], linetimes[p])
end
