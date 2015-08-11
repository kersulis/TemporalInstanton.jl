module TemporalInstanton

#using MAT, MatpowerCases, JLD#, IProfile, HDF5, ProgressMeter

export
    solve_instanton_qcqp, solve_temporal_instanton, LineParams,
    ConductorParams, load_rts96_data, load_polish_data, createY,
    process_instanton_results, InstantonInputData, InstantonOutputData,

    # temporary:
    tmp_inst_Qobj,tmp_inst_A1,tmp_inst_b,tmp_inst_Qtheta,
    return_conductor_params,return_thermal_constants,tmp_inst_A2,
    kernel_rotation,

    # power flow:
    expand_renewable_vector,fixed_wind_A,fixed_wind_b,return_angles,
    return_angle_diffs,

    # plot:
    temperatureTrajectory

type InstantonInputData
    Ridx
    Y
    G0
    D0
    R0
    Sb
    ref
    lines
    res
    reac
    k
    line_lengths
    line_conductors
    Tamb
    T0
    int_length
    time_values
    corr

    #InstantonInputData(Ridx, Y, G0, D0, R0, Sb, ref, lines, res, reac, k, line_lengths, line_conductors) = InstantonInputData(Ridx, Y, G0, D0, R0, Sb, ref, lines, res, reac, k, line_lengths, line_conductors, [], [], [])
end

type InstantonOutputData
    score
    x
    θ
    α
    diffs
    xopt
end
include("DataLoad.jl")
include("PowerFlow.jl")
include("ThermalModel.jl")
include("QCQPMatrixBuilding.jl")
include("manipulations.jl")
include("SolveSecular.jl")
include("plot.jl")

""" Solve the following quadratically-
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

The solution method is due in part to Dr. Dan
Bienstock of Columbia University. It involves
translating and rotating the problem, using
partial KKT conditions, and solving the
resulting secular equation.

* Return `NaN` if there is no intersection
between the secular equation and horizontal line. *
"""
function solve_instanton_qcqp(G_of_x,Q_of_x,A,b,T)
    m,n = size(A)
    Qobj = G_of_x[1]
    c = - Q_of_x[3]

    opt = Array(Vector{Float64},0)

    # Partition A:
    A1,A2,idx1,idx2,idx3 = partition_A(A,Qobj,T)

    # Find translation point:
    x_star = find_x_star(A1,A2,idx1,idx2,n,b)

    # Translate quadratics:
    G_of_y = translate_quadratic(G_of_x,x_star)
    Q_of_y = translate_quadratic(Q_of_x,x_star)

    N = kernel_rotation(A, spqr=true) # take only cols spanning N(A)

    G_of_z = rotate_quadratic(G_of_y,N')
    Q_of_z = rotate_quadratic(Q_of_y,N')

    D,U = eig(full(Q_of_z[1]))
    # won't work because nev cannot be size(Q,1):
    # D,U = eigs(Q_of_z[1],nev=size(Q_of_z[1],1))
    D = round(D,10)

    K = return_K(D)
    Kinv = diagm(1./diag(K))

    G_of_w = rotate_quadratic(G_of_z,(U*Kinv)')
    Q_of_w = rotate_quadratic(Q_of_z,(U*Kinv)')

    B11,B12,B21,B22,b1,b2 = partition_B(G_of_w,Q_of_w)

    Bhat,bhat = return_Bhat(B11,B12,B22,b1,b2)

    eps = 1e-8
    w0 = find_w(0,Bhat,bhat/2)

    if abs((w0'*w0) - c)[1] < eps
        println("v=0 works!")
    end

    solutions, vectors = solve_secular(Bhat,bhat/2,-Q_of_w[3])
    if isempty(solutions)
        return [],NaN
    else
        sol = zeros(length(vectors))
        for i in 1:length(vectors)
            w2 = vectors[i]
            xvec = return_xopt(w2,B11,B12,b1,N,U,K,x_star)
            sol[i] = (xvec'*Qobj*xvec)[1]
            push!(opt,xvec)
        end
        return opt[indmin(sol)],minimum(sol)
    end
end

""" Perform temporal instanton analysis on
many lines in a system at once.

Inputs:

* `Ridx`            Vector: indices of nodes that have wind farms
* `Y`               Admittance matrix
* `G0`              Conventional generation dispatch
* `P0`              Renewable generation forecast
* `D0`              Conventional demand
* `Sb`              System base voltage
* `ref`             Index of system angle reference bus
* `lines`           Vector: tuples (from,to) of lines to loop through
* `res`             Vector: pu resistance for all lines
* `reac`            Vector: pu reactance for all lines
* `k`               Vector: conventional generator participation factors
* `line_lengths`    Vector: line lengths in meters
* `line_conductors` Vector: strings (e.g. "waxwing") indicating conductor types
* `Tamb`            Ambient temperature
* `T0`              Initial line temperature (TODO: compute within)
* `int_length`      Length of each time interval in seconds

Interpreting `score`:

* `score[i]::Float64`: solution found for line `i`
* `score[i]::Inf`: power flow has no effect on temperature (zero resistance)
* `score[i]::NaN`: no secular equation intersection for line `i`
"""
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
    line_conductors,
    Tamb,
    T0,
    int_length,
    corr=[])

    # why does all allocation happen here?
    # (parallel question)
    n = length(k)
    nr = length(Ridx)
    T = round(Int64,length(find(P0))/nr)
    numLines = length(lines)

    # Form objective quadratic:
    Qobj = tmp_inst_Qobj(n,nr,T,corr; pad=true)
    G_of_x = (Qobj,0,0)

    # Create A1 (only A2, the bottom part,
    # changes during line loop):
    A1 = tmp_inst_A1(Ridx,T,Y,ref,k; pad=true)

    b = tmp_inst_b(n,T,G0,P0,D0; pad=true)
    Qtheta = tmp_inst_Qtheta(n,nr,T)

    # Exclude lines with zero length:
    nz_line_idx = find(line_lengths.!=0)

    # loop through lines (having non-zero length)
    results = @parallel (vcat) for idx in nz_line_idx
        if res[idx] == 0.
            # power flow cannot influence line temp, so skip
            xvec,sol = ([],Inf)
        else
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
            Q_of_x = (Qtheta,0,kQtheta)

            # Create A2 based on chosen line:
            A2 = tmp_inst_A2(n,Ridx,T,line,therm_a,int_length)
            # Stack A1 and A2:
            A = [A1; A2]

            # Computationally expensive part: solving QCQP
            #try
            xvec,sol = solve_instanton_qcqp(G_of_x,Q_of_x,A,b,T)
            if isempty(xvec)
                xvec,sol = zeros(size(Qobj,1)),sol
            end
            # this is what will be concatenated into `results`:
        end
        #catch
        #    save("bad_matrices.jld","line",line)
        #end
        # is anything happening? report as each line is finished:
        # println("$(idx)/$(length(lines))")

        xvec,sol
    end

    return results
end

solve_temporal_instanton(inputData::InstantonInputData) = solve_temporal_instanton(
    inputData.Ridx,
    inputData.Y,
    inputData.G0,
    inputData.R0,
    inputData.D0,
    inputData.Sb,
    inputData.ref,
    inputData.lines,
    inputData.res,
    inputData.reac,
    inputData.k,
    inputData.line_lengths,
    inputData.line_conductors,
    inputData.Tamb,
    inputData.T0,
    inputData.int_length,
    inputData.corr
)

function process_instanton_results(results,n,nr,T;return_as_type=false)
    # Store results in more human-readable form:
    score = Float64[]
    α = Array(Vector{Float64},0)
    θ = Array(Array,0)
    x = Array(Array,0)
    diffs = Array(Array,0)
    xopt = Array(Array,0)

    for i in 1:size(results,1)
        xvec,sol = results[i]
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
    end
    if return_as_type
        return InstantonOutputData(score,x,θ,α,diffs,xopt)
    else
        return score,x,θ,α,diffs,xopt
    end
end

end
