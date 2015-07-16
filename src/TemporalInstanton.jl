module TemporalInstanton

using HDF5, JLD, ProgressMeter, IProfile

export
    solve_instanton_qcqp, solve_temporal_instanton, LineModel,
    # temporary:
    tmp_inst_Qobj,tmp_inst_pad_Q,tmp_inst_A,tmp_inst_b,tmp_inst_pad_b,
    tmp_inst_Qtheta,add_thermal_parameters,compute_a,compute_c,
    compute_d,compute_f,tmp_inst_A_scale_new

include("PowerFlow.jl")
include("ThermalModel.jl")
include("QCQPMatrixBuilding.jl")
include("manipulations.jl")
include("SolveSecular.jl")

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
    T = round(Int64,length(find(P0))/nr)

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

    # addprocs(3)
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

        # array of vectors with Float64 values:
        deviations = Array(Vector{Float64},0)
        angles = Array(Vector{Float64},0)
        alpha = Float64[]

        # Create A2 based on chosen line:
        A2 = tmp_inst_A_scale_new(n,Ridx,T,line,therm_a,int_length)
        # Stack A1 and A2:
        A = [A1; A2]

        # Computationally expensive part: solving QCQP
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
    #rmprocs([2,3,4])
    return score,x,θ,α,diffs,xopt
end

end
