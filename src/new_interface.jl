type InputData
    "Fixed injections: nt vectors, each with length nb"
    fixedInj::Vector{Vector{Float64}}
    "Variable injections: can be any length"
    varInj::Vector{VarInj}
    participation::Vector{Float64}
    ysh::Vector{}
    "Vector of line objects"
    lines::Vector{LineParams}
    "Admittance matrix"
    Y::Array{Float64,2}
    Sb::Float64
end

"""
Variable injector object
"""
type VarInj
    bus
    weight
end

type ConductorParams
    "Current limit in A"
    Ilim
    "Heat capacity in "
    mCp
    "Diameter in m"
    D
    "Conductive heat rate coeff."
    eta_c
    "Radiative heat rate coeff."
    eta_r
    "Insolation"
    qs
end

type LineParams
    "Line from bus"
    from::Int64
    "Line to bus"
    to::Int64
    "Per unit resistance"
    r::Float64
    "Per unit reactance"
    x::Float64
    "Line length in m"
    len::Float64
    "Initial line temperature"
    T0::Float64
    "Line temperature limit"
    Tlim::Float64
    "Conductor parameters (uniquely determined by conductor)"
    cp::ConductorParams
end

type TemporalScanParams
    "in seconds"
    interval_length::Float64
    correlation::Array{Float64,2}
    "which lines to scan"
    scan_lines::Vector{Int64}
end

"""
Identify shifts in generation across a network that induce a desired temperature change in each line of `lines`. Solves a quadratically-constrained quadratic program.

Given transmission network data including admittance matrix, line and conductor information, a time horizon, and forecast values for known injections, identify most like (or lowest cost) shifts in a chosen set of generators that can bring a chosen line from its initial temperature to a desired final temperature by the end of the time horizon.
"""
function idenfity_deviation_patterns(d::InputData)
    nb,nt,nv = length(d.fixedInj[1]),length(d.fixedInj),length(d.varInj)
    decision_nodes = unique([Int64(v.bus) for v in d.varInj])

    obj_mat = objective_matrix(nb,nv,nt)
    obj_vec = zeros(size(obj_mat,1))
    obj_const = 0.0
    obj = (obj_mat,obj_vec,obj_const)

    lin_mat_fixed = fixed_linear_matrix(nv,nt,Y,ref,participation)
    kernel_fixed = kernel_basis(lin_mat_fixed[:,1:end-nt])
    lin_vec = linear_vector(nb,nt,fixedInj)

    quad_mat = quadratic_constraint_matrix(nb,nv,nt)
    quad_vec = zeros(size(quad_mat,1))

    # exclude lines with zero resistance or length
    skip_lines = zero_resistance_check(lines,silent)
    skip_lines = union(skip_lines,zero_length_check(lines,silent))
    # exclude lines detached from decision variables
    shift_factors = injection_shift_factors(Y,lines,ref,participation)[:,decision_nodes]
    skip_lines = union(skip_lines,shift_factor_check(shift_factors,silent))

    scan_lines = setdiff(scan_lines,skip_lines)

    # line loop
    results = @parallel (vcat) for idx in scan_lines
        tic()
        line = lines[idx]
        quad_const = quadratic_constraint_constant(line)
        quad = (quad_mat,quad_vec,quad_const)

        lin_mat_var = variable_linear_matrix(nb,varInj,nt,line,t_a,int_length)
        lin_mat = [lin_mat_fixed;lin_mat_var]
        lin = (lin_mat,lin_vec)

        xvec,sol,times = solve_qcqp(obj,quad,lin,kernel_fixed)

        # concatenate into `results`
        xvec,(sol,idx),toq(),times
    end
    process_results(results,nb,nv,nt)
end
