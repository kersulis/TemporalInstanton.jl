module TemporalInstanton
export
    solve_instanton_qcqp, solve_temporal_instanton, LineParams,
    ConductorParams, load_rts96_data, createY,
    process_instanton_results, InstantonInputData, InstantonOutputData,
    mat2tmpinst,

    # temporary:
    tmp_inst_Qobj,tmp_inst_A1,tmp_inst_b,tmp_inst_Qconstr,
    return_conductor_params,return_thermal_constants,tmp_inst_A2,
    kernel_rotation,

    # power flow:
    expand_renewable_vector,fixed_wind_A,fixed_wind_b,return_angles,
    return_angle_diffs,

    # plot:
    temperatureTrajectory

include("DataLoad.jl")
include("PowerFlow.jl")
include("ThermalModel.jl")
include("QCQPMatrixBuilding.jl")
include("manipulations.jl")
include("SolveSecular.jl")
include("plot.jl")
include("solvetmpinst.jl")

end

@doc """
Top-level module for performing temporal instanton analysis. See individual functions in included files for documentation.
""" -> TemporalInstanton
