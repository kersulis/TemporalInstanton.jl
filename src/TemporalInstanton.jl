module TemporalInstanton
export
    solve_instanton_qcqp, solve_temporal_instanton, LineParams,
    ConductorParams, load_rts96_data, InstantonInputData, InstantonOutputData, mat2tmpinst, process_instanton_results,

    # power flow:
    expand_renewable_vector,fixed_wind_A,fixed_wind_b,return_angles,
    return_angle_diffs,

    # collab:
    dict2type,ps2instantoninput,

    testcase

include("PowerSystem.jl")
include("ThermalModel.jl")
include("DataLoad.jl")
include("PowerFlow.jl")
include("QCQPMatrixBuilding.jl")
include("manipulations.jl")
include("PowerSystem.jl")
# include("plot.jl")
include("solvetmpinst.jl")
include("solvetmpinst-collab.jl")
end

@doc """
Top-level module for performing temporal instanton analysis.
See individual functions in included files for documentation.
""" -> TemporalInstanton
