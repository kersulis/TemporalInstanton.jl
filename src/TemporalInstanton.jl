module TemporalInstanton
export
    solve_instanton_qcqp, solve_temporal_instanton, LineParams,
    ConductorParams, load_rts96_data, InstantonInputData, InstantonOutputData, mat2tmpinst,

    # power flow:
    expand_renewable_vector,fixed_wind_A,fixed_wind_b,return_angles,
    return_angle_diffs,

    testcase

include("DataLoad.jl")
include("PowerFlow.jl")
include("ThermalModel.jl")
include("QCQPMatrixBuilding.jl")
include("manipulations.jl")
include("PowerSystem.jl")
# include("plot.jl")
include("solvetmpinst.jl")
end

@doc """
Top-level module for performing temporal instanton analysis.
See individual functions in included files for documentation.
""" -> TemporalInstanton
