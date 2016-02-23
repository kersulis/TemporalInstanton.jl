# single proc
include("../src/TemporalInstanton.jl")
using TemporalInstanton

inputData = load_rts96_data(return_as_type=true)

# Thermal model parameters:
inputData.Tamb = 35.0 # C
inputData.T0 = 60.0 #46.0 # initial line steady-state temp

inputData.time_values = 0.0:30.0:300.0 # five minutes in 30-sec steps
inputData.int_length = 300. # seconds = 5 min

Gp,Dp,Rp = inputData.G0,inputData.D0,inputData.R0
inputData.G0 = [0.7*Gp;0.7*Gp;0.7*Gp;0.7*Gp;0.7*Gp;0.7*Gp]
inputData.D0 = [0.9*Dp;0.9*Dp;0.9*Dp;0.9*Dp;0.9*Dp;0.9*Dp]
inputData.R0 = [Rp;1.1*Rp;1.2*Rp;1.3*Rp;1.4*Rp;1.5*Rp]

@time results = solve_temporal_instanton(inputData);
