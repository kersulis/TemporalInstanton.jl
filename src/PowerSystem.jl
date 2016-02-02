"""
Translation of Jonathon Martin's `Bus` MATLAB structure.
"""
type Bus
    "bus number"
    num::Int64
    "bus type (1 = PQ, 2 = PV, 3 = slack(va), 4 = slack(qa))"
    ty::Int64
    "base voltage (kV)"
    vBase::Float64
    "voltage magnitude (pu)"
    vMag::Float64
    "voltage angle (rad)"
    vAng::Float64
    "minimum voltage magnitude (pu)"
    vMin::Float64
    "maximum voltage magnitude (pu)"
    vMax::Float64
    "shunt conductance (pu)"
    g::Float64
    "shunt susceptance (pu)"
    b::Float64
end

"""
Translation of Jonathon Martin's `reg` MATLAB structure
(nested within `Line`).

___

Contains regulation info for type 2, 3, and 4 lines.

*NOTES: Type 2 lines are assumed to regulate their 'to' bus voltage
magnitude while type 3 and 4 lines regulate the 'from' bus flows.
Regulation targets are midway between `regMin` and `regMax`. Settings in reg should be set to `NaN` for type 0 and 1 lines.*
"""
type Reg
    "minimum of regulation tolerance (in appropriate units)"
    regMin::Float64
    "maximum of regulation tolerance (in appropriate units)"
    regMax::Float64
    "minimum ratio (in tap or phase based on type)"
    ratioMin::Float64
    "maximum ratio (in tap or phase based on type)"
    ratioMax::Float64
    "ratio step size (in tap or radians based on type)"
    stepSize::Float64
end

"""
Translation of Jonathon Martin's `conductorModel` MATLAB structure
(nested within `Line`).
"""
type ConductorModel
    "unit resistance of conductor [ohms/m]"
    R_unit::Float64
    "current limit of network connection [A] (note that this is related to sLim field above. sLim = 0 or very large for voltage base will cause I_max to have a VERY large current rating (~50K) in order to approximate infinity."
    I_max::Float64
    "diameter of complete conductor [mm]"
    D::Float64
    "unit mass of aluminum portion of ACSR [kg/m]"
    Al_m::Float64
    "unit mass of steel portion of ACSR [kg/m]"
    St_m::Float64
    "number of conductors used in the network connection"
    bundle::Int64
end

"""
Translation of Jonathon Martin's `thermalModel` MATLAB structure
(nested within `Line`).
"""
type ThermalModel
    "ambient temperature [C]"
    T_amb::Float64
    "ampacity-limited steady-state temperature [C] (related to sLim and I_max)"
    T_max::Float64
    "conductive heat factor (varies w/ temp) [W/m-C]"
    eta_c::Vector{Float64}
    "radiative heat factor [W/m-C^4]"
    eta_r::Float64
    "linear solar heat gain rate [W/m]"
    eta_s::Float64
    "line heat capacity per unit length (J/m-C)"
    mCp::Float64
    "(unitless) coefficient for temperature dynamics"
    tau::Float64
    "(C-m/W) coefficient for temperature dynamics"
    rho::Float64
    "(unitless) coefficient for temperature dynamics"
    gamma::Float64
    "per-phase active power loss linearization point [W/m]"
    fLossLim::Float64
end

"""
Translation of Jonathon Martin's `Line` MATLAB structure.

___

Unified branch model:
```
        k  tmk:1    p     Zkm      q               m
        |----8 8------O-----ZZZ------O-------------|
  Tap   +             |              |             +   Fixed
  bus   Vk     ykm^sh Z              Z ykm^sh      Vm   bus
        _             |              |             -
       ---------------------------------------------
```
tmk = tap x exp(j x phase) is pu (complex) tap-ratio.

* when tmk = 1: transmission line
* when tmk = real: TCT, tap-ratio (~1)
* when tmk = complex: PST, angle (rad)
"""
type Line
    "from bus number"
    from::Int64
    "to bus number"
    to::Int64
    "line type (0=T-line, 1=fixed Xfrmr, 2=vReg Xfrmr, 3=qReg Xfrmr, 4=phaseReg Xfrmr)"
    ty::Int64
    "resistance (pu)"
    R::Float64
    "reactance (pu)"
    X::Float64
    "shunt susceptance (pu)"
    B::Float64
    "3-phase power rating of line (MVA)"
    sLim::Float64
    "length of line (miles) (Level 3 can estimate lengths if data not available)"
    length::Float64
    "active power flow on line (pu)"
    f::Float64
    "active power losses on line (pu)"
    fLoss::Float64
    "reactive power flow on line (pu)"
    q::Float64
    "reactive power losses on line (pu)"
    qLoss::Float64
    "degrees above thermal line limit (C)"
    delT::Float64
    "tap ratio "
    tap::Complex{Float64}
    "phase ratio (rad)"
    phase::Float64
    "contains regulation info for type 2, 3, and 4 lines"
    reg::Reg
    "**\$**"
    conductorModel::ConductorModel
    "**\$**"
    thermalModel::ThermalModel
    "**\$** width of PWL segments approximating abs(theta_{ij})"
    delTheta::Float64
end

"""
Translation of Jonathon Martin's `Gen` MATLAB structure.
"""
type Gen
    "bus number of generator"
    bus::Int64
    "generator setting (0=offline, 1=online-Vreg, 2=online-noVreg)"
    status::Int64
    "active power injection (pu)"
    Pinj::Float64
    "reactive power injection (pu)"
    Qinj::Float64
    "active power ramping rate (pu/hr)"
    Pramp::Float64
    "voltage setpoint of generator (pu)"
    vSet::Float64
    "min active power injection (pu)"
    Pmin::Float64
    "max active power injection (pu)"
    Pmax::Float64
    "min reactive power injection (pu) (should be <= 0)"
    Qmin::Float64
    "max reactive power injection (pu) (should be >= 0)"
    Qmax::Float64
    "active power ramp rate limit (pu/hr) (for increasing output)"
    pUpLim::Float64
    "active power ramp rate limit (pu/hr) (for decreasing output)"
    pDownLim::Float64
    "reactive power ramp rate limit (pu/hr) (for increasing output)"
    qUpLim::Float64
    "reactive power ramp rate limit (pu/hr) (for decreasing output)"
    qDownLim::Float64
end

"""
Translation of Jonathon Martin's `Load` MATLAB structure.
"""
type Load
    "bus number of load"
    bus::Int64
    "nominal active power demand forecast (pu) 1x(number of time intervals) vector"
    Pload::Float64
    "nominal reactive power demand forecast (pu) 1x(number of time intervals) vector"
    Qload::Float64
    "present active power load (pu)"
    Pdemand::Float64
    "present reactive power load (pu)"
    Qdemand::Float64
    "active power shed from load through demand response (pu)"
    Pshed::Float64
    "reactive power shed from load through demand response (pu)"
    Qshed::Float64
    "3x1 vector specifying [S,I,Z] portion of active load (sum to 1)"
    Psiz::Vector{Float64}
    "3x1 vector specifying [S,I,Z] portion of reactive load (sum to 1)"
    Qsiz::Vector{Float64}
end

"""
Translation of Jonathon Martin's `Wind` MATLAB structure.
"""
type Wind
    "bus number of wind generator"
    bus::Int64
    "forecasted active power production (pu) 1x(number of time intervals) vector"
    Pforecast::Float64
    "actual active power production (pu)"
    Pinj::Float64
    "active power shed from wind forecast through curtailment (pu)"
    Pshed::Float64
end

"""
Translation of Jonathon Martin's `Storage` MATLAB structure.
"""
type Storage
    "bus number of storage"
    bus::Int64
    "active power being consumed by device [measured at bus] (pu)"
    Pdemand::Float64
    "power rating (pu)"
    powerRate::Float64
    "energy rating (pu*hr)"
    energyRate::Float64
    "energy currently stored in device (pu*hr)"
    energyLvl::Float64
    "charging efficiency [between 0 and 1]"
    chargeEff::Float64
    "discharging efficiency [between 0 and 1]"
    dischargeEff::Float64
end

"""
Translation of Jonathon Martin's `Params` MATLAB structure.
"""
type Params
    "3-phase power system base (MVA)"
    sBase::Float64
    "No. of buses in system"
    numOfBus::Int64
    "No. of lines in system"
    numOfLine::Int64
    "No. of generators in system"
    numOfGen::Int64
    "No. of loads in system"
    numOfLoad::Int64
    "No. of wind generators in system"
    numOfWind::Int64
    "No. of storage devices in system"
    numOfStorage::Int64
    "No. of time intervals (for load and wind forecasts)"
    numOfTimeInt::Int64
    "integer indicating how many segments to use in piecewise-linear model of system losses. (0 indicates that PWL model is not used)"
    hasLosses::Int64
    "maximum phase angle difference allowed across a single line in the network [rad]"
    phaseMax::Float64
    "amount of text to display on command window (0  (default) = none, 1 = overall steps, 2 = full detail)"
    displayFlag::Int64
    "slack participation factors (nbus x 1)"
    participation::Vector{Float64}
    "**\$** node arc incidence matrix (nbus x nline, 1 = from -1 = to)"
    NAI
    "**\$** node generator incidence matrix (nbus x ngen, 1 = incidence)"
    NGI
    "**\$** node load incidence matrix (nbus x nload, 1 = incidence)"
    NLI
    "**\$** node wind incidence matrix (nbus x nwind, 1 = incidence)"
    NWI
    "**\$** node storage incidence matrix (nbus x nstorage, 1 = incidence)"
    NSI
    "**\$** PWL approximation of abs(theta_{ij}) for each line"
    phPW_b
end

"""
Translation of Jonathon Martin's `Control` MATLAB structure.
"""
type Control
    "discrete sampling time (for MPC) [sec] (Must set to 60 for now)"
    T_s::Float64
    "receeding prediction horizon for MPC"
    MPChorizon::Int64
    "number of MPC iterations needed without line trip for system to be considered safe"
    MPCsafeTime::Int64
    "maximum number of MPC iterations allowed"
    MPCmaxTime::Int64
    "nline x 1 logical specifying which lines are out of service"
    lineOut::Vector{Bool}
    "scalar specifying tolerance around 0 for which storage devices are allowed to simulcharge"
    storageTol::Float64
    "acceptable line flow overload allowed at safe MPC termination [%]"
    overloadTol::Float64
end

"""
Translation of Jonathon Martin's `params` MATLAB structure
(nested within `Optimize`).
"""
type GurobiParams
    "decide whether GUROBI steps should display on console window (1=yes,0=no)"
    outputflag::Int64
    "max time GUROBI is allowed to operate (sec)"
    timelimit::Float64
    "determine aggregator setting (set to 1 to aggregate)"
    aggregate::Int64
end

"""
Translation of Jonathon Martin's `costFunction` MATLAB structure
(nested within `Optimize`).
"""
type CostFunction
    "penalty for temperature overload [\$/C^2]"
    delThat::Float64
    "(terminal) cost of deviating from slow-timescale (interpolated) energy storage levels [\$/pu^2]"
    Elvl::Float64
    "(terminal) penalizing generator levels from interpolated references [\$/pu^2]"
    fG::Float64
    "generator ramping cost factor that adds to quadratic gen cost term [\$/pu^2]"
    dfG::Float64
    "penalty for spilling wind from wind farms [\$/pu^2]"
    windSpill::Float64
    "penalty for shedding load through demand response [\$/pu^2]"
    loadShed::Float64
    "storage utility cost [\$/pu^2] (i.e. keeping storage level constant is NOT penalized)"
    fQ::Float64
    "phase-shift from PST wrt Level 1/2 reference [\$/rad^2]"
    PSTshift::Float64
    "penalty for voltage limit violations [\$/pu^2]"
    Vm::Float64
end

"""
Translation of Jonathon Martin's `Optimize` MATLAB structure.
"""
type Optimize
    "optimization solver 2 == GUROBI (1 == CPLEX, currently unavailable)"
    QPsolver::Int64
    "GUROBI parameters"
    params::GurobiParams
    costFunction::CostFunction
end

"""
Translation of Jonathon Martin's `Time` MATLAB structure.
"""
type Time
    "scalar current day of year"
    day::Int64
    "scalar current hour of day"
    hour::Int64
    "scalar current minute of hour"
    minute::Int64
end

"""
Translation of Jonathon Martin's `System` MATLAB structure.

___

Description of a power system, which contains
information about the network connections and devices.
Fields which are marked with a dollar sign '**\$**' are initialized
or added to System inside a Level 3 function. Unmarked fields
should be initialized prior to a Level 3 routine being called.

"""
type PowerSystem
    bus::Vector{Bus}
    line::Vector{Line}
    gen::Vector{Gen}
    load::Vector{Load}
    wind::Vector{Wind}
    storage::Vector{Storage}
    params::Params
    control::Control
    optimize::Optimize
    time::Time
end

"""
Convert one of Jon's MATLAB structs (dict in Julia) into an
instance of the PowerSystem type.
"""
function dict2type(sys::Dict{ASCIIString,Any})
    b = sys["Bus"]
    bus = [Bus(b["num"][i],
        b["ty"][i],
        b["vBase"][i],
        b["vMag"][i],
        b["vAng"][i],
        b["vMin"][i],
        b["vMax"][i],
        b["g"][i],
        b["b"][i])
        for i in 1:length(b["num"])]

    r = sys["Line"]["reg"]
    reg = [Reg(
        r[i]["regMin"],
        r[i]["regMax"],
        r[i]["ratioMin"],
        r[i]["ratioMax"],
        r[i]["stepSize"])
        for i in 1:length(r)]

    c = sys["Line"]["conductorModel"]
    conductormodel = Vector{ConductorModel}()
    for i in 1:length(c)
        if sys["Line"]["ty"][i] == 1
            push!(conductormodel,ConductorModel(
            NaN,NaN,NaN,NaN,NaN,typemin(Int64)))
        else
            push!(conductormodel,ConductorModel(
            c[i]["R_unit"],
            c[i]["I_max"],
            c[i]["D"],
            c[i]["Al_m"],
            c[i]["St_m"],
            c[i]["bundle"]))
        end
    end

    t = sys["Line"]["thermalModel"]
    thermalmodel = Vector{ThermalModel}()
    for i in 1:length(t)
        if sys["Line"]["ty"][i] == 1
            push!(thermalmodel,ThermalModel(
            NaN,NaN,[NaN],NaN,NaN,NaN,NaN,NaN,NaN,NaN))
        else
            push!(thermalmodel,ThermalModel(
            t[i]["T_amb"],
            t[i]["T_max"],
            t[i]["eta_c"][:],
            t[i]["eta_r"],
            t[i]["eta_s"],
            t[i]["mCp"],
            t[i]["tau"],
            t[i]["rho"],
            t[i]["gamma"],
            t[i]["fLossLim"]))
        end
    end

    l = sys["Line"]
    line = [Line(
        l["from"][i],
        l["to"][i],
        l["ty"][i],
        l["R"][i],
        l["X"][i],
        l["B"][i],
        l["sLim"][i],
        l["length"][i],
        l["f"][i],
        l["fLoss"][i],
        l["q"][i],
        l["qLoss"][i],
        l["delT"][i],
        l["tap"][i],
        l["phase"][i],
        reg[i],
        conductormodel[i],
        thermalmodel[i],
        NaN)
        for i in 1:length(l["from"])]

    g = sys["Gen"]
        gen = [Gen(
        g["bus"][i],
        g["status"][i],
        g["Pinj"][i],
        g["Qinj"][i],
        g["Pramp"][i],
        g["vSet"][i],
        g["Pmin"][i],
        g["Pmax"][i],
        g["Qmin"][i],
        g["Qmax"][i],
        g["pUpLim"][i],
        g["pDownLim"][i],
        g["qUpLim"][i],
        g["qDownLim"][i])
        for i in 1:length(g["bus"])]

    l = sys["Load"]
    load = [Load(
        l["bus"][i],
        l["Pload"][i],
        l["Qload"][i],
        l["Pdemand"][i],
        l["Qdemand"][i],
        l["Pshed"][i],
        l["Qshed"][i],
        l["Psiz"][i][:],
        l["Qsiz"][i][:])
        for i in 1:length(l["bus"])]

    w = sys["Wind"]
    wind = [Wind(
        w["bus"][i],
        w["Pforecast"][i],
        w["Pinj"][i],
        w["Pshed"][i])
        for i in 1:length(w["bus"])]

    s = sys["Storage"]
    storage = [Storage(
        s["bus"][i],
        s["Pdemand"][i],
        s["powerRate"][i],
        s["energyRate"][i],
        s["energyLvl"][i],
        s["chargeEff"][i],
        s["dischargeEff"][i])
        for i in 1:length(s["bus"])]

    # why isn't pfThreshold included in Params?
    p = sys["Params"]
    params = Params(
    p["sBase"],
    p["numOfBus"],
    p["numOfLine"],
    p["numOfGen"],
    p["numOfLoad"],
    p["numOfWind"],
    p["numOfStorage"],
    typemin(Int64),
    p["hasLosses"],
    p["phaseMax"],
    p["displayFlag"],
    p["participation"][:],
    NaN,
    NaN,
    NaN,
    NaN,
    NaN,
    NaN)

    c = sys["Control"]
    control = Control(
    c["T_s"],
    c["MPChorizon"],
    c["MPCsafeTime"],
    c["MPCmaxTime"],
    c["lineOut"][:],
    c["storageTol"],
    c["overloadTol"])

    g = sys["Optimize"]["params"]
    gurobiparams = GurobiParams(g["outputflag"],NaN,typemin(Int64))

    c = sys["Optimize"]["costFunction"]
    costfunction = CostFunction(
    c["delThat"],
    c["Elvl"],
    c["fG"],
    c["dfG"],
    c["windSpill"],
    c["loadShed"],
    c["fQ"],
    c["PSTshift"],
    c["Vm"])

    optimize = Optimize(
    sys["Optimize"]["QPsolver"],
    gurobiparams,
    costfunction)

    t = sys["Time"]
    pstime = Time(
    typemin(Int64),
    t["hour"],
    t["minute"]
    )

    PowerSystem(
    bus,
    line,
    gen,
    load,
    wind,
    storage,
    params,
    control,
    optimize,
    pstime)
end
