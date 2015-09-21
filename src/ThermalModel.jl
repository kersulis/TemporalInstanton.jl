"""
Line parameters that vary with the line (not inherent to conductor material).
"""
type LineParams
    "Line 'from' node"
    from::Int64
    "Line 'to' node"
    to::Int64
    "pu line resistance"
    rij::Float64
    "pu line reactance"
    xij::Float64
    "Line length in meters"
    length::Float64
end

"""
Line parameters that depend on conductor material.
"""
type ConductorParams
    "[m] conductor diameter"
    D0::Float64
    "[J/m-C] line heat capacity"
    mCp::Float64
    "[A] max. allowable current"
    Ilim::Float64
    "[pu] line resistance"
    r::Float64
    "[C] highest allowable line temperature"
    Tlim::Float64
    "[W/m-C] conductive heat loss rate coefficient"
    ηc::Float64
    "[W/m-C^4] radiative heat loss rate coefficient"
    ηr::Float64
    "[W/m] solar heat gain rate"
    qs::Float64
end

"""
Returns instance of `ConductorParams` filled with
conductor parameters. Data from Mads's MPC paper.

Currently accepts "waxwing" and "dove" as arguments.
"""
function return_conductor_params(conductor::String)
    if conductor == "waxwing"
        return ConductorParams(15.5e-3,383.0,439.0,110e-6,65.0,0.955,2.207e-9,14.4)
    elseif conductor == "dove"
        return ConductorParams(23.5e-3,916.0,753.0,60e-6,69.0,1.179,3.346e-9,21.9)
    end
end

"""
Return (a,c,d,f) thermal constants used in line temperature IVP. Arguments:

* `lp` instance of LineParams
* `cp` instance of ConductorParams
* `Tamb` [C] ambient temperature
"""
function return_thermal_constants(lp,cp,Tamb,Sb,int_length,n,T0)
    a = compute_a(cp.mCp,cp.ηc,cp.ηr,Tamb,cp.Tlim)
    c = compute_c(cp.mCp,lp.rij,lp.xij,Sb,lp.length)
    d = compute_d(cp.mCp,cp.ηc,cp.ηr,Tamb,cp.Tlim,cp.qs)
    f = compute_f(int_length,a,d,n,T0)
    return a,c,d,f
end

"""
Returns constant `a` [1/s]. Arguments:

* `mCp` [J/m-C] is line heat capacity
* `ηc` [W/m-C] is conductive heat loss rate coefficient
* `ηr` [W/m-C^4] is radiative heat loss rate coefficient
* `Tamb` [C] is ambient temperature (of air)
* `Tlim` [C] is highest allowable line temperature
"""
function compute_a(mCp,ηc,ηr,Tamb,Tlim)
    Tmid = (Tamb + Tlim)/2
    return mCp\(-ηc - 4*ηr*(Tmid + 273)^3)
end

"""
Return constant `c` [W/m]. Arguments:

* `mCp` [J/m-C] is line heat capacity
* `r` [pu] is line resistance
* `x` [pu] is line reactance
* `Sb` [W] is system base MVA
* `L` [m] is line length
"""
function compute_c(mCp,r,x,Sb,L)
    return r*Sb/(3*mCp*L*(x^2))
end

"""
Returns constant `d` [W/m]. Arguments:

* `mCp` [J/m-C] is line heat capacity
* `ηc` [W/m-C] is conductive heat loss rate coefficient
* `ηr` [W/m-C^4] is radiative heat loss rate coefficient
* `Tamb` [C] is ambient temperature (of air)
* `Tlim` [C] is highest allowable line temperature
* `q_solar` [W/m] is the solar heat gain rate
"""
function compute_d(mCp,ηc,ηr,Tamb,Tlim,q_solar)
    Tmid = (Tamb + Tlim)/2
    return mCp\(ηc*Tamb - ηr*((Tmid + 273)^4 - (Tamb + 273)^4) + 4*ηr*Tmid*(Tmid+273)^3 + q_solar)
end

"""
Returns constant `f` [C]. Arguments:

* `il` [s] is time interval length
* `a` [1/s] is a constant
* `d` [W/m] is a constant
* `n` [unitless] is the number of time intervals
* `T0` [C] is the initial steady-state line temp
"""
function compute_f(il,a,d,n,T0)
    sum_coeff = sum([(e^(il*a))^i - (e^(il*a))^(i-1) for i in 1:n])
    return (e^(il*a))^n*T0 + (d/a)*sum_coeff
end


"""
Return line's final temperature. Arguments:

* `int_length` [s] is length of each interval
* `a` [1/s] is a constant
* `c` [W/m] is a constant
* `f` [C] is a constant
* `n` [unitless] is the number of time intervals
* `θij` [rad] is the vector of angle differences (sorted by time interval)
"""
function compute_T(int_length,a,c,f,n,θij)
    angle_coeffs = [(e^(int_length*a))^i - (e^(int_length*a))^(i-1) for i in 1:n]
    return f + (c/a)*dot(angle_coeffs,flipud(θij).^2)
end
