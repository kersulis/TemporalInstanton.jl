type LineModel
    from
    to
    rij
    xij
    length
    D0
    mCp
    Ilim
    r
    Tlim
    ηc
    ηr
    qs
end

    """ Assign values to fields of LineModel type instance.
    Uses data from Mads's MPC paper.
    """
function add_thermal_parameters(line_model,conductor_name)
    if conductor_name == "waxwing"
        line_model.D0   = 15.5e-3
        line_model.mCp  = 383.
        line_model.Ilim = 439.
        line_model.r    = 110e-6
        line_model.Tlim = 65.
        line_model.ηc   = 0.955
        line_model.ηr   = 2.207e-9
        line_model.qs   = 14.4
    elseif conductor_name == "dove"
        line_model.D0   = 23.5e-3
        line_model.mCp  = 916.
        line_model.Ilim = 753.
        line_model.r    = 60e-6
        line_model.Tlim = 69.
        line_model.ηc   = 1.179
        line_model.ηr   = 3.346e-9
        line_model.qs   = 21.9
    end
    return line_model
end

""" Returns constant a [1/s]
mCp [J/m-C] is line heat capacity
ηc [W/m-C] is conductive heat loss rate coefficient
ηr [W/m-C^4] is radiative heat loss rate coefficient
Tamb [C] is ambient temperature (of air)
Tlim [C] is highest allowable line temperature
"""
function compute_a(mCp,ηc,ηr,Tamb,Tlim)
    Tmid = (Tamb + Tlim)/2
    return mCp\(-ηc - 4*ηr*(Tmid + 273)^3)
end

""" Return constant c [W/m]
mCp [J/m-C] is line heat capacity
r [pu] is line resistance
x [pu] is line reactance
Sb [W] is system base MVA
L [m] is line length
"""
function compute_c(mCp,r,x,Sb,L)
    return r*Sb/(3*mCp*L*(x^2))
end

""" Returns constant d [W/m]
mCp [J/m-C] is line heat capacity
ηc [W/m-C] is conductive heat loss rate coefficient
ηr [W/m-C^4] is radiative heat loss rate coefficient
Tamb [C] is ambient temperature (of air)
Tlim [C] is highest allowable line temperature
q_solar [W/m] is the solar heat gain rate
"""
function compute_d(mCp,ηc,ηr,Tamb,Tlim,q_solar)
    Tmid = (Tamb + Tlim)/2
    return mCp\(ηc*Tamb - ηr*((Tmid + 273)^4 - (Tamb + 273)^4) + 4*ηr*Tmid*(Tmid+273)^3 + q_solar)
end

""" Returns constant f
int_length [s] is length of each interval
a [1/s] is a constant
d [W/m] is a constant
n [-] is the number of time intervals
T0 [C] is the initial steady-state line temp
"""
function compute_f(int_length,a,d,n,T0)
    sum_coeff = sum([(e^(int_length*a))^i - (e^(int_length*a))^(i-1) for i in 1:n])
    return (e^(int_length*a))^n*T0 + (d/a)*sum_coeff
end


""" Return line's final temperature
int_length [s] is length of each interval
a [1/s] is a constant
c [W/m] is a constant
f [C] is a constant
n [-] is the number of time intervals
θij [rad] is the vector of angle differences (sorted by time interval)
"""
function compute_T(int_length,a,c,f,n,θij)
    angle_coeffs = [(e^(int_length*a))^i - (e^(int_length*a))^(i-1) for i in 1:n]
    return f + (c/a)*dot(angle_coeffs,flipud(θij).^2)
end
