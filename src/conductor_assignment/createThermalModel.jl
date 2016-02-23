"""
This function processes the thermal data of the network of interest in
order to determine the conductor temperature dynamic coefficients for use
in the MPC model.

Inputs:
System: structure containing network information

Output:
model: structure containing thermalModel structure for each line
Fields of thermalModel are
  .T_amb      ambient temperature [C]
  .T_max      ampacity-limited steady-state temperature [C]
  .eta_c      convective heat factor (varies w/ temp) 
  .eta_r      radiative heat factor
  .eta_s      linear solar heat gain rate
  .mCp        line heat capacity per unit length
  .tau        
  .rho
  .gamma

The equation of the temperature dynamics of each transmission line is:
m*C_p dT/dt = Q_l(t) + Q_s - Q_c(T) - Q_r(T^4)
where:
Q_l - ohmic I(t)^2*R losses of conductor (depends on power), [W/m]
Q_s - rate of solar heat gain, [W/m]
Q_c - forced convection (wind) heat loss rate, [W/m]
Q_r - rate of radiated heat loss, [W/m]
"""
function createThermalModel(ps)
	## Set up data constants and initialize coefficients
	# Define whether line resistance is temperature dependent
	isTempDependentResistance = false
	# Tref is the reference temp when resistance is temp dependent
	Tref = 200.

	# Steady-state temperatures of interest to study
	maxTemp = 200.             # set the maximum conductor temp [C] --- if max(T) is too low, some lines have incorrect Tlim
	T       = 1:maxTemp       # conductor temperatures to study [C] 
	T_amb   = 35.              # ambient temp [C]
	nLines  = ps.params.numOfLine    # determine how many network connections there are

	## Initialize storage variables
	# T_max = zeros(nLines,1)
	# mCp   = zeros(nLines,1)
	# eta_r = zeros(nLines,1)
	# eta_c = zeros(nLines,nTemps)
	# eta_s = zeros(nLines,1)

	thermalmodels = Vector()
	for idx = 1:nLines
	    if !isnan(ps.line[idx].conductorModel.I_max)
		# This network connection has a conductor model so it is not a
		# transformer and has realistic flow limit information.
		
		# make local copies of relevant line information
		Al_m = ps.line[idx].conductorModel.Al_m        # aluminum mass of conductor [kg/m]
		St_m = ps.line[idx].conductorModel.St_m        # steel mass of conductor [kg/m]
		D = ps.line[idx].conductorModel.D              # diameter of entire conductor [mm]
		R_unit = ps.line[idx].conductorModel.R_unit    # unit resistance of conductor [ohms/m]
		bundle = ps.line[idx].conductorModel.bundle    # # of conductors in line bundle
		I_max = ps.line[idx].conductorModel.I_max      # max current capacity of line [A]
		T_s = ps.control.T_s                           # sampling time [sec]
		sBase = ps.params.sBase                        # 3 phase power base of system [MVA]
		
		He = 61.        # average conductor height from sealevel, [m] \\assumed avg PJM elevation (200 ft)
		L  = 1.         # unit length [m] / [old]: Length of line (m), NOTE: NEED TO CONFIRM WHAT MADS WAS DOING HERE.  DO I NEED LINE LENGTHS?? ********
		
		# Conservative ambient conditions:
		V_w   = 0.61           # wind speed [m/s] [2 ft/s = 0.61 m/s]
		ang_w = pi/2           # wind/line angle (perpendicular = pi/2; parallel = 0)
		emm   = 0.7            # emmissitivity [0.23, 0.91]
		
		# Solar-specific values:
		# NOTE: dateNow is inactivated due to days365 function not being
		# available below.
		dateNow = "10-June-2012"   # date of interest (e.g. 'today' = date);
		hourNow = 12               # hour of the day (1:24)
		Lat     = 40*pi/180        # latitude of line (on Earth)
		Z_l     = pi/2             # line azimuth (e.g. east-to-west = pi/2; South-to-north = 0);
		alpha   = 0.9              # solar absorbtivity [0.23, 0.91]

		## m*C_p - heat capacity :
		# Assumes ACSR cables are used in ThermalModel_data_(network)
		# Cp: {Alu, Steel, Copper, Alu-Clad Steel} = {955, 476, 390, 530} [J/kg-C] CONFIRM THESE VALUES********
		Al_Cp = 955.
		St_Cp = 476.
		# mass = aluminum (strand) + steel (core)
		# mCp = mCp{Alu} + mCp{Steel}
		# L*J/m-C = 1*J/m-C = (J/m-C)
		mCp = L*(Al_m*Al_Cp + St_m*St_Cp)
		
		## Q_c - forced convection heat loss:
		# Q_c = eta_c*(T - T_amb), eta depends on wind, temp, diameter, height
		
		# eta_c (low wind) = [1.01 + 0.0372*(D*p_f*V_w/ mu_f)^.52]*k_f*K_angle
		# D - diameter [mm]
		# p_f - air density [kg/m^3], range [0.574, 1.293], depends on temp & height [see Table 2, std-738]
		# mu_f - dynamic viscosity [Pa-s], range [1.72e-5, 2.17e-5], depends on temp [see Table 2, std-738]
		# k_f - thermal cond of air [W/m-C], range [0.0242,0.0317], depends on temp [see Table 2]
		# V_w - wind speed [m/s], range [0.2, 2.2]
		# K_angle - wind direction factor = 1.194-cos(a)+0.194*cos(2*a)+0.368*sin(2*a), where
		# a - angle between wind direction and conductor axis
		
		# calculated values:
		T_film  = 0.5*(T+T_amb) # film temperature, (T+Ta)/2, [C]
		p_f     = (1.293 - 1.525e-4*He + 6.379e-9*He^2)./(1+0.00367*T_film)
		mu_f    = 1.458e-6*(T_film +273).^1.5./(T_film+383.4)
		k_f     = 2.424e-2 + 7.477e-5*T_film - 4.407e-9*T_film.^2
		K_angle = 1.194-cos(ang_w)+0.194*cos(2*ang_w)+0.368*sin(2*ang_w)
		
		eta_c = L*max(0.0119*(D.*p_f.*V_w./ mu_f).^.60, (1.01 + 0.0372*(D.*p_f.*V_w./ mu_f).^.52)).*k_f*K_angle # L*W/m-C = 1*W/m-C = W/m-C
		q_c = eta_c.*(T - T_amb) # [W/m]
		
		## Q_r - radiated heat loss:
		# depends mainly on solar output
		# Q_r = eta_r*((T+273)^4 - (T_amb+273)^4;
		
		# eta_r = 0.0178*D*e/100^4
		# e = emmissivity, range [0.23, 0.91], unity
		
		# Calculate values:
		eta_r = L*0.0178*D*emm/100^4 # L*W/m-C^4 = 1*W/m-C^4 = W/C^4
		
		q_r = eta_r*((T+273).^4 - (T_amb+273)^4) # [W/m]
		
		## Q_s - rate of solar gain
		# q_s = eta_s
		
		# eta_s = alpha*Q_se*sin(b)*A'
		# alpha - solar absorbtivity, range [0.23, 0.91]
		# Q_se - solar+sky radiated heat flux rate= = K_solar*Q_s , (W/m^2)
		# K_solar - solar heat multiplying factor = 1 + 1.148*He - 1.108*He^2;
		# Q_s - heat flux received at sea lvl
		#    = -42.2391+63.8044*Hc-1.9220*Hc^2+3.46921e-2*Hc^3-3.61118e-4*Hc^4+1.94318e-6*Hc^5-4.07608e-9*Hc^6
		# b - angle of incidence of the sun rays = acos(cos(Hc)*cos(Zc?Zl))
		# Hc - altitude of sun (degrees)=asin(cos(Lat)cos(c)cos(w)+sin(Lat)sin(c))
		# Lat - degrees of latitude (degrees)
		# c - solar declination angle = 23.4583*sin(0.9863*(284+N))
		# N - day of the year (i.e. Jan 21 = 21 and Feb 12 = 43). (integer)
		# w - number of hours from noon x 15 degrees (i.e. 10AM = -30) (degrees)
		# Zc - azimuth of sun (degrees) = C_c + atan(X_c)
		# C_c - solar azimuth constant: {0, 180, 360} depends on w and X_c
		# X_c - solar azimuth variable = sin(w)/(sin(Lat)*cos(w)?cos(Lat)*tan(c))
		# Zl - azimuth of line (degrees)
		# A_l - projected area of conductor per unit length, (m^2/m) = D/1000
		
		# Determine daily calculations:
		#N = days365('1-Jan-2012',dateNow);  % unique idx for each day
		N = 160 #NOTE CHANGED due to days365 function not being available.  Based on June 10, 2012 as 'dateNow'********
		c = 23.4583*sin(0.9863*(284+N)*pi/180)*pi/180 # unique for each day;
		# determine hourly values:
		w = (hourNow-12)*15*pi/180 # may have multiple values (if hourNow=1:24)
		
		# calculate values:
		Hc = asin(cos(Lat)*cos(c)*cos(w) + sin(Lat)*sin(c)) # may have mult vals
		X_c = sin(w)./(sin(Lat)*cos(w)-cos(Lat)*tan(c)) # may have mult vals
		# compute
		C_c = zeros(size(X_c))
		for i=1:length(X_c)
		    if X_c[i] >= 0
			    w[i] < 0 ? C_c[i] = 0 : C_c[i] = pi
		    else
			    w[i] < 0 ? C_c[i] = pi : C_c[i] = 2*pi
		    end
		end
		Z_c = C_c + atan(X_c) # may take vector form
		b = acos(cos(Hc).*cos(Z_c-Z_l)) # may take vector form
		
		K_solar = 1 + 1.148e-4*He - 1.108e-8*He^2
		Hc_deg = Hc*180/pi
		Q_s = -42.2391 + 63.8044*Hc_deg - 1.9220*Hc_deg.^2 + 3.46921e-2*Hc_deg.^3 - 3.61118e-4*Hc_deg.^4 + 1.94318e-6*Hc_deg.^5 - 4.07608e-9*Hc_deg.^6
		Q_se = K_solar*Q_s
		
		A_l = D/1000 # m^2/m = m
		
		eta_s = L*alpha.*sin(b).*A_l.*Q_se # L*W/m = 1*W/m
		eta_s = eta_s[1]
		
		q_s = eta_s # [W/m]
		
		## Find ohmic losses and determine temp limit
		# Ohmic losses (per-phase):
		# Define resistance in terms of temperature:
		if isTempDependentResistance
		    rT = R_unit[idx]/bundle[idx]*(1+0.0039*(T-Tref))
		    q_L = rT.*I_max^2 # W/m
		else
		    q_L = R_unit/bundle*I_max^2 # [W/m]
		end
		# Find steady-state T_max of each conductor:
		# (per unit distance: W/m, Ohm/m)
		E_total = (q_L + q_s - q_c - q_r)/(sBase*1e6) # 100 MW = 1 p.u.
		E_zero  = minabs((q_L + q_s - q_c - q_r))/(sBase*1e6)
		T_max   = minimum([T[E_total.==E_zero];T[E_total.==-E_zero]])
		
		## Determine dynamic coefficients
		# line parameters
		gamma_c_bar = eta_c[round(Int64,T_max)] + 4*eta_r*(T_max + 273)^3     # linearization term
		gamma_a_bar = eta_c[round(Int64,T_max)] + 4*eta_r*(T_amb + 273).^3    # linearization term

		# System dynamics:
		# dT[k+1] = tau*dT[k] + rho*df_loss[k] + rho*dq_s[k] + gamma*dT_amb[k],
		# where:
		tau   = 1 - T_s*gamma_c_bar/mCp             # unitless
		gamma = T_s*gamma_a_bar/mCp                 # unitless
		rho   = T_s/mCp                             # s/(J/C) = C/(W/m)
		
		
		## store the results for output
		# model.thermalModel(idx) = struct('T_amb',T_amb,'T_max',T_max,'eta_c',eta_c(T_max),'eta_r',eta_r,'eta_s',eta_s,'mCp',mCp,'tau',tau,'rho',rho,'gamma',gamma,'fLossLim',q_L)
		push!(thermalmodels,ThermalModel(T_amb,T_max,eta_c[round(Int64,T_max)],eta_r,eta_s,mCp,tau,rho,gamma,q_L))
	    else
		# This connection does not have a conductor model so it is either a
		# transformer or has unrealistic flow limits.
		
		# model.thermalModel(idx) = struct('T_amb',[],'T_max',[],'eta_c',[],'eta_r',[],'eta_s',[],'mCp',[],'tau',1,'rho',0,'gamma',0,'fLossLim',0)
		push!(thermalmodels,NaN)
	    end
	end
	return thermalmodels
end
