using GraphLayout

"""
Return temperature trajectories corresponding to instances of
`InstantonInputData` and `InstantonOutputData`.
`eidx` is the event index, with 1 denoting the instanton.
"""
function temperatureTrajectory(
    i::InstantonInputData,
    o::InstantonOutputData,
    cp::ConductorParams,
    eidx,
    fixed_wind=[];
    )

    n = size(i.Y,1)
    nr = length(i.Ridx)
    numSteps = convert(Int64,length(i.G0)/n)

    # Thermal parameters:
    Tmid = (i.Tamb + cp.Tlim)/2
    temp_eq(t,T0,a,b) = (T0 + b/a).*exp(a*t) - b/a # solution to approx. heat balance IVP
    therm_a = cp.mCp\(-cp.ηc - 4*cp.ηr*(Tmid + 273)^3) # Fixed wrt power flow

    # Line parameters:
    idx = find(o.score.==sort(o.score)[eidx])[1]
    line = i.lines[i.line_lengths.!=0][idx]
    from,to = line
    r_ij = i.res[i.line_lengths.!=0][idx]
    x_ij = i.reac[i.line_lengths.!=0][idx]
    L_ij = i.line_lengths[i.line_lengths.!=0][idx]

    temp_trajectory = Vector()
    angle_dump = Vector()
    diffs_dump = Vector()

    if isempty(fixed_wind)
        fixed_wind = Array(Float64,0)
        for t in 1:numSteps
            append!(fixed_wind,o.x[idx][t])
        end
    end

    fixed_A = fixed_wind_A(numSteps,i.Y,i.ref,i.k)
    fixed_P = expand_renewable_vector(fixed_wind,i.Ridx,n,numSteps)
    fixed_b = fixed_wind_b(n,numSteps,i.G0,i.R0+fixed_P,i.D0)
    fixed_x = fixed_A\fixed_b
    angles,alpha = return_angles(fixed_x,n,numSteps)
    push!(angle_dump,angles)
    fixed_diffs = return_angle_diffs(angles,line)
    push!(diffs_dump,fixed_diffs)

    temp_values = [i.T0]
    power_flow = Float64[]

    T0 = i.T0
    for θij in fixed_diffs
        f_loss_pu = r_ij*(θij/x_ij)^2 # pu
        f_loss_si = f_loss_pu*i.Sb/(3*L_ij) # W/m
        push!(power_flow,(i.Sb/1e6)*θij/x_ij)
        therm_b = cp.mCp\( f_loss_si + cp.ηc*i.Tamb - cp.ηr*((Tmid + 273)^4 -
            (i.Tamb+273)^4) + 4*cp.ηr*Tmid*(Tmid + 273)^3 + cp.qs )
            temp_values = [temp_values;temp_eq(i.time_values,T0,therm_a,therm_b)[2:end]]
        T0 = temp_values[end]
    end
    return temp_values
end

"""
Use GraphLayout.jl to compute a spring layout
for the graph corresponding to a provided
adjacency matrix. Returns Void; useful only
for quickly obtaining network graphs in an
interactive environment.

Optional `C` controls spring layout tightness.
"""
function simple_graph_adj(adj,C=2.0)
    loc_x, loc_y = layout_spring_adj(adj,C=C)
    layout_tree
    labels = collect(1:size(adj,1))

    draw_layout_adj(adj,loc_x,loc_y,
    labels=labels,labelc="#000000",
    arrowlengthfrac=0.0,
    nodefillc="#FFFFFF",
    labelsize=2.0)
end
