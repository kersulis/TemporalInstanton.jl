"""
Convenience method for converting an instance of PowerSystem to InstantonInputData format.
"""
function ps2instantoninput(ps::PowerSystem)
    # Admittance matrix
    lineOut = ps.control.lineOut
    from = [Int64(l.from) for l in ps.line]
    to = [Int64(l.to) for l in ps.line]
    r = [Float64(l.R) for l in ps.line]
    x = [Float64(l.X) for l in ps.line]
    b = [Float64(l.B) for l in ps.line]
    tap = Vector{Complex{Float64}}([l.tap for l in ps.line])
    ysh = Vector{Complex{Float64}}([Complex(b.g,b.b) for b in ps.bus])

    # RTS-96 bus tags are sorted integers. Fortunately.
    tags = sort(unique([from;to]))
    # map bus tags to indices.
    f = tags2indices(from,tags)
    t = tags2indices(to,tags)

    # admittance matrix
    Y = imag(getY(lineOut,f,t,r,x,b,tap,ysh))

    # Renewable indices
    r = convert(Vector{Int64},[w.bus for w in ps.wind])
    Ridx = tags2indices(r,tags)

    # G0, D0, R0
    ps.params.numOfTimeInt = 6

    # G0
    g = [Int64(g.bus) for g in ps.gen]
    gidx = tags2indices(g,tags)
    ginj = [Float64(g.Pinj) for g in ps.gen]
    gfull = zeros(length(ps.bus))
    gfull[gidx] = ginj
    G0 = repmat(gfull,ps.params.numOfTimeInt)

    # D0
    d = [Int64(l.bus) for l in ps.load]
    lidx = tags2indices(d,tags)
    linj = [Float64(l.Pload) for l in ps.load]
    lfull = zeros(length(ps.bus))
    lfull[lidx] = linj
    D0 = repmat(lfull,ps.params.numOfTimeInt)

    # R0
    rinj = [Float64(r.Pinj) for r in ps.wind]
    rfull = zeros(length(ps.bus))
    rfull[Ridx] = rinj
    R0 = repmat(rfull,ps.params.numOfTimeInt)

    # easy translations
    Sb = ps.params.sBase*1e6
    ref = findfirst([Int64(b.ty) for b in ps.bus],3)
    type0lines = find([l.ty for l in ps.line].==0)
    lines = [Tuple{Int64,Int64}((f[i],t[i])) for i in 1:length(f)]
    res = [Float64(l.R) for l in ps.line]
    reac = [Float64(l.X) for l in ps.line]

    # remove type 1 lines from instanton analysis:
    lines = lines[type0lines]
    res = res[type0lines]
    reac = reac[type0lines]

    k = ps.params.participation
    line_lengths = [Float64(l.length) for l in ps.line]
    line_lengths = line_lengths[type0lines]

    # line_conductors needs work. I can't assume all lines are
    # same conductor type any more. I need to use Jon's conductor info
    # to build instances of LineParams and ConductorParams so my code can
    # infer thermal constant values.

    # I effectively need to move thermal model generation outside the loop.
    # This is probably something I should have done before anyway. It makes
    # no sense to penalize my algorithm by forcing it to generate thermal models
    # that in a practical implementation would simply be read in.

    # this one's pretty easy:
    Tamb = maximum(unique([Float64(l.thermalModel.T_amb) for l in ps.line]))

    T_max = [Float64(l.thermalModel.T_max) for l in ps.line]
    T_max = T_max[type0lines]
    delT = [Float64(l.delT) for l in ps.line]
    delT = delT[type0lines]
    # This is a vector, not a scalar. My algorithm should require this anyway.
    # requires a slight change to my code.
    T0 = T_max + delT

    int_length = 300.
    time_values = 0:30.:300
    corr = Array{Float64,2}()

    # construct instance of LineParams for each line
    line_params_vec = Vector{LineParams}()
    for l in ps.line
        push!(line_params_vec,LineParams(
        l.from,
        l.to,
        l.R,
        l.X,
        l.length*1609.34
        ))
    end
    line_params_vec = line_params_vec[type0lines]

    # construct instance of ConductorParams for each line
    conductor_params_vec = Vector{ConductorParams}()
    for l in ps.line
        c = l.conductorModel
        t = l.thermalModel
        push!(conductor_params_vec,ConductorParams(
        c.D*1e-3,
        t.mCp,
        c.I_max,
        c.R_unit,
        t.T_max,
        t.eta_c,
        t.eta_r,
        t.eta_s,
        ))
    end
    conductor_params_vec = conductor_params_vec[type0lines]

    return InstantonInputData(
    Ridx,
    Y,
    G0,
    D0,
    R0,
    Sb,
    ref,
    lines,
    res,
    reac,
    k,
    line_lengths,
    (line_params_vec,conductor_params_vec),
    Tamb,
    T0,
    int_length,
    time_values,
    corr
    )
end

"""
    (Y,Yf,Yt) = getY(lineOut,f,t,r,x,b,tap,ysh)
This function builds admittance matrices.

INPUTS:

*  `lineOut`: [nline x 1, logical] which lines are out of service
*  `f`: nline Vector of from buses (range 1 to nbus)
*  `t`: nline Vector of to buses (range 1 to nbus)
*  `r`: nline Vector of series resistances (pu)
*  `x`: nline Vector of series reactances (pu)
*  `b`: nline Vector of line charging shunt susceptances (pu)
*  `tap`: nline vector of complex tap ratios (pu)
*  `ysh`: nbus vector of complex bus shunt admittances (pu)

for each branch, compute the elements of the branch admittance matrix where
```
  | If |   | Yff  Yft |   | Vf |
  |    | = |          | * |    |
  | It |   | Ytf  Ytt |   | Vt |
```
*Assumes the grid is connected (i.e. no island nodes).*
"""
function getY(
    lineOut::Vector{Bool},
    from::Vector{Int64},
    to::Vector{Int64},
    r::Vector{Float64},
    x::Vector{Float64},
    b::Vector{Float64},
    tap=fill(Complex(1.),length(f))::Vector{Complex{Float64}},
    ysh=fill(0.im,length(unique([f;t])))::Vector{Complex{Float64}}
    )
    # remap bus tags to 1:nbus indices
    tags = sort(unique([from;to]))
    f = tags2indices(from,tags)
    t = tags2indices(to,tags)

    inService = find(!lineOut)
    tap = tap[inService]
    fis = f[inService]
    tis = t[inService]

    Ys = (1./(r + x*im))[inService]
    Bc = b[inService]
    Ytt = Ys + Bc*im/2
    Yff = Ytt./(tap.*conj(tap))
    Yft = -Ys./conj(tap)
    Ytf = -Ys./tap

    # build connection matrices
    nbus,nline = length(ysh),length(lineOut)
    # connection matrix for line and from buses
    Cf = sparse(1:nline,f,ones(nline),nline,nbus)
    # connection matrix for line and to buses
    Ct = sparse(1:nline,t,ones(nline),nline,nbus)

    # build Yf and Yt such that Yf*V is the vector of complex branch currents injected.
    i = [inService;inService]
    Yf = sparse(i,[fis;tis],[Yff;Yft],nline,nbus)
    Yt = sparse(i,[fis;tis],[Ytf;Ytt],nline,nbus)

    # build Ybus
    Y = Cf'*Yf + Ct'*Yt + sparse(1:nbus,1:nbus,ysh,nbus,nbus)
    return Y #,Yf,Yt
end

"""
Given a vector `v` containing bus tags, and a complete and sorted
vector of bus tags, return a version of `v` mapped to bus indices.

Example: if there are three buses in the network called 101, 202, and 303,
one would let `tags = [101;202;303]`. Corresponding indices are simply `[1;2;3]`. So this function would convert a vector `v = [202;303]` into `vidx = [2;3]`.
"""
function tags2indices(
    v::Vector{Int64},
    tags::Vector{Int64}
    )
    indices = collect(1:length(tags))
    vidx = Vector{Int64}()
    for vi in v
        push!(vidx,indices[findfirst(tags,vi)])
    end
    return vidx
end
