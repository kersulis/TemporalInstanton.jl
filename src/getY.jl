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
    indices = collect(1:length(tags))
    f = Vector{Int64}()
    for fi in from
        push!(f,indices[findfirst(tags,fi)])
    end
    t = Vector{Int64}()
    for ti in to
        push!(t,indices[findfirst(tags,ti)])
    end

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
