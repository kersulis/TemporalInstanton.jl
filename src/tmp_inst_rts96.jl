# Begin by importing RTS-96 data as with previous project.

using MAT

type powerSystemData
    Sb
    f
    t
    r
    x
    b
    Y
    bustype
    Gp
    Gq
    Dp
    Dq
    Rp
    Rq
    Pmax
    Pmin
    Qmax
    Qmin
    Plim
    Vg
    Vceiling
    Vfloor
    busIdx
    N
    Nr
    Ng
    k
end

function createY(f,t,r,x,b,s)
    # Create an admittance matrix from three vectors: from, to, and adm. value
    # Note: 'f' and 't' must be integer vectors.
    # DC if s == true

    if s == true
        y = 1./x
        return sparse([f; t; t; f],[t; f; t; f],[-y; -y; y; y])
    else
        G = 1./r
        G[G.==Inf] = 0
        B = 1./x
        y = complex(G,B)
        return sparse([f; t; t; f],[t; f; t; f],[-y; -y; y + b./2; y + b./2])
    end
end

function unpack_psDL(psDL)
# code used to unpack psDL type instance:
return (psDL.Sb,
    psDL.f,
    psDL.t,
    psDL.r,
    psDL.x,
    psDL.b,
    psDL.Y,
    psDL.bustype,
    psDL.Gp,
    psDL.Gq,
    psDL.Dp,
    psDL.Dq,
    psDL.Rp,
    psDL.Rq,
    psDL.Pmax,
    psDL.Pmin,
    psDL.Qmax,
    psDL.Qmin,
    psDL.Plim,
    psDL.Vg,
    psDL.Vceiling,
    psDL.Vfloor,
    psDL.busIdx,
    psDL.N,
    psDL.Nr,
    psDL.Ng,
    psDL.k)
end

include("readRTS96Data.jl")
# psData = psDataLoad()