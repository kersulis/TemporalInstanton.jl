# Begin by importing RTS-96 data as with previous project.
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
