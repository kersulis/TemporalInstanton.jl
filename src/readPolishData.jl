using MatpowerCases

function load_polish_data()
    mpc = loadcase("case3120sp",describe=false)
    bus_i = mpc["bus"][:,1]
    genBus = mpc["gen"][:,1]
    Sb = mpc["baseMVA"]
    busType = mpc["bus"][:,2]
    Gp_long = mpc["gen"][:,2]
    #Gq_long = mpc["gen"][:,3]

    #Rp = []
    #Rq = []

    Dp = mpc["bus"][:,3]./Sb
    #Dq = mpc["bus"][:,4]./Sb

    #Pmax = mpc["gen"][:,9]./Sb
    #Pmin = mpc["gen"][:,10]./Sb
    #Qmax = mpc["gen"][:,4]./Sb
    #Qmin = mpc["gen"][:,4]./Sb

    Vm = mpc["bus"][:,8]
    Va = zeros(length(Vm))

    f = round(Int64,mpc["branch"][:,1]) # "from bus" ...
    t = round(Int64,mpc["branch"][:,2]) # ... "to bus"
    r = mpc["branch"][:,3]              # resistance, pu
    x = mpc["branch"][:,4]              # reactance, pu
    b = mpc["branch"][:,5]              # susceptance, pu

    Gp = zeros(length(bus_i))
    for i in bus_i
        Gp[i] = sum(Gp_long[find(genBus.==i)])/Sb
    end

    # convert generators smaller than 82 MW into wind farms
    Rp = zeros(length(Gp))
    for i in 1:length(Gp)
        if Gp[i] < 0.82
            Rp[i] = Gp[i]
            Gp[i] = 0
        end
    end

    #N = length(bus_i)
    #Nr = []
    #Ng = length(find(Gp))

    Y = createY(f,t,r,x,b,true)

    Ridx = find(Rp)

    Sb = Sb*1e6 # convert from MW to W

    ref = 1

    lines = [(f[i],t[i]) for i in 1:length(f)]

    res = r
    reac = x

    # Allow each generator to participate equally in droop response.
    # Note: this only applies to analysis types with droop response!
    k = Float64[]
    for i = 1:length(Gp)
        if Gp[i] != 0
            push!(k, 1/length(find(Gp)))
        else
            push!(k,0)
        end
    end

    # return  Sb, f, t, r, x, b, Y, bustype,
    #         Gp, Gq, Dp, Dq, Rp, Rq,
    #         Pmax, Pmin, Qmax, Qmin,
    #         Plim, Vg, Vceiling, Vfloor,
    #         bus_i, N, Nr, Ng

    # use RTS-96 line lengths to generate similar line lengths
    line_lengths = load("../data/polish_line_lengths.jld","line_lengths")

    # temporary (re-use rts-96 line conductor parameters)
    line_conductors = fill("waxwing",length(line_lengths))

    return  Ridx, Y,
            Gp, Dp, Rp,
            Sb, ref, lines,
            res, reac, k,
            line_lengths, line_conductors
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
