using MatpowerCases, MAT
include("PowerSystem.jl")
include("solvetmpinst-collab.jl")

function network_data(casename::String)
	if !in(casename,casenames())
		error("Argument is not a Matpower case name or path to .mat file")
	end
	mpc = loadcase(casename,describe=false)
	tags = Vector{Int64}(mpc["bus"][:,1])
	genBus = tags2indices(Vector{Int64}(mpc["gen"][:,1]),tags)
	try
		Sb = mpc["baseMVA"]
	catch
		warn("Base MVA missing from mpc data. Using 100 MVA default.")
		Sb = 100.0
	end
	f = round(Int64,mpc["branch"][:,1]) # "from bus" ...
	f = tags2indices(f,tags)
	t = round(Int64,mpc["branch"][:,2]) # ... "to bus"
	t = tags2indices(t,tags)
	r = mpc["branch"][:,3]              # resistance, pu
	x = mpc["branch"][:,4]              # reactance, pu
	b = mpc["branch"][:,5]              # susceptance, pu
	Y = createY(f,t,x)
	G_long = mpc["gen"][:,2]
	G = zeros(length(tags))
	for i in 1:length(tags)
		G[i] = sum(G_long[find(genBus.==i)])/Sb
	end
	D = mpc["bus"][:,3]./Sb
	Sb = Sb*1e6 # convert from MW to W
	ref = 1
	# Allow each generator to participate equally in droop response.
	k = Vector{Float64}()
	ngen = length(find(G))
	for i = 1:length(G)
		if G[i] != 0
		    push!(k, 1/ngen)
		else
		    push!(k,0)
		end
	end
	n = Dict(
	"f"=>f,
	"t"=>t,
	"r"=>r,
	"x"=>x,
	"b"=>b,
	"G"=>G,
	"k"=>k,
	"D"=>D,
	"Y"=>Y,
	"Sb"=>Sb)
	return n
end

function network_data(ps::PowerSystem)
	from = [Int64(l.from) for l in ps.line]
	to = [Int64(l.to) for l in ps.line]
	r = [Float64(l.R) for l in ps.line]
	x = [Float64(l.X) for l in ps.line]
	b = [Float64(l.B) for l in ps.line]
	lineOut = ps.control.lineOut
	tap = [Complex(l.tap) for l in ps.line]
	ysh = [Complex(b.b) for b in ps.bus]

	# RTS-96 bus tags are sorted integers. Fortunately.
	tags = sort(unique([from;to]))
	# map bus tags to indices.
	f = tags2indices(from,tags)
	t = tags2indices(to,tags)

	# admittance matrix
	Y = imag(getY(lineOut,f,t,r,x,b,tap,ysh))
	# G0
	g = [Int64(g.bus) for g in ps.gen]
	gidx = tags2indices(g,tags)
	ginj = [Float64(g.Pinj) for g in ps.gen]
	G = zeros(length(ps.bus))
	G[gidx] = ginj
	# D0
	d = [Int64(l.bus) for l in ps.load]
	lidx = tags2indices(d,tags)
	linj = [Float64(l.Pload) for l in ps.load]
	D = zeros(length(ps.bus))
	D[lidx] = linj
	# R0
	r = convert(Vector{Int64},[w.bus for w in ps.wind])
	Ridx = tags2indices(r,tags)
	rinj = [Float64(r.Pinj) for r in ps.wind]
	R = zeros(length(ps.bus))
	R[Ridx] = rinj
	D -= R
	# easy translations
	Sb = ps.params.sBase*1e6
	ref = findfirst([Int64(b.ty) for b in ps.bus],3)
	k = ps.params.participation

	n = Dict(
	"f"=>f,
	"t"=>t,
	"r"=>r,
	"x"=>x,
	"b"=>b,
	"G"=>G,
	"k"=>k,
	"D"=>D,
	"Y"=>Y,
	"Sb"=>Sb)
	return n
end

function thermal_data(casename::String)
	
end

function thermal_data(ps::PowerSystem)
	D = [Float64(l.conductorModel.D*1e-3) for l in ps.line]
	Ilim = [Float64(l.conductorModel.I_max) for l in ps.line]
	mCp = [Float64(l.thermalModel.mCp) for l in ps.line]
	eta_c = [Float64(l.thermalModel.eta_c) for l in ps.line]
	eta_r = [Float64(l.thermalModel.eta_r) for l in ps.line]
	eta_s = [Float64(l.thermalModel.eta_s) for l in ps.line]
	length = [Float64(l.length) for l in ps.line]
	t = Dict(
	"D"=>D,
	"Ilim"=>Ilim,
	"mCp"=>mCp,
	"eta_c"=>eta_c,
	"eta_r"=>eta_r,
	"eta_s"=>eta_s,
	"length"=>length
	)
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

"""
    createY(f,t,x [,r,b]) -> Y
Create an admittance matrix for AC power flow.
All inputs are real.

* `f`,`t`: vectors encoding all lines (fi,ti)
* `x`: per-unit reactance xi for all lines
* `r`: per-unit resistance ri for all lines
* `b`: per-unit susceptance bi for all lines
"""
function createY(
    f::Vector{Int64},
    t::Vector{Int64},
    x::Vector{Float64},
    r=0.0::Union{Vector{Float64},Float64},
    b=0.0::Union{Vector{Float64},Float64}
    )
    z = r + x*1im
    y = 1./z
    b = b*1im
    Y = sparse([f; t; t; f],[t; f; t; f],[-y; -y; y + b./2; y + b./2])

    # for DC power flow, we typically want a matrix with real entries:
    if r == 0
        return imag(Y)
    else
        return Y
    end
end

function thermal_data(n)

end
