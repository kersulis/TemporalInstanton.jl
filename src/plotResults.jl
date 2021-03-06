type temporalInstantonResults
	x
	θ
	α
end

function generateDot(name,
	psData,
	tmpInstRes,
    idx,tidx,
    G0,R0,D0,
    dim)

	x = tmpInstRes.x
	θ = tmpInstRes.θ
	α = tmpInstRes.α
	n = length(psData.Gp)

    renewableInj = zeros(n)
    renewableInj[find(psData.Rp)] = x[idx][tidx]
    renewableInj += R0[(tidx-1)*n+1:tidx*n]

    busInj = (G0+ - D0)[(tidx-1)*n+1:tidx*n] + α[idx][tidx]*psData.k + renewableInj

    lineFlow = Float64[]
    f = psData.f
    t = psData.t
    x = psData.x
    for i = 1:length(psData.f)
        push!(lineFlow, (θ[idx][tidx][f[i]]-θ[idx][tidx][t[i]])/x[i])
    end

    writeDot(name,
        psData.busIdx,
        busInj,
        renewableInj,
        psData.f,
        psData.t,
        lineFlow,
        psData.Plim,
    	dim
		)
    #return busInj,renewableInj
end

function writeDot(name, busIdx, busInj, renGen, f, t, lineFlow, lineLim, size=(11,14))
    # This function generates a graph that richly expresses the RTS96 system state.
    # name              a name for the graph and resulting dot file
    # busIdx            bus names (could be text) in order of busInj
    # busInj            injection at each bus
    # renGen            renewable generation at each bus (0 for non-wind buses)
    # f                 "from" node for each line
    # t                 "to" node for each line
    # lineFlow          flow on each line
    # lineLim           limit for each line
    # size              size of graphical output

    busInj = round(busInj,2)
    lineFlow = round(lineFlow,2)

    font = "Helvetica"

    # Open the dot file, overwriting anything already there:
    dotfile = open("$(name).dot","w")

    # Begin writing the dot file:
    write(dotfile, "digraph $(name) {\nnewrank=true;\n")

    # Set graph properties:
    write(dotfile,
    "graph [fontname=\"$(font)\", tooltip=\" \", overlap=false, size=\"$(size[1]),$(size[2])\", ratio=fill, orientation=\"portrait\",layout=dot];\n")

    # Set default node properties:
    write(dotfile, "node [fontname=\"$(font)\", shape=square, style=filled, fontsize=25, color=\"#cccccc\"];\n")

    # Set default edge properties:
    write(dotfile, "edge [fontname=\"$(font)\", penwidth=4, fontsize=25];\n")

    # Write bus data to dot file:
    write(dotfile,
    "subgraph cluster_a1 {\nlabel=\"Area 1: High Wind\";\nfontcolor=\"#000000\";\nfontname=\"$(font)\";\nfontsize=35;\ncolor=\"#ffffff\";\nlabeljust=\"c\";\n")

    for i = 1:24
        write(dotfile,
        "$(i) [label=$(int(busIdx[i])), tooltip=\"Inj = $(busInj[i])\"") # bus label and tooltip

        # Represent renewable nodes with blue circles:
        if union(find(renGen),i) == find(renGen)
            write(dotfile, ", shape=circle, color=\"#CCE6FF\"")
        end

        write(dotfile, "];\n")
    end
    write(dotfile, "}\n")

    write(dotfile,
    "subgraph cluster_a2 {\nlabel=\"Area 2: Moderate Wind\";\nfontcolor=\"#000000\";\nfontname=\"$(font)\";\nfontsize=35;\ncolor=\"#ffffff\";\nlabeljust=\"l\";\n")

    for i = 25:48
        write(dotfile,
        "$(i) [label=$(int(busIdx[i])), tooltip=\"Inj = $(busInj[i])\"") # bus label and tooltip

        # Represent renewable nodes with blue circles:
        if union(find(renGen),i) == find(renGen)
            write(dotfile, ", shape=circle, color=\"#CCE6FF\"")
        end

        write(dotfile, "];\n")
    end
    write(dotfile, "}\n")

    write(dotfile,
    "subgraph cluster_a3 {\nlabel=\"Area 3: Low Wind\";\nfontcolor=\"#000000\";\nfontname=\"$(font)\";\nfontsize=35;\ncolor=\"#ffffff\";\nlabeljust=\"r\";\n")
    for i = 49:length(busIdx)
        write(dotfile,
        "$(i) [label=$(int(busIdx[i])), tooltip=\"Inj = $(busInj[i])\"") # bus label and tooltip

        # Represent renewable nodes with blue circles:
        if union(find(renGen),i) == find(renGen)
            write(dotfile, ", shape=circle, color=\"#CCE6FF\"")
        end

        write(dotfile, "];\n")
    end
    write(dotfile, "}\n")

    # Write line data to file:

    for i = 1:length(f)

        normStrain = abs(lineFlow[i])/lineLim[i] # normalized strain on line i

        # Use flow direction to determine arrow direction,
        # label with flow,
        # and color according to strain
        if lineFlow[i] > 0
            write(dotfile,
            "$(f[i]) -> $(t[i]) [label=$(lineFlow[i])")
        else
            write(dotfile,
            "$(t[i]) -> $(f[i]) [label=$(-lineFlow[i])")
        end
        write(dotfile,
        ", tooltip=\" \", labeltooltip=\"Flow = $(int(normStrain*100))%\", color=\"$(round((1 - normStrain)/3,3)) 1.000 0.700\"];\n")
    end

    # Cap off the dot file and close it:
    write(dotfile, "}")
    close(dotfile)

    #println("$(name).dot generated.\nUse \";dot -Tsvg $(name).dot -o $(name).svg\" to create an SVG.")
end
