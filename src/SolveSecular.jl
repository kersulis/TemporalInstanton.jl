""" Solve the secular equation via binary search.
"""
function solve_secular(D,d,c)
    tinynumber = 1e-8
    solutions = Float64[]
    vectors = Array(Vector{Float64},0)
    poles = sort(unique(round(diag(D),10)))

    # Each diagonal element is a pole.
    for i in 1:length(poles)

        # Head left first:
        high = poles[i]
        if length(poles) == 1
            low = high - high
        elseif i == 1
            low = high - abs(poles[i] - poles[i+1])
        else
            low = high - abs(poles[i] - poles[i-1])/2
        end

        # Initialize v:
        v = (high + low)/2
        w = find_w(v,D,d)
        currentdiff = (w'*w)[1] - c
        olddiff = 0
        stall = false
        while abs(currentdiff) > tinynumber
            if currentdiff == olddiff
                stall = true
                break
            end
            if diff > 0
                high = v
            else
                low = v
            end
            v = (high + low)/2
            w = find_w(v,D,d)
            olddiff = currentdiff
            currentdiff = (w'*w)[1] - c
        end
        if !stall
            push!(solutions,v)
            push!(vectors,w)
        end

        # Now head right:
        high = poles[i]
        if length(poles) == 1
            low = high + high
        elseif i == length(poles)
            low = high + abs(poles[i] - poles[i-1])
        else
            low = high + abs(poles[i] - poles[i+1])/2
        end

        v = (high + low)/2
        w = find_w(v,D,d)
        currentdiff = (w'*w)[1] - c
        olddiff = 0
        stall = false
        while abs(currentdiff) > tinynumber
            if currentdiff == olddiff
                stall = true
                break
            end
            if currentdiff > 0
                high = v
            else
                low = v
            end
            v = (high + low)/2
            w = find_w(v,D,d)
            olddiff = currentdiff
            currentdiff = (w'*w)[1] - c
        end
        if !stall
            push!(solutions,v)
            push!(vectors,w)
        end
    end
    return solutions,vectors
end
