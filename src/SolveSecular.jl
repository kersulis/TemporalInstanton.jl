"""
Return vector w satisfying
w = (num)./(v - poles)
"""
function find_w(
    v::Float64,
    D::Array{Float64,2},
    d::Vector{Float64}
    )
    if v == 0.0
        w = float([-d[i]/(D[i,i]) for i in 1:length(d)])
    else
        w = float([d[i]/(v - D[i,i]) for i in 1:length(d)])
    end
    return w
end

"""
Solve the secular equation by returning `v` satisfying
s(v) = w'\*w = c,

where `w = num./(v - poles)`

Inputs:

* `num` is a diagonal matrix with numerator
* `poles` is a vector of poles (no need to sort them)
* `c` is the target constant; we return solutions s(v) = c.

Output:

* `vals` is a vector of solutions `v` satisfying s(v) = c
* `vecs` is a vector of vectors `w` corresponding to `vals`
"""
function solvesecular(
    num::Vector{Float64},
    poles::Vector{Float64},
    c::Float64
    )
    s(vv::Vector{Float64}) = [sumabs2(num./(v - poles)) for v in vv]
    s(v::Float64) = sumabs2(num./(v - poles))
    w(vv::Vector{Float64}) = [float(num./(v - poles)) for v in vv]
    w(v::Float64) = float(num./(v - poles))
    tiny = c/1e8

    vals = Vector{Float64}()
    sortedpoles = sort(poles)

    for i in 1:length(sortedpoles)
        upper = sortedpoles[i]

        # head left
        dist = finddist(i,sortedpoles,true)
        lower = upper - dist/2
        intersect = secularintersect(c,upper,lower,s,tiny)

        isinf(intersect) && continue

        # there must be at least one solution for leftmost pole
        if i == 1
            ex = 0
            while isnan(intersect)
                ex += 1
                if ex == 10
                    intersect = NaN
                    break
                end
                # warn("expanding $ex")
                lower -= 100*dist
                intersect = secularintersect(c,upper,lower,s,tiny)
            end
        end
        push!(vals,intersect)

        # head right
        dist = finddist(i,sortedpoles,false)
        lower = upper + dist/2
        intersect = secularintersect(c,upper,lower,s,tiny)
        # there must be at least one solution for rightmost pole
        if i == length(sortedpoles)
            ex = 0
            while isnan(intersect)
                ex += 1
                if ex == 10
                    intersect = NaN
                    break
                end
                # warn("expanding $ex")
                lower += 100*dist
                intersect = secularintersect(c,upper,lower,s,tiny)
            end
        end
        push!(vals,intersect)
    end
    filter!(x -> !isnan(x),vals)
    return vals, w(vals)
end

"""
Given current pole index, vector of sorted poles, and desired search direction (`left` or right==`!left`), return appropriate initial search distance for `secularintersect`.
"""
function finddist(
    i::Int64,
    sortedpoles::Vector{Float64},
    left::Bool
    )
    p = sortedpoles[i]
    if length(sortedpoles) == 1
        dist = p
    elseif i == 1 && left
        dist = abs(p - sortedpoles[2])
    elseif i == length(sortedpoles) && !left
        dist = abs(p - sortedpoles[end-1])
    else
        dist = left ? abs(p - sortedpoles[i-1]) : abs(p - sortedpoles[i+1])
    end
    return dist
end

"""
Given c, upper bound, lower bound, function s, and a "tiny" number, return the value v where s(v) = c (if it exists).
"""
function secularintersect(
    c::Float64,
    upper::Float64,
    lower::Float64,
    s::Function,
    tiny::Float64
    )
    if c == 0.0
        error("c > 0 required; s(v) = 0 has no solution")
    end
    if s(lower) > 1.0/tiny
        warn("c is too small; s(v) = c has no reasonable solution")
        return Inf
    end
    ubound = upper
    olddiff = 0.0
    v = mean([upper;lower])
    currentdiff = s(v) - c

    while abs(currentdiff) > tiny
        if abs(currentdiff - olddiff) < tiny
            #||            (currentdiff > 0 && olddiff > 0 &&  currentdiff > olddiff)
            return NaN
        end
        if currentdiff > 0
            upper = v
        else
            lower = v
        end
        v = mean([upper;lower])
        olddiff = currentdiff
        currentdiff = s(v) - c
    end
    return v
end
