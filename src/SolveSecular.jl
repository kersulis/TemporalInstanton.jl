"""
Solve the explicit secular equation
by returning the smallest `v` satisfying

f(v) = w'\*w - s2 = 0,

where `w = d./(D - v)`.

Inputs:

* `d` is a vector with numerator values
* `D` is a vector of poles (no need to sort them)
* `s2` is the target constant; we return solutions s(v) = s2.

Output:

* `val` is the optimal `v` satisfying s(v) = c
* `vecs` is a vector of vectors `w` corresponding to `vals`
"""
function solvesecular(
    d::Vector{Float64},
    D::Vector{Float64},
    s2::Float64
    )

    # sort eigenvalues by magnitude:
    p = sortperm(D)
    D = D[p]
    d = d[p]

    # explicit secular eq
    f(v::Float64) = sumabs2(d./(D-v)) - s2

    # derivative:
    fprime(v::Float64) = sum((2*d.^2)./((D-v).^3))

    # vector in terms of v
    w(v::Float64) = float(d./(v - D))[invperm(p)]

    # index of first nonzero of d
    k = minimum(find(d))

    iterates = Vector{Float64}()
    v0 = D[k] - abs(d[k])/sqrt(s2)
    push!(iterates,v0)

    stop = false
    while !stop
        v = iterates[end]
        vnew = v - 2*((f(v) + s2)/fprime(v))*(sqrt((f(v)+s2)/s2) - 1)
        push!(iterates,vnew)
        stop = vnew >= v
    end
    v = iterates[end]
    return v, w(v)
end
