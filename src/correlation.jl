"Similar to MATLAB's polyfit; from Rosetta Code"
function polyfit(x, y, n)
  A = [ float(x[i])^p for i = 1:length(x), p = 0:n ]
  A \ y
end

"Generate n random coordinates in a square given by (lower,upper)"
function returncoordinates(n,lower,upper)
    c = rand(collect(lower:upper),2*n)
    return [(c[i],c[i+1]) for i in 1:n]
end

"Turn a set of geographic points into a distance matrix."
function distancematrix(pts)
    n = length(pts)
    D = zeros(n,n)
    for i in 1:n
        for j in 1:n
            D[i,j] = hypot(pts[i][1] - pts[j][1],pts[i][2] - pts[j][2])
        end
    end
    return D
end

"Use freris2008 Figure 2.10 to map distance to correlation coefficient."
function dist2cor(dist)
    x = [0;200;400;600;800;1000;1200.]
    y = [0.88;0.65;0.45;0.30;0.20;0.12;0.08]
    coef = polyfit(x,y,2)
    f(x) = coef[1] + coef[2].*x + coef[3].*x.^2
    D = f(dist)
    for i in 1:size(D,1)
        D[i,i] = 1.0
    end
    return D
end
