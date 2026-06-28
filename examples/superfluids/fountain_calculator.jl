module fc
const T = 1.65
const S = 335.0 
const d = 1.55e-3
const rho = 145.2352 
const g = 9.8
const A = pi*d*d/4

function power2speed(W)
    return W/(T*S*rho*A)
end

function power2height(W)
    v = power2speed(W)
    return 0.5*v*v/g
end

function height2power(h)
    v = sqrt(2*h*g)
    return v*T*S*rho*A
end

end