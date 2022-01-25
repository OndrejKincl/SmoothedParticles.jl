module FixPA
module FixPA

include("../src/SPHLib.jl")
using .SPHLib

export rev_add, rev_sub, static

const FixPA_eps = 1/(2^30)

function nom(x::Float64)::Int64
    return Int64(round(x/FixPA_eps))
end

function rev_add(x::Float64, y::Float64)::Float64
    return FixPA_eps*(nom(x) + nom(y))
end

function rev_add(x::RealVector, y::RealVector)::RealVector
    v1 = FixPA_eps*(nom(x[1]) + nom(y[1]))
    v2 = FixPA_eps*(nom(x[2]) + nom(y[2]))
    v3 = FixPA_eps*(nom(x[3]) + nom(y[3]))
    return RealVector(v1,v2,v3)
end

end