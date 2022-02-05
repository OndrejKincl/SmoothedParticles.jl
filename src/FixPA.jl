"""
    FixPA_eps

    Epsilon in the fixed-point arithmetic
"""
const FixPA_eps = 1/(2^30)

"""
    nom(x::Float64)::Int64

Value of ``x`` in the fixed-point arithmetic.
"""
function nom(x::Float64)::Int64
    return Int64(round(x/FixPA_eps))
end

"""
    rev_add(x::Float64, y::Float64)::Float64    

Reversible addition of two floats.
"""
function rev_add(x::Float64, y::Float64)::Float64
    return FixPA_eps*(nom(x) + nom(y))
end


"""
    rev_add(x::RealVector, y::RealVector)::RealVector

Reversible addition of two vectors.
"""
function rev_add(x::RealVector, y::RealVector)::RealVector
    v1 = FixPA_eps*(nom(x[1]) + nom(y[1]))
    v2 = FixPA_eps*(nom(x[2]) + nom(y[2]))
    v3 = FixPA_eps*(nom(x[3]) + nom(y[3]))
    return RealVector(v1,v2,v3)
end