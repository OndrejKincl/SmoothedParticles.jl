const EPS = 1.4901161193847656e-8

#STRUCTURAL KERNELS
#------------------

@fastmath function wendland2h(h::Float64, r::Float64)::Float64
    x = r/h
    return x < 1.0 ? 14.0*(1.0 - x)^3*(14.0*x^2 - 3.0*x - 1.0)/(pi*h^2) : 0.
end

@fastmath function rDwendland2h(h::Float64, r::Float64)::Float64
    x = r/h
    return x < 1.0 ? 140.0*(1.0 - x)^2*(4.0 - 7.0*x)/(pi*h^4) : 0.
end

@fastmath function wendland3h(h::Float64, r::Float64)::Float64
    x = r/h
    return x < 1.0 ? 21.0*(1.0 - x)^3*(32.0*x^2 - 9.0*x - 3.0)/(2.0*pi*h^3) : 0.
end

@fastmath function rDwendland3h(h::Float64, r::Float64)::Float64
    x = r/h
    return x < 1.0 ? 210.0*(1.0 - x)^2*(5.0 - 8.0*x)/(pi*h^5) : 0.
end

#ALGEBRAIC TOOLS
#---------------

@inbounds function outer(x::RealVector, y::RealVector)::RealMatrix
    return RealMatrix(
        x[1]*y[1], x[2]*y[1], x[3]*y[1], 
        x[1]*y[2], x[2]*y[2], x[3]*y[2],
        x[1]*y[3], x[2]*y[3], x[3]*y[3]
    )
end

@inbounds function trace(G::RealMatrix)::Float64
    return G[1,1] + G[2,2] + G[3,3]
end

@inbounds function dev(G::RealMatrix)::RealMatrix
    return G - 1/3*trace(G)*MAT1
end

@inbounds function subinv(A::RealMatrix)::RealMatrix
    idet = 1.0/(A[1,1]*A[2,2] - A[2,1]*A[1,2])
    return RealMatrix(
        +idet*A[2,2], -idet*A[2,1], 0., 
        -idet*A[1,2], +idet*A[1,1], 0.,
        0., 0., 0.
    )
end

#DATA SAVING
#-----------

function vec2string(a::Vector)::String
    out = ""
    for i in 1:length(a)-1
        out = out*string(a[i])*","
    end
    if length(a) > 0
        out = out*string(a[end])
    end
    out = out*"\n"
end

#ARGMIN FUNCTION
#---------------

function find_minimizer(f::Function, sys::ParticleSystem)::Particle
    p = sys.particles[1]
    pval = f(p)
    for q in sys.particles
          qval = f(q)
        if qval < pval
              p = q
              pval = qval
          end
    end
    return p
 end

 #SAVE INVERSE
 #------------

 @fastmath function saveinv(x::Float64)::Float64
    return x/(x*x + EPS)
 end