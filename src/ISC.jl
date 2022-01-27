using SparseArrays

function add_element!(I::Vector{Int64}, J::Vector{Int64}, V::Vector{Float64}, i::Int64, j::Int64, v::Float64)
    if v == 0.0
        return
    end
    push!(I, i)
    push!(J, j)
    push!(V, v)
end

function lhs(sys::ParticleSystem, rDker::Function)::SparseMatrixCSC{Float64}
	N = length(sys.particles)
    I = Int64[]
    J = Int64[]
    V = Float64[]
    for i in 1:N
        p = sys.particles[i]
        S1 = SPHLib.sum(sys, (p,q,r) -> -2*rDker(sys.h,r)*(p.x[1]-q.x[1]), p)
        S2 = SPHLib.sum(sys, (p,q,r) -> -2*rDker(sys.h,r)*(p.x[2]-q.x[2]), p)
        S = RealVector(S1, S2, 0.)
        for j in 1:N
            q = sys.particles[j]
            r = dist(p,q)
            if r > sys.h
                continue
            end
            tmp0 = Float64(p==q)
            tmp1 = rDker(sys.h,r)*(p.x - q.x)
            tmp2 = 0.5*Float64(p==q)*S
            #UL block = I
            add_element!(I,J,V,i+0N,j+0N,tmp0)
            add_element!(I,J,V,i+1N,j+1N,tmp0)
            #UR block = Grad
            add_element!(I,J,V,i+0N,j+2N, tmp1[1]-tmp2[1])
            add_element!(I,J,V,i+1N,j+2N, tmp1[2]-tmp2[2])
            #DL block = Div
            add_element!(I,J,V,i+2N,j+0N, tmp1[1]+tmp2[1])
            add_element!(I,J,V,i+2N,j+1N, tmp1[2]+tmp2[2])
        end
    end
    return sparse(I, J, V, 3N, 3N)
end

function rhs(sys::ParticleSystem, ker::Function, rho0::Float64)
    N = length(sys.particles)
    b = zeros(3N)
    for i in 1:N
        p = sys.particles[i]
        b[i+2N] = SPHLib.sum(sys, (p,q,r) -> ker(sys.h,r), p) - rho0
	end
    return b
end

function renormalize!(sys::ParticleSystem, dr::Float64; 
    tol = 1e-6, max_steps = 10,
    ker::Function = wendland2, rDker::Function = rDwendland2)
    @info("renormalizing, please wait...")
    rho0 = 1/dr^2
    for p in sys.particles
        p.x += 0.3*dr*RealVector(rand() - 1.0, rand() - 1.0, 0.)
    end
    for _ in 1:max_steps
        create_cell_list!(sys)
        b = rhs(sys, ker, rho0)
        err = maximum(x -> abs(x), b)
        @show err
        if err < tol
            break
        end
        A = lhs(sys, rDker)
        y = A\b
        N = length(sys.particles)
        for i in 1:N
            sys.particles[i].x += RealVector(y[i], y[i+N], 0.)
        end
	end
end
