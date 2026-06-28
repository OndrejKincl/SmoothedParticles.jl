module invertibility

const dr = 1/100
const h = 3*dr
const c = 20.0

using SmoothedParticles
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using IncompleteLU

mutable struct Particle <: AbstractParticle
    x::RealVector
end


function wab(p::Particle, q::Particle, r::Float64)::Float64
	return h^2*wendland2(h,r)
end

function main()
    dt = 0.1*h/c
    grid = Grid(dr, :hexagonal)
    geo = Rectangle(0., 0., 1., 1.)
    dom = Rectangle(-0.2, -0.2, 1.2, 1.2)
    sys = ParticleSystem(Particle, dom, h)
    generate_particles!(sys, grid, geo, x -> Particle(x))
    push!(sys.particles, Particle(sys.particles[1].x))
    create_cell_list!(sys)
    A = assemble_matrix(sys, wab)
    (N, ) = size(A)
    b = ones(N)
    #@show A
    #=
    @time begin
        @info "mldiv"
        x1 = A\b
        #@show x
    end
    =#
    @time begin
        @show dt
        @info "cg"
        #Pl = ilu(A, τ = 0.01)
        Pl = Identity()
        (x2, ch) = cg(A,b; Pl = Pl, log = true)
        @show ch.isconverged
        #@show ch[:resnorm]
        #@show x
    end
    #@show x1-x2
end

end