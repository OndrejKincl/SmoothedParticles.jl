module test_IO
using Test
include("../src/SmoothedParticles.jl")
using .SmoothedParticles
using Parameters

const dr = 1/100
const h = 2.0*dr

@with_kw mutable struct Particle <: AbstractParticle
    x::RealVector
    s::Float64    = 0.    # scalar variable
    v::RealVector = VEC0  # vector variable
    M::RealMatrix = MAT0  # matrix variable
end

# make some random data
function get_vars(x::RealVector)::Tuple{Float64, RealVector, RealMatrix}
    s = x[2]
    v = x[2]*VECX - x[1]*VECY
    M = x[1]*RealMatrix(0.,1.,2.,3.,4.,5.,6.,7.,8.)
    return (s,v,M)
end

# generate particle system
function make_sys()::ParticleSystem
    dom = Circle(0., 0., 1.)
    sys = ParticleSystem(Particle, dom, 0.1)
    return sys
end

function main()

    sys = make_sys()
    grid = Grid(dr, :hexagonal)
    generate_particles!(sys, grid, Circle(0., 0., 1.), x -> Particle(x=x))
    
    # assign some data
    for p in sys.particles
        (p.s, p.v, p.M) = get_vars(p.x)
    end

    @testset "save data to vtk" begin
        out = new_pvd_file("test_IO")
        save_frame!(out, sys, :s, :v, :M)
        save_pvd_file(out)
        @test ispath("test_IO/frame0.vtp")
        @test ispath("test_IO/result.pvd")
    end

    @testset "read data from vtk" begin
        sys_ = make_sys()
        import_particles!(sys_, "test_IO/frame0.vtp", x -> Particle(x=x))
        @test length(sys.particles) == length(sys_.particles)
        @test all((p.s, p.v, p.M) == get_vars(p.x) for p in sys_.particles)
        # import even more particles 
        import_particles!(sys_, "test_IO/frame0.vtp", x -> Particle(x=x))
        @test 2*length(sys.particles) == length(sys_.particles)
        @test all((p.s, p.v, p.M) == get_vars(p.x) for p in sys_.particles)
    end

    # clean up
    rm("test_IO/", recursive=true)
    @info("deleted folder: test_IO")
end


end