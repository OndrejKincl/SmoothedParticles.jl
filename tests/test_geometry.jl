module test_geometry
include("../src/SmoothedParticles.jl")
using .SmoothedParticles
using Test
#using Images

const rtol_2d = 0.01
const rtol_3d = 0.03
const N = 200
const dA = 1/(N*N)
const dV = 1/(N*N*N)

mutable struct Particle <: AbstractParticle
    x::RealVector
    testno::Float64
end

function test_area(sys::ParticleSystem, testno, area)
    _area = 0.0
    for p in sys.particles
        if p.testno == testno
            _area += dA
        end
    end
    err = abs(_area/area - 1.0)
    @test err < rtol_2d
    if err > rtol_2d
        println("error in test no ", round(Int,testno))
        println("computed value = ", _area)
        println("correct result = ", area)
    end
end

function test_volume(sys::ParticleSystem, testno, vol)
    _vol = 0.0
    for p in sys.particles
        if p.testno == testno
            _vol += dV
        end
    end
    err = abs(_vol/vol - 1.0)
    @test err < rtol_3d
    if err > rtol_3d
        println("error in test no ", round(Int,testno))
        println("computed value = ", _vol)
        println("correct result = ", vol)
    end
end

function rotmat(x::Float64)::RealMatrix
    return RealMatrix(
         cos(x),  sin(x), 0.,
        -sin(x),  cos(x), 0.,
             0.,      0., 1.
    )
end

function main()
    @testset "area tests" begin
        dom = Rectangle(-10., -10., 10., 10.)
        sys = ParticleSystem(Particle, dom, 1.0)
        grid1 = Grid(1/N, :square)
        grid2 = Grid(1/N, :hexagonal)
        grid3 = Grid(1/N, :vogel)
        out = new_pvd_file("area_tests")

        s1 = Circle(0.0, 0.0, 1.0)
        generate_particles!(sys, grid1, s1, x -> Particle(x, 1.0))
        test_area(sys, 1.0, pi)

        s2 = Rectangle(0., -1., 2.0, 5.0)
        generate_particles!(sys, grid2, s2, x -> Particle(x, 2.0))
        test_area(sys, 2.0, 12.0)

        s3 = Ellipse(0., 0., 4.0, 1.0)
        generate_particles!(sys, grid3, s3, x -> Particle(x, 3.0))
        test_area(sys, 3.0, 4.0*pi)

        tool1 = Rectangle(0.0, -1.0, 4.0, 1.0)
        s4 = s3 - tool1
        generate_particles!(sys, grid1, s4, x -> Particle(x, 4.0))
        test_area(sys, 4.0, 2.0*pi)

        s5 = s3 * tool1
        generate_particles!(sys, grid2, s5, x -> Particle(x, 5.0))
        test_area(sys, 5.0, 2.0*pi)

        s6 = s4 + s5
        generate_particles!(sys, grid3, s6, x -> Particle(x, 6.0))
        test_area(sys, 6.0, 4.0*pi)

        tool2 = Rectangle(-4.0, -1.0, 4.0, 1.0)
        s7 = Specification(tool2, x -> (x[2] < cos(pi*x[1])))
        generate_particles!(sys, grid1, s7, x -> Particle(x, 7.0))
        test_area(sys, 7.0, 8.0)

        s8 = Transform(s2, A=rotmat(pi/7), b=RealVector(-2.0, 0.0, 0.))
        generate_particles!(sys, grid2, s8, x -> Particle(x, 8.0))
        test_area(sys, 8.0, 12.0)

        s9 = Polygon((-1.0, 0.0), (2.0, 0.0), (0.0, 3.0))
        generate_particles!(sys, grid3, s9, x -> Particle(x, 9.0))
        test_area(sys, 9.0, 4.5)

        #save_frame!(out, sys, :testno)
        #save_pvd_file(out)
    end

    @testset "volume tests" begin
        dom = Box(-1., -1., -1, 1., 1., 1.)
        sys = ParticleSystem(Particle, dom, 1.0)
        grid1 = Grid(1/N, :cubic)
        grid2 = Grid(1/N, :facecentered)
        grid3 = Grid(1/N, :bodycentered)
        grid4 = Grid(1/N, :diamond)
        out = new_pvd_file("volume_tests")

        s1 = Box(-0.7, -0.6, -0.5, 0.7, 0.6, 0.5)
        generate_particles!(sys, grid1, s1, x -> Particle(x, 1.0))
        test_volume(sys, 1.0, 1.4*1.2*1.0)

        s2 = Ball(0.0, 0.0, 0.0, 0.8)
        generate_particles!(sys, grid2, s2, x -> Particle(x, 2.0))
        test_volume(sys, 2.0, 4/3*pi*0.8^3)

        s3 = Ellipsoid(0.0, 0.0, 0.0, 0.8, 0.5, 0.3)
        generate_particles!(sys, grid3, s3, x -> Particle(x, 3.0))
        test_volume(sys, 3.0, 4/3*pi*0.8*0.5*0.3)

        s4 = Cone(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.6, 0.3)
        generate_particles!(sys, grid4, s4, x -> Particle(x, 4.0))
        test_volume(sys, 4.0, pi*(2*0.6*0.6 - 0.3*0.3)/3)

        tool1 = Polygon((0.0, 0.0), (0.6, 0.0), (0.0, 0.7))
        s5 = RevolutionBody(tool1)
        generate_particles!(sys, grid1, s5, x -> Particle(x, 5.0))
        test_volume(sys, 5.0, pi/3*0.6*0.6*0.7)

        #save_frame!(out, sys, :testno)
        #save_pvd_file(out)
    end

    return
end


end
