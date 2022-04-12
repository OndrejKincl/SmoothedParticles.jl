module quadratic
using SmoothedParticles
using Parameters
using Plots
using ForwardDiff
import LinearAlgebra as linalg
using CSV
using DataFrames
import StaticArrays as stars

const ExtendedVector = stars.SVector{5,Float64}
const ExtendedMatrix = stars.SMatrix{5,5,Float64}
const EMAT0 = zero(ExtendedMatrix)
const EVEC0 = zero(ExtendedVector)


#ALGEBRAIC TOOLS

@inbounds function outer(x::RealVector, y::RealVector)::RealMatrix
    return RealMatrix(
        x[1]*y[1], x[2]*y[1], 0., 
        x[1]*y[2], x[2]*y[2], 0.,
        0., 0., 0.
    )
end

@inbounds function det(A::RealMatrix)::Float64
    return A[1,1]*A[2,2] - A[1,2]*A[2,1]
end

@inbounds function inv(A::RealMatrix)::RealMatrix
    idet = 1.0/det(A)
    return RealMatrix(
        +idet*A[2,2], -idet*A[2,1], 0., 
        -idet*A[1,2], +idet*A[1,1], 0.,
        0., 0., 0.
    )
end


@with_kw mutable struct Particle <: AbstractParticle
    x::RealVector
    u::Float64 = 0.
    #partial derivatives
    du::ExtendedVector = EVEC0
    e::Float64 = 0.
    #renormalization matrix
    R::ExtendedMatrix = EMAT0
    vol::Float64
    char::Float64
end

function make_system(h::Float64, dr::Float64)::ParticleSystem
    grid = Grid(dr, :hexagonal)
    #grid.center = RealVector(0.5, 0.5, 0.)
    box = Rectangle(0., 0., 1., 1.)
    wall = BoundaryLayer(box, grid, h)
    dom = box + wall
    sys = ParticleSystem(Particle, dom, h)
    generate_particles!(sys, grid, box, x -> Particle(x=x, u=fun1(x), vol = dr*dr, char = 1.0))
    generate_particles!(sys, grid, wall, x -> Particle(x=x, u=fun1(x), vol = dr*dr, char = 0.0))
    #add noise
    for p in sys.particles
        p.x += 0.1*dr*RealVector(rand() - 0.5, rand() - 0.5, 0.)
        p.u = fun1(p.x)
    end
    create_cell_list!(sys)
    return sys
end

function fun1(x)
    #return sin(pi*x[1])*sin(pi*x[2])
    return x[1]*x[2] + exp(-x[1]^2-x[2]^2)/10
end

function compute_error!(sys::ParticleSystem, fun::Function)::Float64
    dfun = x -> ForwardDiff.gradient(fun, x)
    for p in sys.particles
        p.e = linalg.norm(narrow(p.du) - dfun(p.x))
    end
    E = sqrt(sum(p -> p.vol*p.char*p.e^2, sys.particles))
    @show E
    return E
end

function extend(x::RealVector, h::Float64)::ExtendedVector
    return ExtendedVector(x[1], x[2], x[1]^2/h, x[1]*x[2]/h, x[2]^2/h)
end

function narrow(x::ExtendedVector)::RealVector
    return RealVector(x[1], x[2], 0.)
end

function findgrad!(p::Particle, q::Particle, r::Float64, h::Float64)
    p.du += wendland2(h,r)*(p.u - q.u)*extend(p.x - q.x, h)
    p.R  += wendland2(h,r)*extend(p.x - q.x, h)*extend(p.x - q.x, h)'
end

function renormalize!(p::Particle)
    if p.char == 1
        p.du = Base.inv(p.R)*p.du
    end
end

function test!(out, h::Float64, dr::Float64)::Float64
    sys = make_system(h, dr)
    @show dr
    @show h
    apply!(sys, (p,q,r) -> findgrad!(p,q,r,h))
    apply!(sys, renormalize!)
    E = compute_error!(sys, fun1)
    save_frame!(out, sys, :u, :du, :e, :char)
    println()
    return E
end

function main()
    out = new_pvd_file("quadratic")
    Es = []
    Ns = range(20, 300, 10)
    drs = []
    for N in Ns
        dr = 1.0/N
        E = test!(out, 2.4dr, dr)
        push!(Es, E)
        push!(drs, dr)
    end
    save_pvd_file(out)
    p = plot(
        drs, Es,
        xlabel = "dr", 
        ylabel = "e", 
        legend = :none,
        markershape = :star5
    )
    savefig(p, "convergence_curve_g2.pdf")
    df = DataFrame(dr = drs, E = Es)
	CSV.write("convergence_curve_g2.csv", df)
end

end