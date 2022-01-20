module ICR

using Printf
include("../src/SPHLib.jl")
using .SPHLib
using Parameters
using SparseArrays
using LinearAlgebra
using IterativeSolvers
using IncompleteLU

#=
import Base: eltype, size
import Base.LinearAlgebra: mul!
using Base.Threads

struct SparseMul{Tv,Ti}
    A::SparseMatrixCSC{Tv,Ti}
    B::SparseMatrixCSC{Tv,Ti}
    tmp::Vector{Tv}
    I::Int
    J::Int
    K::Int
    SparseMul{Tv,Ti}(A,B) = begin
        (I, J) = size(A)
        (_J, K) = size(B)
        @assert(J == _J) 
        tmp = zeros(J)
        return new{Tv,Ti}(A,B,tmp, I, J, K)
    end
end

function mul!(y::AbstractVector, M::SparseMul, x::AbstractVector)
    mul!(M.tmp, M.B, x)
    mul!(y, M.A, M.tmp)
end

eltype(M::SparseMul) = eltype(M.A)
size(M::SparseMul) = (M.I, M.K)
=#

#=
Declare constant parameters
=#

const water_column_width = 1.0
const water_column_height = 2.0
const box_height = 3.0
const box_width = 4.0

##physical
const dr = 2.0e-2          #average particle distance (decrease to make finer simulation)
const h = 2.2dr           #size of kernel support
const rho0 = 1000.   	   #fluid density
const m = rho0*dr^2        #particle mass
const vol = dr^2

##temporal
const steps = 20

#const Lmin = 0*kernel(h,0.)*(pi - (dr/h)^2)  #free particles are those that satisfy L < Lmin

#artificial
const eps = 1e-12
const kernel = SPHLib.wendland2
const rDkernel = SPHLib.rDwendland2

@with_kw mutable struct Particle <: AbstractParticle
	x::RealVector #position
	rho::Float64 = 0. #density
    S::RealVector = VEC0
end

#=
Define geometry and make particles
=#

function make_system()
	grid = Grid(dr, :vogel)
    domain = Circle(0., 0., 0.4)
    #domain = Rectangle(-0.3, -0.3, 0.3, 0.3)
    sys = ParticleSystem(Particle, BoundaryLayer(domain, grid, 10h), h)
	generate_particles!(sys, grid, domain, x -> Particle(x = x))
    for p in sys.particles
        #p.x += RealVector(dr*(0.15 - 0.3rand()), dr*(0.15 - 0.3rand()), 0.0)
    end
	return sys
end

#=
Define particle interactions
=#

function find_S!(p::Particle, q::Particle, r::Float64)
	p.S += -2.0*vol*vol*rDkernel(h,r)*(p.x - q.x)
end

@inbounds function find_rho!(p::Particle, q::Particle, r::Float64)
    p.rho += m*kernel(h,r)
end

function reset!(p::Particle)
    p.rho = 0.0
    p.S = VEC0
end

function energy(p::Particle)::Float64
	internal =  0.5*m*(p.rho - rho0)^2/rho0^2
	return internal
end

#=
Functions to build the linear system
=#

function add_element!(I::Vector{Int64}, J::Vector{Int64}, V::Vector{Float64}, i::Int64, j::Int64, v::Float64)
    if v == 0.0
        return
    end
    push!(I, i)
    push!(J, j)
    push!(V, v)
end

function lhs(sys::ParticleSystem)::SparseMatrixCSC{Float64}
	N = length(sys.particles)
    I = Int64[]
    J = Int64[]
    V = Float64[]
    for i in 1:N, j in 1:N
        p = sys.particles[i]
        q = sys.particles[j]
        r = dist(p,q)
        if r > h
            continue
        end
        tmp0 = Float64(p==q)
        tmp1 = m*rDkernel(h,r)*(p.x - q.x)
        tmp2 = 0.5*Float64(p==q)*p.rho^2/m*p.S
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
    return sparse(I, J, V, 3N, 3N)
end

function rhs(sys::ParticleSystem)
    N = length(sys.particles)
    b = zeros(3N)
    for i in 1:N
        p = sys.particles[i]
        b[i+2N] = p.rho - rho0
	end
    return b
end

#=
Put everything into a time loop
=#



function main()
	sys = make_system()
	out = new_pvd_file("results/newton_test")
    P = ParticleField(sys, :P)
    for k in 1:(steps-1)
        create_cell_list!(sys)
        apply!(sys, reset!)
        apply!(sys, find_rho!, self = true)
        apply!(sys, find_S!)
        N = length(sys.particles)
        E = N > 0 ? sum(energy, sys.particles) : 0.0
        @show N,E
        if E < 1e-12
            break
        end
        save_frame!(out, sys, :rho)
        A = lhs(sys)
        b = rhs(sys)
        y = A\b
        for i in 1:N
            sys.particles[i].x += RealVector(y[i], y[i+N], 0.)
        end
	end
    save_frame!(out, sys, :rho, :S)
	save_pvd_file(out)
end ## function main


end ## module