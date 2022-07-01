module test_collision_2d
# two colliding circles
# test conservation of energy

include("../src/SmoothedParticles.jl")
using .SmoothedParticles
using Parameters
using Plots
using Test
using Statistics

#=
Declare constant parameters
=#

##physical
const dr = 2.0e-2          #average particle distance (decrease to make finer simulation)
const h = 2.4*dr           #size of kernel support
const rho0 = 1000.   	   #fluid density
const m = rho0*dr^2        #particle mass
const c = 20.0
const v0 = 1.0*VECX

##geometrical
const circ_rad = 0.4
const dom_len = 20.0
const dom_wid = 20.0
const deltaX = 1.0
const deltaY = 0.2

##temporal
const dt = 0.1*h/c
const t_end = 1.0
const dt_frame = t_end/10

@with_kw mutable struct Particle <: AbstractParticle
	x::RealVector #position
	v::RealVector = VEC0 #velocity
	a::RealVector = VEC0 #acceleration
	P::Float64 = 0. #pressure
	rho::Float64 = 0. #density
    rho0::Float64 = 0.
end

#=
Define geometry and make particles
=#

function make_system()
	grid = Grid(dr, :square)
    circ1 = Circle(-0.5*deltaX, -0.5*deltaY, circ_rad)
    circ2 = Circle( 0.5*deltaX,  0.5*deltaY, circ_rad)
	domain = Rectangle(-0.5*dom_len, -0.5*dom_wid, 0.5*dom_len, 0.5*dom_wid)
	sys = ParticleSystem(Particle, domain, h)
	generate_particles!(sys, grid, circ1, x -> Particle(x = x, v =  v0))
	generate_particles!(sys, grid, circ2, x -> Particle(x = x, v = -v0))
	return sys
end

#=
Define particle interactions
=#

@inbounds function find_rho!(p::Particle, q::Particle, r::Float64)
	p.rho += m*wendland2(h,r)
end

@inbounds function find_rho0!(p::Particle, q::Particle, r::Float64)
	p.rho0 += m*wendland2(h,r)
end

function find_pressure!(p::Particle)
	p.P = c^2*(p.rho - p.rho0)
end

@inbounds function internal_force!(p::Particle, q::Particle, r::Float64)
	ker = m*rDwendland2(h,r)
	p.a += -ker*(p.P/rho0^2 + q.P/rho0^2)*(p.x - q.x)
end

function reset_a!(p::Particle)
    p.a = VEC0
end

function reset_rho!(p::Particle)
    p.rho = 0.0
end

function move!(p::Particle)
	p.x += dt*p.v
end

function accelerate!(p::Particle)
	p.v += 0.5*dt*p.a
end

function energy(p::Particle)::Float64
	kinetic = 0.5*m*dot(p.v, p.v)
	internal =  0.5*m*c^2*(p.rho - p.rho0)^2/rho0^2
	return kinetic + internal
end

#=
Put everything into a time loop
=#

function verlet_step!(sys::ParticleSystem)
    apply!(sys, accelerate!)
    apply!(sys, move!)
    create_cell_list!(sys)
    apply!(sys, reset_rho!)
    apply!(sys, find_rho!, self = true)
    apply!(sys, find_pressure!)
    apply!(sys, reset_a!)
    apply!(sys, internal_force!)
    apply!(sys, accelerate!)
end

function main()
    #out = new_pvd_file("test_collision_2d")
    println("running test_collision_2d.jl")
	sys = make_system()
    #initialization
    create_cell_list!(sys)
    apply!(sys, find_rho0!, self = true)
    apply!(sys, find_rho!, self = true)
    apply!(sys, find_pressure!)
    apply!(sys, internal_force!)
    t = Float64[]
    N = Float64[]
    E = Float64[]
	for k = 0 : Int64(round(t_end/dt))
        verlet_step!(sys)
        if (k % Int64(round(dt_frame/dt)) == 0)
            t = k*dt
            @show t
            push!(N, length(sys.particles))
            push!(E, sum(p -> energy(p), sys.particles))
            #save_frame!(out, sys, :v)
        end
	end
    @testset "count particles" begin
        @test all(N .== N[1])
    end
    @testset "energy conservation" begin
        err = maximum(E/E[1] .- 1.0)
        @test err < 1e-2
    end
    #save_pvd_file(out)
end ## function main

end ## module