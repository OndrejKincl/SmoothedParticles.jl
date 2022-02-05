#=

# 7: Water collapse 3d

=#

module collapse3d

using Printf
include("../src/SPHLib.jl")
using .SPHLib

#=
Declare constant parameters
=#

##physical
const dr = 5.0e-3          #average particle distance (decrease to make finer simulation)
const h = 2.0*dr           #size of kernel support
const rho0 = 1000.   	   #fluid density
const m = rho0*dr^3        #particle mass
const c = 50.0             #numerical speed of sound
const g = -9.8*VECZ        #gravitational acceleration
const mu = 8.4e-4          #dynamic viscosity of water
const nu = 1.0e-4          #pressure stabilization

##geometrical
const water_column_width = 0.142
const water_column_height = 0.293
const box_height = 0.35
const box_width = 0.584
const box_depth = 0.15
const wall_width = 2.5*dr

##temporal
const dt = 0.1*h/c
const t_end = 0.5
const dt_frame = t_end/200


##particle types
const FLUID = 0.
const WALL = 1.

#=
Declare variables to be stored in a Particle
=#

mutable struct Particle <: AbstractParticle
	x::RealVector #position
	v::RealVector #velocity
	a::RealVector #acceleration
	P::Float64 #pressure
	rho::Float64 #density
	Drho::Float64 #rate of density
	type::Float64 #particle_type
	Particle(x, type) = new(
		x, VEC0, VEC0,
		0.,
		rho0, 0.,
		type
	)
end

#=
Define geometry and make particles
=#

function make_system()
	grid = Grid(dr, :cubic)
	box = Box(0., 0., 0., box_width, box_height, box_depth)
	fluid = Box(0., 0., 0., water_column_width, water_column_height, box_depth)
	walls = BoundaryLayer(box, grid, wall_width)
	walls = Specification(walls, x -> (x[2] < box_height))
	domain = SPHLib.boundarybox(walls)
	sys = ParticleSystem(Particle, domain, h)
	generate_particles!(sys, grid, fluid, x -> Particle(x, FLUID))
	generate_particles!(sys, grid, walls, x -> Particle(x, WALL))
	return sys
end

#=
Define particle interactions
=#

@inbounds function balance_of_mass!(p::Particle, q::Particle, r::Float64)
	ker = m*rDwendland3(h,r)
	p.Drho += ker*(dot(p.x-q.x, p.v-q.v) + 2*nu*(p.rho-q.rho))
end

function find_pressure!(p::Particle)
	p.rho += p.Drho*dt
	p.Drho = 0.0
	p.P = c^2*(p.rho - rho0)
end

@inbounds function internal_force!(p::Particle, q::Particle, r::Float64)
	if p.type == FLUID
		ker = m*rDwendland3(h,r)
		p.a += -ker*(p.P/rho0^2 + q.P/rho0^2)*(p.x - q.x)
		p.a += +2*ker*mu/rho0^2*(p.v - q.v)
	end
end

function move!(p::Particle)
	p.a = VEC0
	if p.type == FLUID
		p.x += dt*p.v
	end
end

function accelerate!(p::Particle)
	if p.type == FLUID
		p.v += 0.5*dt*(p.a + g)
	end
end

function energy(p::Particle)::Float64
	kinetic = 0.5*m*dot(p.v, p.v)
	potential = -m*dot(g, p.x)
	internal =  0.5*m*c^2*(p.rho - rho0)^2/rho0^2
	return kinetic + potential + internal
end

#=
Put everything into a time loop
=#
function main()
	sys = make_system()
	out = new_pvd_file("results/collapse3d")
    println("# of parts = ", length(sys.particles))
	#a modified Verlet scheme
	@time for k = 0 : Int64(round(t_end/dt))
	#move particles
		apply!(sys, move!)
		create_cell_list!(sys)
		apply!(sys, balance_of_mass!)
		apply!(sys, find_pressure!)
		apply!(sys, internal_force!)
		apply!(sys, accelerate!)
		#save data at selected frames
		if (k % Int64(round(dt_frame/dt)) == 0)
			@printf("t = %.6e\n", k*dt)
			@printf("E = %.6e\n", sum(energy, sys.particles))
			@printf("\n")
			save_frame!(out, sys, :v, :P, :rho, :type)
		end
		#accelerate
		apply!(sys, accelerate!)
	end
	save_pvd_file(out)
end ## function main

end ## module
