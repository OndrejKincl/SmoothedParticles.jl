#=

# 2: Water collapse (explicit)

```@raw html
	<img src='../assets/collapse_exp.png' alt='missing' width="50%" height="50%" /><br>
```

Simulation of a water column collapsing under its own weight onto dry bottom.
This is, where SPH is more useful than typical mesh-based methods
=#

module collapse_dry

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
const m = rho0*dr^2        #particle mass
const c = 50.0             #numerical speed of sound
const g = Vec2(0.0, -9.8)  #gravitational acceleration
const mu = 8.4e-4          #dynamic viscosity of water
const nu = 1.0e-4           #pressure stabilization

##geometrical
const water_column_width = 0.142
const water_column_height = 0.293
const box_height = 0.35
const box_width = 0.584
const wall_width = 2.5*dr

##temporal
const dt = 0.1*h/c
const t_end = 1.0
const dt_frame = t_end/100


##particle types
const FLUID = 0.
const WALL = 1.

mutable struct Particle <: AbstractParticle
	x::Vec2 #position
	v::Vec2 #velocity
	a::Vec2 #acceleration
	P::Float64 #pressure
	rho::Float64 #density
	Drho::Float64 #rate of density
	type::Float64 #particle_type
	Particle(x, type) = new(
		x, zero(Vec2), zero(Vec2),
		0.,
		rho0, 0.,
		type
	)
end

#=
Define geometry and make particles
=#

function make_system()
	grid = Grid(dr, :square)
	box = Rectangle(0., 0., box_width, box_height)
	fluid = Rectangle(0., 0., water_column_width, water_column_height)
	walls = BoundaryLayer(box, grid, wall_width)
	walls = Specification(walls, x -> (x[2] < box_height))
	domain = Rectangle(-box_width, -box_width, 2*box_width, 3*box_height)
	sys = ParticleSystem(Particle, domain, h)
	generate_particles!(sys, grid, fluid, x -> Particle(x, FLUID))
	generate_particles!(sys, grid, walls, x -> Particle(x, WALL))
	return sys
end

#=
Define particle interactions
=#

@inbounds function balance_of_mass!(p::Particle, q::Particle, r::Float64)
	ker = m*rDwendland2(h,r)
	p.Drho += ker*(dot(p.x-q.x, p.v-q.v) + 2*nu*(p.rho-q.rho))
end

function find_pressure!(p::Particle)
	p.rho += p.Drho*dt
	p.Drho = 0.0
	p.P = c^2*(p.rho - rho0)
end

@inbounds function internal_force!(p::Particle, q::Particle, r::Float64)
	if p.type == FLUID
		ker = m*rDwendland2(h,r)
		p.a += -ker*(p.P/rho0^2 + q.P/rho0^2)*(p.x - q.x)
		p.a += +2*ker*mu/rho0^2*(p.v - q.v)
	end
end

function move!(p::Particle)
	p.a = zero(Vec2)
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
	internal =  0.25*m*c^2*(p.rho - rho0)^2/rho0^2
	return kinetic + potential + internal
end

#=
Put everything into a time loop
=#
function main()
	sys = make_system()
	out = new_pvd_file("results/collapse_dry")
    #a modified Verlet scheme
	for k = 0 : Int64(round(t_end/dt))
        #move particles
        apply!(sys, move!)
        create_cell_list!(sys)
		apply!(sys, balance_of_mass!)
		apply!(sys, find_pressure!)
        apply!(sys, internal_force!)
        apply!(sys, accelerate!)
        #save data at selected frames
        if (k %  Int64(round(dt_frame/dt)) == 0)
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
