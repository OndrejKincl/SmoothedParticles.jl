#=

# 1: Static container

```@raw html
	<img src='../assets/static_container.png' width="50%" height="50%" alt='missing' /><br>
```

This is hello world of SPH. 
Simulates motionless fluid in uniform gravitational field.
Things that may be tested by this benchmark are: 
* walls,
* zero point stability of time integrator,
* stability of free surface.

Let us begin by declaring a module and importing some stuff.
=#

module static_container

using Printf
include("../src/SPHLib.jl")
using .SPHLib

#=
Declare constant parameters
=#

const dr = 3.0e-3       #average particle distance
const h = 1.8*dr        #size of kernel support
const rho0 = 1000.0     #fluid density
const m = rho0*dr^2     #particle mass
const c = 40.0          #numerical speed of sound
const g = Vec2(0., -1.) #gravitational acceleration
const mu = 8.4e-4       #dynamic viscosity of water

const water_depth = 0.14
const box_height = 0.18
const box_width = 0.14
const wall_width = 2*dr

##temporal parameters
const dt = 0.2*h/c
const dt_frame = 0.1
const t_end = 0.5

##particle types
const FLUID = 0.
const WALL = 1.

mutable struct Particle <: AbstractParticle
	x::Vec2
	v::Vec2
	a::Vec2
	rho::Float64
	type::Float64
	Particle(x::Vec2, type::Float64) = new(
		x,
		zero(Vec2),
		zero(Vec2),
		rho0,
		type
	)
end


##dependance of pressure on density
@inline function pressure(p::Particle)
	return c^2*(p.rho - rho0)
end

##fluid identier
function isfluid(p::Particle)::Float64
	return Float64(p.type == FLUID)
end

#=
Define geometry and make particles
=#
function make_system()
	box = Rectangle(0., 0., box_width, box_height)
	fluid = Rectangle(0., 0., box_width, water_depth)
	walls = BoundaryLayer(box, dr, wall_width)
	domain = box + walls
	sys = ParticleSystem(Particle, domain, h)
	generate_particles!(sys, fluid, x -> Particle(x, FLUID), dr)
	generate_particles!(sys, walls, x -> Particle(x,  WALL), dr)
	return sys
end

#=
Define particle interactions
=#
@inline function balance_of_mass!(p::Particle, q::Particle, r::Float64)
	p.rho += dt*dot(p.x-q.x, p.v-q.v)*m*rDwendland2(h,r)
end

@inline function internal_force!(p::Particle, q::Particle, r::Float64)
	if p.type == FLUID
		ker = m*rDwendland2(h,r)
		#pressure force
		p.a += -ker*(pressure(p)/p.rho^2 + pressure(q)/q.rho^2)*(p.x - q.x)
		#viscous force
		p.a += ker*2*mu/(p.rho*q.rho)*(p.v - q.v)
	end
end

@inline function update!(p::Particle)
	p.v += dt*isfluid(p)*(p.a + g)
	p.x += dt*isfluid(p)*p.v
	p.a = zero(Vec2)
end


#=
Put everything into a time loop
=#
function main()
	sys = make_system()
	@show sys.key_max
	out = new_pvd_file("static_container")
	for k = 0 : Int64(round(t_end/dt))
		if (k %  Int64(round(dt_frame/dt)) == 0)
			println("# of particles = ", length(sys.particles))
			@printf("t = %.6e\n", k*dt)
			save_frame!(sys, out, :rho, :type, :v)
		end
		create_cell_list!(sys)
		apply!(sys, balance_of_mass!)
		apply!(sys, internal_force!)
		apply!(sys, update!)
	end
	save_pvd_file(out)
end ##function main()

end ##module