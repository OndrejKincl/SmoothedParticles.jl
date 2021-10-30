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

const dr = 3.0e-3     #average particle distance
const h = 1.8*dr      #size of kernel support
const rho0 = 1000.0   #fluid density
const m = rho0*dr^2   #particle mass
const c = 40.0        #numerical speed of sound
const g = 1.0         #gravitational acceleration
const mu = 8.4e-4     #dynamic viscosity of water

const water_depth = 0.14
const box_height = 0.18
const box_width = 0.14
const wall_width = 2*dr


const dt = 0.2*h/c
const dt_frame = 0.1
const t_end = 0.5

##particle types
const FLUID = 0.
const WALL = 1.

#=
Declare fields (unknowns)

* `Dx`, `Dy` = velocity
* `DDx`, `DDy` = acceleration
* `rho` = density
* `type` = particle type
=#

SPHLib.@define_particle Particle Dx Dy DDx DDy rho type

##dependance of pressure on density
@inline function pressure(p::Particle)
	return c^2*(p.rho - rho0)
end

#=
Define geometry and make particles
=#
function main()

	box = Rectangle((0., box_width), (0., box_height))
	fluid = Rectangle((0., box_width), (0., water_depth))
	walls = BoundaryLayer(box, dr, wall_width)

	xrange = (-wall_width, box_width  + wall_width)
	yrange = (-wall_width, box_height + wall_width)
	sys = ParticleSystem(Particle, xrange, yrange, dr, h)
	generate_particles!(sys, fluid, (x,y) -> Particle(x,y; rho = rho0, type = FLUID))
	generate_particles!(sys, walls, (x,y) -> Particle(x,y; rho = rho0, type = WALL))

#=
Define particle interactions
=#

	@inline function balance_of_mass!(p::Particle, q::Particle, r::Float64)
		p.rho += dt*((p.x - q.x)*(p.Dx - q.Dx) + (p.y - q.y)*(p.Dy - q.Dy))*m*rDwendland2(h,r)
	end

	@inline function internal_force!(p::Particle, q::Particle, r::Float64)
		if p.type == FLUID
			pressure_term = pressure(p)/p.rho^2 + pressure(q)/q.rho^2
			viscous_term = 2.0*mu/(p.rho*q.rho)
			kernel = m*rDwendland2(h,r)
			p.DDx += kernel*(-pressure_term*(p.x - q.x) + viscous_term*(p.Dx - q.Dx))
			p.DDy += kernel*(-pressure_term*(p.y - q.y) + viscous_term*(p.Dy - q.Dy))
		end
	end

	@inline function update!(p::Particle)
		if p.type == FLUID
			p.Dx += p.DDx*dt
			p.Dy += (p.DDy - g)*dt
			p.x += p.Dx*dt
			p.y += p.Dy*dt
		end
		p.DDx = 0.0
		p.DDy = 0.0
	end


#=
Put everything into a time loop
=#

	type = ScalarField(sys, :type, "type")
	rho = ScalarField(sys, :rho, "density")
	v = VectorField(sys, (:Dx, :Dy), "velocity")

	out = new_pvd_file("static_container")
	for k = 0 : Int64(round(t_end/dt))
		if (k %  Int64(round(dt_frame/dt)) == 0)
			@printf("t = %.6e\n", k*dt)
			save_frame!(sys, out, rho, v, type)
		end
		create_cell_list!(sys)
		apply!(sys, balance_of_mass!)
		apply!(sys, internal_force!)
		apply!(sys, update!)
	end
	save_pvd_file(out)
end ##function main()

end ##module
