# 2: Water collapse (explicit)

```@raw html
	<img src='../assets/collapse_exp.png' alt='missing' width="50%" height="50%" /><br>
```

Simulation of a water column collapsing under its own weight onto dry bottom.
This is, where SPH is more useful than typical mesh-based methods

````julia
module collapse_dry

using Printf
include("../src/SPHLib.jl")
using .SPHLib
````

Declare constant parameters

````julia
##physical
const dr = 5.0e-3       #average particle distance (decrease to make finer simulation)
const h = 2.0*dr        #size of kernel support
const rho0 = 1.0     #fluid density
const m = rho0*dr^2     #particle mass
const c = 50.0          #numerical speed of sound
const g = 1.0           #gravitational acceleration
const mu = 0.0         #dynamic viscosity of water
const nu = 1.0         #pressure stabilization

##geometrical
const water_column_width = 0.142
const water_column_height = 0.293
const box_height = 0.35
const box_width = 0.584
const wall_width = 2.0*dr

##temporal
const dt = 0.0001 #0.2*h/c
const dt_frame = dt
const t_end = 50*dt


##particle types
const FLUID = 0.
const WALL = 1.
````

Declare fields (unknowns)

* `Dx`, `Dy` = velocity
* `Fx`, `Fy` = force
* `rho` = density
* `P` = pressure
* `Drho` = rate of density
* `type` = particle type

````julia
SPHLib.@define_particle Particle Dx Dy Fx Fy rho Drho P type
````

Define geometry and make particles

````julia
function main()
	container = Rectangle((0., box_width), (0., box_height))
	water_column = Rectangle((0., water_column_width), (0., water_column_height))
	walls = BoundaryLayer(container, dr, wall_width)
	walls = Specification(walls, (x,y) -> (y < box_height))

	xrange = (-100*wall_width, box_width + 100*wall_width)
	yrange = (-100*wall_width, 3*box_height)
	sys = ParticleSystem(Particle, xrange, yrange, dr, h)
	generate_particles!(sys, water_column, (x,y) -> Particle(x, y; type = FLUID, rho = rho0))
	generate_particles!(sys, 		walls, (x,y) -> Particle(x, y; type =  WALL, rho = rho0))
````

Define particle interactions

````julia
	function balance_of_mass!(p::Particle, q::Particle, r::Float64)
        p.Drho += ( (p.x - q.x)*(p.Dx - q.Dx) + (p.y - q.y)*(p.Dy - q.Dy) #divergence of velocity
````

````julia
				  )*m*rDwendland2(h,r)
	end

	function find_pressure!(p::Particle)
		p.rho += p.Drho*dt
		p.Drho = 0.0
		p.P = c^2*(p.rho - rho0)
	end

	function pressure_force!(p::Particle, q::Particle, r::Float64)
		temp = -(p.P/rho0^2 + q.P/rho0^2)*m*m*rDwendland2(h,r)
		p.Fx += temp*(p.x - q.x)
		p.Fy += temp*(p.y - q.y)
	end

   function viscous_force!(p::Particle, q::Particle, r::Float64)
        temp = 2.0*mu/(rho0^2)*m*m*rDwendland2(h,r)
        p.Fx += temp*(p.Dx - q.Dx)
		p.Fy += temp*(p.Dy - q.Dy)
	end

	function move!(p::Particle)
		p.Fx = 0.0
		p.Fy = 0.0
        #if p.type != WALL
            p.x += p.Dx*dt
            p.y += p.Dy*dt
        #end
	end

	function gravity!(p::Particle)
		if p.type != WALL
			p.Fy -= g
		end
	end

	function accelerate!(p::Particle)
        #if p.type == FLUID
            p.Dx += 0.5*p.Fx*dt/m
            p.Dy += 0.5*p.Fy*dt/m
        #end
    end

	function energy(p::Particle)::Float64
		return m*(0.5*p.Dx*p.Dx + 0.5*p.Dy*p.Dy + g*p.y + 0.5*c^2*(p.rho - rho0)^2/(2*rho0^2))
	end
````

Put everything into a time loop

````julia
	type = ScalarField(sys, :type, "type")
	P = ScalarField(sys, :P, "pressure")
	rho = ScalarField(sys, :rho, "rho")
	v = VectorField(sys, (:Dx, :Dy), "velocity")
	out = new_pvd_file("collapse_dry")

    #a modified Verlet scheme
	for k = 0 : Int64(round(t_end/dt))
        #move particles
        apply!(sys, move!)
        create_cell_list!(sys)
        #compute forces
		apply!(sys, gravity!)
        #apply!(sys, viscous_force!)
		#apply!(sys, balance_of_mass!)
		#apply!(sys, find_pressure!)
		#apply!(sys, pressure_force!)
        #accelerate
        apply!(sys, accelerate!)
        #save data at selected frames
        if (k %  Int64(round(dt_frame/dt)) == 0)
            @printf("t = %.6e\n", k*dt)
			@printf("E = %.6e\n", sum(energy, sys.particles))
			@printf("\n")
            save_frame!(sys, out, v, P, rho, type)
        end
        #accelerate
		apply!(sys, accelerate!)
	end
	save_pvd_file(out)
end ## function main

end ## module
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

