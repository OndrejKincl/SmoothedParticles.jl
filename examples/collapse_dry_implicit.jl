#=

# 3: Water collapse (implicit)

```@raw html
	<img src='../assets/collapse_exp.png' width="50%" height="50%" alt='missing' /><br>
```

Simulation of a water column collapsing under its own weight onto dry bottom.
Here with strictly incompressible approach (Projection method).
=#

module collapse_dry_implicit

using Printf
include("../src/SPHLib.jl")
using .SPHLib
using LinearAlgebra
using IterativeSolvers
using IncompleteLU

#=
Declare constant parameters
=#

##kernel functions
const kernel = spline23
const Dkernel = Dspline23
const rDkernel = rDspline23

const dr = 4.0e-3       #average particle distance (decrease to make finer simulation)
const h = 2.8*dr        #size of kernel support
const rho0 = 1000.0     #fluid density
const g = 9.8           #gravitational acceleration
const mu = 8.4e-4       #dynamic viscosity
const m = dr^2*rho0     #particle mass

const Lmin = 3*kernel(h,0.)/rho0*(pi - (dr/h)^2)  #free particles are those that satisfy L < Lmin

##geometry parameters
const water_column_width = 0.142
const water_column_height = 0.293
const box_height = 0.35
const box_width = 0.584
const nlayers = 3 #number of wall layers
const wall_width = nlayers*dr

##temporal parameters
const dt = h/20.0
const dt_frame = 0.01
const t_end = 0.5

##labels for particle types
const FLUID = 0.
const  WALL = 1.
const DUMMY = 2.

#=
Declare fields (unknowns)

* `Dx`, `Dy` = velocity
* `DDx`, `DDy` = acceleration
* `type` = particle type
* `P` = pressure
* `div` = divergence of velocity
* `L` = value determining whether particle lies on a free surface
=#

SPHLib.@define_particle Particle Dx Dy DDx DDy type P div L

#=
Define geometry and create particles
=#
function main()

	container = Rectangle((0., box_width), (0., box_height))
	water_column = Rectangle((0., water_column_width), (0., water_column_height))
	wall = BoundaryLayer(container, dr, dr)
	wall = Specification(wall, (x,y) -> (y < box_height))
	dummy_wall = BoundaryLayer(container, dr, wall_width)
	dummy_wall = Specification(dummy_wall, (x,y) -> (y < box_height))
	dummy_wall = BooleanDifference(dummy_wall, wall)

	xrange = (-wall_width, box_width + wall_width)
	yrange = (-wall_width, 3*box_height)
	sys = ParticleSystem(Particle, xrange, yrange, dr, h)
	generate_particles!(sys, water_column, (x,y) -> Particle(x, y; type = FLUID))
	generate_particles!(sys, 		 wall, (x,y) -> Particle(x, y; type =  WALL))
	generate_particles!(sys,   dummy_wall, (x,y) -> Particle(x, y; type = DUMMY))

#=
Particle interactions
=#

	function initialize!(p::Particle)
		if p.type == FLUID
			p.x += dt*p.Dx
			p.y += dt*p.Dy
			##gravity
			p.Dy += -g*dt
		end
		p.div = 0.
		p.L = 0.
	end

	function viscous_force!(p::Particle, q::Particle, r::Float64)
		temp = 2.0*m*mu*rDkernel(h,r)/rho0^2
		p.DDx += temp*(p.Dx - q.Dx)
		p.DDy += temp*(p.Dy - q.Dy)
	end

	function find_div_and_L!(p::Particle, q::Particle, r::Float64)
		p.div += -((p.x - q.x)*(p.Dx - q.Dx) + (p.y - q.y)*(p.Dy - q.Dy))*m*rDkernel(h,r)/rho0
		p.L += -2.0*m*rDkernel(h,r)/rho0^2
	end

	@fastmath function internal_force!(p::Particle, q::Particle, r::Float64)
		temp =  (p.P + q.P)/rho0^2
		temp *= m*rDkernel(h,r)
		p.DDx += -temp*(p.x - q.x)
		p.DDy += -temp*(p.y - q.y)
	end

	function accelerate!(p::Particle)
		if p.type == FLUID
			p.Dx += dt*p.DDx
			p.Dy += dt*p.DDy
		end
		p.DDx = 0.
		p.DDy = 0.
	end

#=
Functions to build the linear system
=#

	@fastmath function minus_laplace(p::Particle, q::Particle, r::Float64)::Float64
		if p == q
			return p.type == DUMMY ? p.L : max(p.L, Lmin)
		end
		return 2.0*m*rDkernel(h,r)/rho0^2
	end

	@fastmath function rhs(p::Particle)::Float64
		return -p.div/dt
	end


#=
Time iteration
=#

	v = VectorField(sys, (:Dx, :Dy), "velocity")
	P = ScalarField(sys, :P, "pressure")
	type = ScalarField(sys, :type, "type")
	out_pvd = new_pvd_file("collapse_dry_implicit")
	out_txt = open("collapse_dry_implicit/data.txt", "w")
	N0 = length(sys.particles)

	for k = 0 : Int64(round(t_end/dt))
		if (k %  Int64(round(dt_frame/dt)) == 0)
			@printf("t = %.6e\n", k*dt)
			save_frame!(sys, out_pvd, v, P, type)
			dimless_time = string(k*dt*sqrt(g/water_column_height))
			leading_edge = maximum(p -> (p.type == FLUID ? p.x - water_column_width : 0.), sys.particles)/water_column_height
			write(out_txt, string(dimless_time)*" "*string(leading_edge)*"\n")
		end
		apply!(sys, initialize!)
		create_cell_list!(sys)
		apply!(sys, viscous_force!)

		##assemble linear system and solve for pressure
		apply!(sys, find_div_and_L!)
		A = assemble_matrix(sys, minus_laplace)
		b = assemble_vector(sys, rhs)
		try
			P .= cg(A, b; Pl = ilu(A; τ = 5.0))
		catch
			save_pvd_file(out_pvd)
			close(out_txt)
			error("Unable to solve linear system.")
		end
		apply!(sys, internal_force!)
		apply!(sys, accelerate!)
	end
	save_pvd_file(out_pvd)
	close(out_txt)
end

end
