# 4: Lid-driven cavity

````julia
module cavity_flow

using Printf
include("../src/SPHLib.jl")
using .SPHLib
````

Declare const parameters (all dims in SI)

````julia
##geometrical/physical parameters
const Re = 100.0                #Reynolds number
const llid = 0.2                #length of the lid
const mu = 8.4e-4               #viscosity of water
const rho0 = 1.0e+3             #density of water
const vlid = mu*Re/(rho0*llid)     #flow speed of the lid
const dr = llid/70 		        #average particle distance
const h = 2.2*dr		        #size of kernel support
const m = rho0*dr^2             #particle mass
const c = 10*vlid			#numerical speed of sound
const wwall = h

##temporal parameters
const dt = 0.2*h/c              #numerical time-step
const dt_frame = 20.           #how often save data
const t_end = 2000.            #end of simulation

##particle types
const FLUID = 0.
const WALL = 1.
const LID = 2.
````

Declare fields (unknowns)

* `Dx`, `Dy` = velocity
* `DDx`, `DDy` = acceleration
* `rho` = density
* `P` = pressure
* `type` = particle type

````julia
@define_particle Particle Dx Dy DDx DDy rho type P

function main()
````

Define geometry and create particles

````julia
	sys = ParticleSystem(Particle, (-2*wwall, llid + 2*wwall), (-2*wwall, llid + 2*wwall), dr, h)
	box = Rectangle((0., llid), (0., llid))
	walls = BoundaryLayer(box, dr, wwall)
	lid   = Specification(walls, (x,y) -> y > llid)
	walls = Specification(walls, (x,y) -> y <= llid)

	generate_particles!(sys, box, (x,y) -> Particle(x,y; type = FLUID, rho = rho0))
	generate_particles!(sys, lid, (x,y) -> Particle(x,y; type = LID, rho = rho0, Dx = vlid))
	generate_particles!(sys, walls, (x,y) -> Particle(x,y; type = WALL, rho = rho0))
````

Define interactions between particles

````julia
	function balance_of_mass!(p::Particle, q::Particle, r::Float64)
    	if p.type == FLUID
      		p.rho += dt*((p.x - q.x)*(p.Dx - q.Dx) + (p.y - q.y)*(p.Dy - q.Dy))*m*rDwendland2(h,r)
    	end
  	end

	function internal_force!(p::Particle, q::Particle, r::Float64)
		pressure_term = p.P/p.rho^2 + q.P/q.rho^2
		viscous_term = 2.0*mu/(p.rho*q.rho)
		kernel = m*rDwendland2(h,r)
		p.DDx += kernel*(-pressure_term*(p.x - q.x) + viscous_term*(p.Dx - q.Dx))
		p.DDy += kernel*(-pressure_term*(p.y - q.y) + viscous_term*(p.Dy - q.Dy))
	end

	function find_pressure!(p::Particle)
		p.P = c^2*(p.rho - rho0)
	end

	function update!(p::Particle)
		if p.type == FLUID
			p.Dx += p.DDx*dt
			p.Dy += p.DDy*dt
			p.x += p.Dx*dt
			p.y += p.Dy*dt
		end
		p.DDx = 0.0
		p.DDy = 0.0
	end
````

Time iteration

````julia
	out = new_pvd_file("cavity_flow")
	P = ScalarField(sys, :P, "pressure")
	v = VectorField(sys, (:Dx, :Dy), "velocity")
	type = ScalarField(sys, :type, "type")

	println(count(p->(p.type == LID), sys.particles), " lid particles\n")
	println(count(p->(p.type == WALL), sys.particles), " wall particles\n")
	println(count(p->(p.type == FLUID), sys.particles), " fluid particles\n")
	@time for k = 0 : Int64(round(t_end/dt))
		if (k % Int64(round(dt_frame/dt)) == 0) #save the frame
			@printf("t = %.6e\n", k*dt)
			save_frame!(sys, out, P, v, type)
		end
		create_cell_list!(sys)
		apply!(sys, balance_of_mass!)
		apply!(sys, find_pressure!)
		apply!(sys, internal_force!)
		apply!(sys, update!)
	end
	save_pvd_file(out)
end

end
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

