#=

# 4: Lid-driven cavity

=#

module cavity_flow

using Printf
include("../src/SPHLib.jl")
using .SPHLib


#=
Declare const parameters (all dims in SI)
=#

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

#=
Declare variables to be stored in a Particle
=#

mutable struct Particle <: AbstractParticle
	x::Vec2	#position
	v::Vec2 #velocity
	a::Vec2 #acceleratation
	rho::Float64 #density
	Drho::Float64 #rate of density
	P::Float64 #pressure
	type::Float64 #particle type
	Particle(x::Vec2, type::Float64) = begin
		v1 = (type == LID ? vlid : 0.)
		return new(x, Vec2(v1, 0.0), zero(Vec2), rho0, 0., 0., type)
	end
end

#=
Define geometry and create particles
=#

function make_system()
	grid = Grid(dr, :hexagonal)
	box = Rectangle(0., 0., llid, llid)
	wall = BoundaryLayer(box, grid, wwall)
	sys = ParticleSystem(Particle, box + wall, h)

	lid   = Specification(wall, x -> x[2] > llid)
	wall = Specification(wall, x -> x[2] <= llid)

	generate_particles!(sys, grid, box, x -> Particle(x, FLUID))
	generate_particles!(sys, grid, lid, x -> Particle(x, LID))
	generate_particles!(sys, grid, wall, x -> Particle(x, WALL))
	return sys
end

#=
Define interactions between particles
=#

#Define interactions between particles

@inbounds function balance_of_mass!(p::Particle, q::Particle, r::Float64)
	if p.type == FLUID
		ker = m*rDwendland2(h,r)
		p.Drho += ker*(dot(p.x-q.x, p.v-q.v))
	end
end

function find_pressure!(p::Particle)
	p.rho += p.Drho*dt
	p.Drho = 0.0
	p.P = c^2*(p.rho-rho0)
end

@inbounds function internal_force!(p::Particle, q::Particle, r::Float64)
	ker = m*rDwendland2(h,r)
	p.a += -ker*(p.P/rho0^2 + q.P/rho0^2)*(p.x - q.x)
	p.a += +2*ker*(mu/rho0^2)*(p.v - q.v)
end

function move!(p::Particle)
	p.a = zero(Vec2)
	if p.type == FLUID
		p.x += dt*p.v
	end
end

function accelerate!(p::Particle)
	if p.type == FLUID
		p.v += 0.5*dt*p.a
	end
end

#=
Time iteration
=#
function main()
	sys = make_system()
	out = new_pvd_file("results/cavity_flow")
	@time for k = 0 : Int64(round(t_end/dt))
		apply!(sys, move!)
		create_cell_list!(sys)
		apply!(sys, balance_of_mass!)
		apply!(sys, find_pressure!)
		apply!(sys, internal_force!)
		apply!(sys, accelerate!)
		if (k % Int64(round(dt_frame/dt)) == 0) #save the frame
			t = k*dt
			@show t
			save_frame!(out, sys, :P, :v, :type)
		end
		apply!(sys, accelerate!)
	end
	save_pvd_file(out)
end

end
