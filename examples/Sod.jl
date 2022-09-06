#=

# 9: Sod schock test

```@raw html
	<img src='../assets/Sod.png' width="50%" height="50%" alt='missing' /><br>
```
 
```@raw html
A 2D simulation of the Sod schock test.
 <a href="https://en.wikipedia.org/wiki/Sod_shock_tube">Wiki</a>
```

=#

module Sod

using Printf
using SmoothedParticles

const folder_name = "results/Sod"

#=
Declare constants
=#

#physical parameters
const rho0 = 1.0		#referential fluid density
const rhoL = 1.0		#left fluid density
const rhoR = 0.125		#right fluid density
const pL = 1.0		#left pressure
const pR = 0.1		#right pressure
const gamma = 5/3
csL = sqrt(gamma * pL / rhoL) #speed of sound left
csR = sqrt(gamma * pR / rhoR) #speed of sound right
@show csL
@show csR
const c = max(csL, csR)	#numerical speed of sound
#const mu = 1.0e-3		#dynamic viscosity of water
#const nu = 1.0e-3		#pressure stabilization
#const P0 = 1.2          #anti-clump term

#geometry parameters
const chan_l = 1.0           #length of the channel
const chan_w = chan_l/10     #width of the channel
const dr = chan_w / 20 		     #average particle distance for the reference density
@show dr
const m = rho0*dr^2		#particle mass
const drL = dr * sqrt(rho0 / rhoL) 		     #average particle distance (decrease to make finer simulation)
const drR = dr * sqrt(rho0 / rhoR)		     #average particle distance (decrease to make finer simulation)
@show drL
@show drR
const h = 2.5*max(drL, drR)		     #size of kernel support
const wall_w = 2.5*dr        #width of the wall
const LRboundary = chan_l/10 #boundary position between L and R

#temporal parameters
const dt = 0.2*h/c      #time step
@show dt
const t_end = 1      #end of simulation
const dt_frame = t_end/10    #how often data is saved
@show dt_frame

#particle types
const FLUID = 0.
const WALL = 1.

#=
Declare variables to be stored in a Particle
=#

mutable struct Particle <: AbstractParticle
    x::RealVector #position
    v::RealVector #velocity
    a::RealVector #acceleration
    rho::Float64 #density
    Drho::Float64 #rate of density
    P::Float64 #pressure
    type::Float64 #particle type
    Particle(x,type) = begin
        return new(x, VEC0, VEC0, rho0, 0., 0., type)
    end
end

function make_system()
    grid = Grid(dr, :square)
	box = Rectangle(0., 0., chan_l, chan_w)
	boxL = Rectangle(0., 0., LRboundary, chan_w)
	boxR = Rectangle(LRboundary, 0., chan_l, chan_w)
	wall = BoundaryLayer(box, grid, wall_w)
	sys = ParticleSystem(Particle, box + wall, h)
    gridL = Grid(drL, :square)
    gridR = Grid(drR, :square)

	#wall = Specification(wall, x -> x[2] <= llid)

	generate_particles!(sys, grid, wall, x -> Particle(x, WALL))
	generate_particles!(sys, gridL, boxL, x -> Particle(x, FLUID))
	generate_particles!(sys, gridR, boxR, x -> Particle(x, FLUID))

	#push!(sys.particles, Particle(x=particle_position, v=particle_velocity, type = FLUID))	

    return sys
end

#Define interactions between particles

@inbounds function balance_of_mass!(p::Particle, q::Particle, r::Float64)
	ker = m*rDwendland2(h,r)
	#p.Drho += ker*(dot(p.x-q.x, p.v-q.v) + 2*nu*(p.rho-q.rho))
	p.Drho += ker*(dot(p.x-q.x, p.v-q.v))
end

function find_pressure!(p::Particle)
	p.rho += p.Drho*dt
	p.Drho = 0.0
	p.P = rho0*c^2*((p.rho/rho0)^7 - 1.0)/7
end

@inbounds function internal_force!(p::Particle, q::Particle, r::Float64)
	ker = m*rDwendland2(h,r)
	p.a += -ker*(p.P/rho0^2 + q.P/rho0^2)*(p.x - q.x)
	#p.a += +2*ker*mu/rho0^2*(p.v - q.v)
    ker = m*rDwendland2(h/2,r)
    #p.a += -2*ker*P0/rho0^2*(p.x - q.x)
end

function move!(p::Particle)
	p.a = VEC0
	if p.type == FLUID 
		p.x += dt*p.v
	end
end

function accelerate!(p::Particle)
	if p.type == FLUID
		p.v += 0.5*dt*p.a
	end
end

function  main()
    sys = make_system()
	out = new_pvd_file(folder_name)

    #a modified Verlet scheme
	for k = 0 : Int64(round(t_end/dt))
        t = k*dt
        apply!(sys, move!)

        create_cell_list!(sys)
		apply!(sys, balance_of_mass!)
        apply!(sys, find_pressure!)
        apply!(sys, internal_force!)
        apply!(sys, accelerate!)
        #save data at selected frames
        if (k %  Int64(round(dt_frame/dt)) == 0)
            @show t
            save_frame!(out, sys, :v, :P, :type)
        end
		apply!(sys, accelerate!)
	end
	save_pvd_file(out)
end

end
