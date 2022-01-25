#=

# 2: Water collapse (explicit)

```@raw html
	<img src='../assets/collapse_exp.png' alt='missing' width="50%" height="50%" /><br>
```

Simulation of a water column collapsing under its own weight onto dry bottom.
This is, where SPH is more useful than typical mesh-based methods
=#

module symplectic

using Printf
include("../src/SPHLib.jl")
using .SPHLib
using Parameters
using Plots
using Ipopt
using JuMP
include("FixPA.jl")

#=
Declare constant parameters
=#

##physical
const dr = 3.0e-3          #average particle distance (decrease to make finer simulation)
const h = 3.0*dr           #size of kernel support
const rho0 = 1000.   	   #fluid density
const m = rho0*dr^2        #particle mass
const g = -9.8*VECY  #gravitational acceleration
const mu = 0.0#8.4e-4          #dynamic viscosity of water

##geometrical
const water_column_width = 0.142
const water_column_height = 0.293
const box_height = 0.35
const box_width = 0.584
const wall_width = 2.5*dr


##artificial
const c = 50.0             #numerical speed of sound
const dr_wall = 0.95*dr
const E_wall = 10*norm(g)*water_column_height
const eps = 1e-16

##temporal
const dt = 0.1*h/c
const t_end = 0.5
const dt_frame = t_end/100

##particle types
const FLUID = 0.
const WALL = 1.

@with_kw mutable struct Particle <: AbstractParticle
	x::RealVector #position
	v::RealVector = VEC0 #velocity
	a::RealVector = VEC0 #acceleration
	P::Float64 = 0. #pressure
	rho::Float64 = 0. #density
    rho0::Float64 = 0.
	type::Float64 #particle_type
end

#=
Define geometry and make particles
=#

function make_system()
	grid = Grid(dr, :square)
	box = Rectangle(0., 0., box_width, box_height)
	fluid = Rectangle(0., 0., water_column_width, water_column_height)
	walls = BoundaryLayer(box, grid, wall_width)
	#walls = Specification(walls, x -> (x[2] < box_height))
	domain = Rectangle(-box_width, -box_width, 2*box_width, 3*box_height)
	sys = ParticleSystem(Particle, domain, h)
	generate_particles!(sys, grid, fluid, x -> Particle(x = x, type = FLUID))
	generate_particles!(sys, grid, walls, x -> Particle(x = x, type = WALL))
	return sys
end

#=
Define particle interactions
=#

@inbounds function find_rho!(p::Particle, q::Particle, r::Float64)
    if p.type == FLUID && q.type == FLUID
		p.rho += m*wendland2(h,r)
	end
end

@inbounds function find_rho0!(p::Particle, q::Particle, r::Float64)
    if p.type == FLUID && q.type == FLUID
		p.rho0 += m*wendland2(h,r)
	end
end

function find_pressure!(p::Particle)
	p.P = c^2*(p.rho - p.rho0)
end

@inbounds function internal_force!(p::Particle, q::Particle, r::Float64)
	if p.type == FLUID && q.type == FLUID
		ker = m*rDwendland2(h,r)
		p.a += -ker*(p.P/rho0^2 + q.P/rho0^2)*(p.x - q.x)
		#p.a += +2*ker*mu/rho0^2*(p.v - q.v)
	elseif p.type == FLUID && q.type == WALL && r < dr_wall
		s = dr_wall/(r + eps)
		p.a += -E_wall/(r + eps)^2*(s^2 - s^4)*(p.x - q.x)
	end	
end

function reset_a!(p::Particle)
    p.a = zero(RealVector)
end

function reset_rho!(p::Particle)
    p.rho = 0.0
end

function move!(p::Particle)
	if p.type == FLUID
		p.x = FixPA.rev_add(p.x, dt*p.v)
	end
end

function accelerate!(p::Particle)
	if p.type == FLUID
		p.v = FixPA.rev_add(p.v, 0.5*dt*(p.a + g))
	end
end

function LJ_potential(p::Particle, q::Particle, r::Float64)::Float64
	if q.type == WALL && p.type == FLUID && r < dr_wall
		s = dr_wall/(r + eps)
		return m*E_wall*(0.5s^2 - 0.25s^4 -0.25)
	else
		return 0.0
	end
end

function energy(sys::ParticleSystem, p::Particle)::Float64
	kinetic = 0.5*m*dot(p.v, p.v)
	internal =  0.5*m*c^2*(p.rho - p.rho0)^2/rho0^2
	gravity_potential = -m*dot(g, p.x)
	wall_potential = SPHLib.sum(sys, LJ_potential, p)
	return kinetic + internal + gravity_potential + wall_potential
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

function save_results!(out::SPHLib.DataStorage, sys::ParticleSystem, k::Int64)
    if (k %  Int64(round(dt_frame/dt)) == 0)
        @printf("t = %.6e\n", k*dt)
        #energy
        E = sum(p -> energy(sys,p), sys.particles)
        @show E
        println("# of part. = ", length(sys.particles))
        println()
        save_frame!(out, sys, :v, :a, :P, :rho, :rho0)
    end
end


function main()
	sys = make_system()
	out = new_pvd_file("results/collapse_fixpa")
    #initialization
    create_cell_list!(sys)
    apply!(sys, find_rho0!, self = true)
    apply!(sys, find_rho!, self = true)
    apply!(sys, find_pressure!)
    apply!(sys, internal_force!)
	for k = 0 : Int64(round(t_end/dt))
        verlet_step!(sys)
        save_results!(out, sys, k)
	end
    #revert velocities
    for p in sys.particles
        p.v = -p.v
    end
    for k = Int64(round(t_end/dt)):-1:0
        verlet_step!(sys)
        save_results!(out, sys, k)
	end
	save_pvd_file(out)
end ## function main

function boltzmann(beta, e)::Float64
	return beta*exp(-e*beta)
end

function plot_energy_distr(path::String)
	N = 100
	domain = Rectangle(-box_width, -box_width, 2*box_width, 3*box_height)
	sys = ParticleSystem(Particle, domain, h)
	read_vtk!(sys, path, x -> Particle(x = x, type = 0.0))
	v_max = sqrt(2*norm(g)*water_column_height)
	dv = v_max/N
	vs = 0.:dv:v_max
	ns = zeros(length(vs))
	for k in 1:length(sys.particles)
		v = norm(sys.particles[k].v)
		n = Int64(round(v/dv))
		if 1 <= n <= length(ns)
			ns[n] += 1.0/(dv*length(sys.particles))
		end
	end
	model = Model(Ipopt.Optimizer)
	
	@variable(model, beta)
	@NLobjective(
            model,
            Min,
            sum((ns[i] - m*beta*vs[i]*exp(-0.5*m*beta*vs[i]^2))^2 for i in 1:length(ns)),
        )
    optimize!(model)
	beta = value(beta)
	ns_boltz = zeros(length(vs))
	for i in 1:length(ns_boltz)
		ns_boltz[i] = m*beta*vs[i]*exp(-0.5*m*beta*vs[i]^2)
	end
	p = plot(vs, [ns ns_boltz], label = ["data" "boltz"])
	savefig(p, "energy_distr.pdf")
	@show(beta)
	T = 1/(beta*1.380649e-23)
	@show(T)
end

end ## module

