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
#using ReadVTK 
#using VTKDataIO

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
const t_end = 1.00
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

function energy_kinetic(sys::ParticleSystem)::Float64
	return sum(p -> 0.5*m*dot(p.v, p.v), sys.particles)
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

function main(;revert = true)
	sys = make_system()
	out = new_pvd_file("results/collapse_fixpa")
    #initialization
    create_cell_list!(sys)
    apply!(sys, find_rho0!, self = true)
    apply!(sys, find_rho!, self = true)
    apply!(sys, find_pressure!)
    apply!(sys, internal_force!)

	N_of_particles = length(sys.particles)
	@show(N_of_particles)
	@show(m)
	

	step_final = Int64(round(t_end/dt))
	times = Float64[] #time instants
	Ss = Float64[] # Entropy values
	Ekin = Float64[] # Kinetic energy values
	for k = 0 : step_final
        verlet_step!(sys)
        save_results!(out, sys, k)
    	if k % round(step_final/100) == 0 # store a number of entropy values
			distr = velocity_histogram(sys, N = 100)
			S = entropy_2D_MB(distr)
			push!(times, k*dt)
			push!(Ss, S)
			push!(Ekin, energy_kinetic(sys))
			@show(S)
		end
	end

	# Plotting the velocity distribution in comparison with Maxwell-Boltzmann
	T = plot_velocity_distr(sys, "energy_distribution_middle.pdf")

	# Plotting the entropy in time
	Sred_eq_E = [(1+log(Ekin[k]/(m*length(sys.particles)))) for k in 1:length(Ss)]
	Sred_eq_T= (1+log(kB*T/m))*ones(Float64, length(Ss))
	p = plot(times, [Ss Sred_eq_T Sred_eq_E], label = ["entropy" "S_eq(T)" "S_eq(E)"],legend=:bottomright)
	savefig(p, "entropy_middle.pdf")

	if revert
		#revert velocities
		println("--------------------")
		println("Reverting velocities")
		println("--------------------")
		for p in sys.particles
			p.v = -p.v
		end
		Ss_rev = Float64[]
		for k = step_final:-1:0
			verlet_step!(sys)
			save_results!(out, sys, k)
			if k % round(step_final/100) == 0 # store a number of entropy values
				distr = velocity_histogram(sys, v_max = sqrt(2*norm(g)*water_column_height), N = 100)
				S = entropy_2D_MB(distr)
				push!(Ss_rev, S)
				@show(S)
			end
		end
		plot_velocity_distr(sys, "energy_distribution_final.pdf")

		# Plotting the entropy in time
		p = plot(times, [Ss Ss_rev Sred_eq_T Sred_eq_E], label = ["entropy forward" "entropy backward" "S_eq(T)" "S_eq(E)"], legend=:bottomright)
		savefig(p, "entropy_final.pdf")
	end

	save_pvd_file(out)

end ## function main

#function boltzmann(beta, e)::Float64
#	return beta*exp(-e*beta)
#end

"""
	Histogram(xs::Vector{Float64}, ys::Vector{Float64},N,dx)

Histogram structure storing ``N`` x values ``xs`` with uniform bin width ``dx`` and ``N`` y values ``ys``.
"""
struct Histogram
	xs::Vector{Float64}
	ys::Vector{Float64}
	N::Int64
	dx::Float64
end

"""
	velocity_histogram(sys::ParticleSystem; v_max = 0, N = 10)

Building the histrogram of 2D velocities (norms) with ``v_max`` the maximum velocity in the histogram and ``N`` bins.
"""
function velocity_histogram(sys::ParticleSystem; v_max = 0.0, N = 100)::Histogram
	if v_max == 0.0 # if v_max = 0, find the maximum velocity of the particles
		for k in 1:length(sys.particles) 
			v = norm(sys.particles[k].v)
			if v > v_max
				v_max = v
			end
		end
	end

	# Find the heights of the histogram bins
	dv = v_max/N # velocity increment between the bins
	vs = 0.:dv:v_max
	ns = zeros(length(vs))
	for k in 1:length(sys.particles)
		v = norm(sys.particles[k].v)
		n = Int64(round(v/dv))
		if 1 <= n <= length(ns)
			ns[n] += 1.0/(dv*length(sys.particles))
		end
	end

	return Histogram(vs,ns,100,dv)
end



"""
	kB

Boltzmann constant (in the SI units)
"""
const kB = 1.380649e-23

"""
	entropy(fMB::Histogram)::Float64

Calculate Boltzmann entropy of a 2D Maxwell-Boltzmann distribution approximated by an ``fMB`` histogram.
"""
function entropy_2D_MB(fMB::Histogram)::Float64
	@assert(fMB.xs[1] == 0) # Assuming that the histogram starts at zero velocity

	S = 0.0

	# Approximating the reduced entropy near v=0, where a numerical singularity could appear
	fMBder = (fMB.ys[2]-fMB.ys[1])/fMB.dx
	if fMBder > 0
		S = - fMB.ys[1] * (log(fMBder)*fMB.dx - fMBder*(fMB.dx^3)/6)
	end

	# Approximating the rest of entropy
	for k in 2:length(fMB.xs)
		if fMB.xs[k] != 0
			if fMB.ys[k] > 0
				S += -fMB.ys[k] * log(fMB.ys[k]/fMB.xs[k]) * fMB.dx
			end
		end
	end

	return S
end

""" 
	plot_velocity_distr(sys::ParticleSystem, name::String)

Plots the distribution of velocity magnitudes among the particles of ``sys`` and saves the resulting pdf,
together with a fit of the Maxwell-Boltzmann distribution, to file ``name``.
Returns the temperature.
"""
function plot_velocity_distr(sys::ParticleSystem, name::String; v_max = 0.0)::Float64
	distr = velocity_histogram(sys, v_max=v_max, N = 100)
	
	# fitting the histogram to a 2D Maxwell-Boltzmann distribution
	model = Model(Ipopt.Optimizer)
	@variable(model, beta)
	@NLobjective(
			model,
			Min,
			sum((distr.ys[i] - m*beta*distr.xs[i]*exp(-0.5*m*beta*distr.xs[i]^2))^2 for i in 1:length(distr.xs)),
		) 
	optimize!(model)
	beta = value(beta)

	# Plotting both the actual histogram and the fitted Maxwell-Boltzmann distribution
	ns_boltz = zeros(length(distr.xs))
	for i in 1:length(ns_boltz)
		ns_boltz[i] = m*beta*distr.xs[i]*exp(-0.5*m*beta*distr.xs[i]^2)
	end
	@show(beta)
	T = 1/(beta*kB)
	@show(T)

	p = plot(distr.xs, [distr.ys ns_boltz], label = ["data" "Maxwell-Boltzmann, T="*string(T)])
	savefig(p, name)
	return T
end

"""
	plot_velocity_distr(path::String, name:: String)

Plot velocity distribution of particles loaded from a VTK file.
"""
function plot_velocity_distr(path::String, name::String)
	@error "Not working, sorry"
	domain = Rectangle(-box_width, -box_width, 2*box_width, 3*box_height) 
	sys = ParticleSystem(Particle, domain, h)
	read_vtk!(sys, path, x -> Particle(x = x, type = 0.0))
	plot_velocity_distr(sys, name)  #not yet working
end

end ## module

