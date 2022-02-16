#=

# 9: Kepler vortex (explicit, symplectic, and reversible)

```@raw html
	<img src='../assets/Kepler.png' alt='missing' width="50%" height="50%" /><br>
```
Kepler vortex ... TODO
=#

module Kepler_vortex

using Printf
using SmoothedParticles
using Parameters
using Plots
using DataFrames # to store the csv file
using CSV # to store the csv file
using QuadGK #To calculate the integral over Σ
using Interpolations
using Roots
include("utils/FixPA.jl")
include("utils/entropy.jl")
using .FixPA
using .entropy

#=
Declare constant parameters
=#
##graviational
const r0 = 10.0 #central ring radius
const GM = 1000.

function vphi_r(r::Float64)::Float64
	return sqrt(GM)/sqrt(r)
end
const vphi0 = vphi_r(r0)
@show vphi0
const omega0= vphi0/r0
@show omega0

function Σ(r::Float64)::Float64
	return 2*pi*r*exp(-30*(1-r/r0)^2)
end
const denominator = quadgk(Σ, 0, 40, rtol=1.0e-06)[1]
@show denominator
function f(r::Float64)::Float64
	return quadgk(Σ, 0, r, rtol=1.0e-03)[1]/denominator
end

const rs = 0.:0.5:25.0
const f_table = [f(r) for r in rs]
const f_interpolation = interpolate(f_table, BSpline(Cubic(Line(OnGrid()))))
const f_interpolation_scaled = scale(f_interpolation, rs)
function r_f(F::Float64)::Float64
	return find_zero(x -> (f_interpolation_scaled(x)-F), r0)
end
@show r_f(0.5)

## computational
#const rs_in_vortex = 7.:dr:14. # equidistant rings in the vortex. Boundaries chosen as the edges of the histogram for 10^4 points over the interval [0,1]
const N_rings = 25
@show N_rings
const rs_in_vortex = [r_f(u) for u in 0.01:(0.99-0.01)/N_rings:0.99] # Gaussian rings
dr = r_f(0.25+1/N_rings)-r_f(0.25)
@show dr

const h = 3.0*dr           #size of kernel support
const rho0 = 1.   	   #fluid density
const m = rho0*dr^2        #particle mass
#const mu = 0.0#8.4e-4          #dynamic viscosity of water

##geometrical
const box_width = 4*r0
const wall_width = 2.5*dr

##artificial
const c = 0.01             #numerical speed of sound
const dr_wall = 0.95*dr
const E_wall = GM/r0 
const eps = 1e-16

##temporal
const dt = 0.0001*h/c
@show dt
const t_end = 100 * 2 * pi / omega0 #ten revolutions
@show t_end
const dt_frame = t_end/200

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

function generate_circle!(sys::ParticleSystem, r::Float64, dφ::Float64; φ_start=0.0, x_center=0.0, y_center=0.0, vφ= 0.0)
	φ = φ_start
	while φ < 2*pi+φ_start
		x = cos(φ) #x position on unit circle
		y = sin(φ) #y position on unit circle
		particle_position = RealVector(x_center+r*x, y_center+r*y, 0.)
		particle_velocity = RealVector(-vφ*y, vφ*x, 0.)
		push!(sys.particles, Particle(x=particle_position, v=particle_velocity, type = FLUID))	
		φ += dφ
	end
end



function make_system()
    domain = Rectangle(-box_width, -box_width, box_width, box_width)
    sys = ParticleSystem(Particle, domain, h)

	dφ = rs_in_vortex[2]/rs_in_vortex[1]-1.0
	n_particles = 0
	for i in 1:length(rs_in_vortex)-1 # not the last one because dphi would be unknown
		r = rs_in_vortex[i]
		generate_circle!(sys, r, dφ; vφ = vphi_r(r))	
		dφ = (rs_in_vortex[i+1]-r)/r
		n_particles += Int64(round(2*pi/dφ))
	end
	@show n_particles

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
		p.x = rev_add(p.x, dt*p.v)
	end
end

function accelerate!(p::Particle)
	if p.type == FLUID
		p.v = rev_add(p.v, 0.5*dt*(rev_add(p.a, -GM/(norm(p.x)^3)*p.x)))
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

function energy_potential(sys::ParticleSystem)::Float64
	return sum(p -> -GM*m/(norm(p.x)), sys.particles)
end

function energy_internal(sys::ParticleSystem)::Float64
	return sum(p -> 0.5*m*c^2*(p.rho - p.rho0)^2/rho0^2, sys.particles)
end


function energy(sys::ParticleSystem, p::Particle)::Float64
	kinetic = 0.5*m*dot(p.v, p.v)
	internal =  0.5*m*c^2*(p.rho - p.rho0)^2/rho0^2
	gravity_potential = - GM*m/(norm(p.x))
	wall_potential = SmoothedParticles.sum(sys, LJ_potential, p)
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

function save_results!(out::SmoothedParticles.DataStorage, sys::ParticleSystem, k::Int64)
    if (k %  Int64(round(dt_frame/dt)) == 0)
        @printf("t = %.6e\n", k*dt)
        #energy
        Etot = sum(p -> energy(sys,p), sys.particles)
        @show Etot
        Ekin = energy_kinetic(sys)
        @show Ekin
        Epot = energy_potential(sys)
        @show Epot
        Eint = energy_internal(sys)
        @show Eint
        println("# of part. = ", length(sys.particles))
        println()
        save_frame!(out, sys, :v, :a, :P, :rho, :rho0)
    end
end

function main(;revert = true) #if revert=true, velocities are inverted at the end of the simulation and the simulation then goes backward
	sys = make_system()
	out = new_pvd_file("results/Kepler_vortex")
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
#	times = Float64[] #time instants
#	Ss = Float64[] # Entropy values
#	Ekin = Float64[] # Kinetic energy values
	for k = 0 : step_final
        verlet_step!(sys)
        save_results!(out, sys, k)
#    	if k % round(step_final/100) == 0 # store a number of entropy values
#			distr = velocity_histogram(sys, N = 100)
#			S = entropy_2D_MB(distr)
#			push!(times, k*dt)
#			push!(Ss, S)
#			push!(Ekin, energy_kinetic(sys))
#			@show(S)
#        	println()
#		end
	end

	# Plotting the velocity distribution in comparison with Maxwell-Boltzmann
	#T = plot_velocity_distr(sys, m, "results/Kepler_vortex/energy_distribution_middle.pdf")

	# Plotting the entropy in time
#	sred_eq_E = [(1+log(Ekin[k]/(m*length(sys.particles)))) for k in 1:length(Ss)]
#	sred_eq_T= (1+log(kB*T/m))*ones(Float64, length(Ss))
#	p = plot(times, [Ss Sred_eq_T Sred_eq_E], label = ["entropy" "S_eq(T)" "S_eq(E)"],legend=:bottomright)
#	savefig(p, "results/Kepler_vortex/entropy_middle.pdf")
#	df = DataFrame(time_steps = times, S_Boltzmann = Ss, S_eq_T = Sred_eq_T, S_eq_E = Sred_eq_E)
#	CSV.write("results/Kepler_vortex/entropy_middle.csv", df)

#	if revert
#		#revert velocities
#		println("--------------------")
#		println("Reverting velocities")
#		println("--------------------")
#		for p in sys.particles
#			p.v = -p.v
#		end
#		Ss_rev = Float64[]
#		for k = step_final:-1:0
#			verlet_step!(sys)
#			save_results!(out, sys, k)
#			if k % round(step_final/100) == 0 # store a number of entropy values
#				distr = velocity_histogram(sys, v_max = sqrt(2*norm(g)*water_column_height), N = 100)
#				S = entropy_2D_MB(distr)
#				push!(Ss_rev, S)
#				@show(S)
#				println()
#			end
#		end
#		plot_velocity_distr(sys, m, "results/Kepler_vortex/energy_distribution_final.pdf")
#
#		# Plotting the entropy in time
#		p = plot(times, [Ss Ss_rev Sred_eq_T Sred_eq_E], label = ["entropy forward" "entropy backward" "S_eq(T)" "S_eq(E)"], legend=:bottomright)
#		savefig(p, "results/Kepler_vortex/entropy_final.pdf")
#		df = DataFrame(time_steps = times, S_Boltzmann = Ss, S_eq_T = Sred_eq_T, S_eq_E = Sred_eq_E)
#		CSV.write("results/Kepler_vortex/entropy_final.csv", df)
#	end

	save_pvd_file(out)

end ## function main

end ## module

