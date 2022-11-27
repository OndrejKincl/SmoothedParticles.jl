#=

# 8: Water collapse (explicit, symplectic, and reversible)

```@raw html
	<img src='../assets/fixpa.png' alt='missing' width="50%" height="50%" /><br>
```

Simulation of a water column collapsing under its own weight onto dry bottom.
Here we use a symplectic scheme and get a reversible simulation. 
At the end of the simulation, the velocities are reverted and the simulation goes back to its initial conditions.
Despite the reversibility, Boltzmann entropy grows and attains its maximum value just before the velocities are reverted.
=#

module two_parts

using Printf
using SmoothedParticles
using Parameters
using Plots
using DataFrames # to store the csv file
using CSV# to store the csv file
include("utils/FixPA.jl")
include("utils/entropy.jl")
include("utils/ICR.jl")
using .FixPA
using .entropy
using LaTeXStrings #for better legends
using Random, Distributions
using LsqFit


#using ReadVTK  #not implemented
#using VTKDataIO

const path = "results/two_parts"
#=
Declare constant parameters
=#

##physical
const dr = 1.0e-2          #average particle distance (decrease to make finer simulation)
const h = 3.0*dr           #size of kernel support
const rho0 = 1000.   	   #fluid density
const m = rho0*dr^2        #particle mass
const g = -9.8*VECY  #gravitational acceleration
const mu = 0.0#8.4e-4          #dynamic viscosity of water

##geometrical
const box_height = 1.0
const box_width = 1.0
const wall_width = 2.5*dr
const slit_height = box_height/10

##artificial
const c = 50.0             #numerical speed of sound
const dr_wall = 0.95*dr
const E_wall = 10*norm(g)
const eps = 1e-6

##temporal
const dt = 0.1*h/c
const t_end = 40.0
const dt_frame = t_end/300

##particle types
const FLUID = 0.
const WALL = 1.
const EMPTY = 2.

mutable struct Particle <: AbstractParticle
	x::RealVector	#position
	v::RealVector #velocity
	a::RealVector #acceleratation
	rho::Float64 #density
	rho0::Float64 #density
	Drho::Float64 #rate of density
	P::Float64 #pressure
	type::Float64 #particle type
	Particle(x::RealVector, type::Float64) = begin
		return new(x, VEC0, VEC0, rho0, rho0, 0., 0., type)
	end
end

#=
Define geometry and make particles
=#

function make_system()
	grid = Grid(dr, :square)
	boxL = Rectangle(0., 0., box_width, box_width)
	boxR = Rectangle(box_width, 0., 2*box_width, box_width)
	wallL = BoundaryLayer(boxL, grid, wall_width)
	wallR = BoundaryLayer(boxR, grid, wall_width)
	sys = ParticleSystem(Particle, boxL + wallL+wallR, h)

	#wallL = Specification(wallL, x -> x[2] <= wall_width)
	#wallR = Specification(wallR, x -> x[2] <= wall_width)

	generate_particles!(sys, grid, boxL, x -> Particle(x, FLUID))
	#generate_particles!(sys, grid, boxR, x -> Particle(x, FLUID))
	generate_particles!(sys, grid, wallL, x -> Particle(x, WALL))
	generate_particles!(sys, grid, wallR, x -> Particle(x, WALL))

	Random.seed!(42) #
	dist = MvNormal(2, 1.0/sqrt(2)) #nondimensional distribution of velocities
	for i in 1:length(sys.particles)
		p = sys.particles[i]
		if p.type == WALL
			if (p.x[1] >= box_width - wall_width) && (p.x[1]<= box_width+wall_width) && (p.x[2] >= box_height/2 - slit_height) && (p.x[2]<= box_height/2 + slit_height)
				p.type = EMPTY
			end
		elseif p.type == FLUID
			vxy = rand(dist, 1)
			p.v = RealVector((vxy[1], vxy[2], 0.0))
		end
	end

	return sys
end

#=
Define particle interactions
=#

@inbounds function internal_force!(p::Particle, q::Particle, r::Float64)
	if p.type == FLUID && q.type == FLUID
		#ker = m*rDwendland2(h,r)
		#p.a += -ker*(p.P/rho0^2 + q.P/rho0^2)*(p.x - q.x)
		#p.a += +2*ker*mu/rho0^2*(p.v - q.v)
	elseif p.type == FLUID && q.type == WALL && r < dr_wall
		s2 = (dr_wall^2 + eps^2)/(r^2 + eps^2)
		p.a += -E_wall/(r^2 + eps^2)*(s2 - s2^2)*(p.x - q.x)
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
		#p.x = rev_add(p.x, dt*p.v)
		p.x = p.x + dt*p.v
	end
end

function accelerate!(p::Particle)
	if p.type == FLUID
		#p.v = rev_add(p.v, 0.5*dt*p.a)
		p.v = p.v + 0.5*dt*p.a
	end
end

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

function LJ_potential(p::Particle, q::Particle, r::Float64)::Float64
	if q.type == WALL && p.type == FLUID && r < dr_wall
		s2 = (dr_wall^2 + eps^2)/(r^2 + eps^2)
		return m*E_wall*(0.25*s2^2 - 0.5*s2 + 0.25)
	else
		return 0.0
	end
end

function energy_kinetic(sys::ParticleSystem)::Float64
	return sum(p -> 0.5*m*dot(p.v, p.v), sys.particles)
end

function left(sys::ParticleSystem)::Float64
	left_number = 0
	for p in sys.particles
		if p.x[1] <= box_width
			left_number += 1
		end
	end
	return left_number
end

function energy(sys::ParticleSystem)
	(E_kin, E_int, E_gra, E_wal, E_tot) = (0., 0., 0., 0., 0.)
	for p in sys.particles
		E_kin += 0.5*m*dot(p.v, p.v)
		#E_int +=  0.5*m*c^2*(p.rho - p.rho0)^2/rho0^2
		#E_gra += -m*dot(g, p.x)
		E_wal += SmoothedParticles.sum(sys, LJ_potential, p)
	end
	#E_tot = E_kin + E_int + E_gra + E_wal
	E_tot = E_kin + E_wal
	return (E_tot, E_kin, E_wal)
end

#=
Put everything into a time loop
=#

function verlet_step!(sys::ParticleSystem)
    apply!(sys, accelerate!)
    apply!(sys, move!)
    create_cell_list!(sys)
    #apply!(sys, reset_rho!)
    #apply!(sys, find_rho!, self = true)
    #apply!(sys, find_pressure!)
    apply!(sys, reset_a!)
    apply!(sys, internal_force!)
    apply!(sys, accelerate!)
end

function save_results!(out::SmoothedParticles.DataStorage, sys::ParticleSystem, k::Int64, E0::Float64)
    if (k %  Int64(round(dt_frame/dt)) == 0)
        @printf("t = %.6e\n", k*dt)
        #energy
        (E_tot, E_kin, E_wal) = energy(sys)
		@show E_tot
		@show E_wal	
		E_err = E_tot - E0
		@show E_err	
        println("# of part. = ", length(sys.particles))
        println()
        save_frame!(out, sys, :v, :a, :type, :P, :rho, :rho0)
    end
end

function main(;revert = false) #if revert=true, velocities are inverted at the end of the simulation and the simulation then goes backward
	sys = make_system()
	out = new_pvd_file(path)
    #initialization
    create_cell_list!(sys)
    #apply!(sys, find_rho0!, self = true)
    apply!(sys, find_rho!, self = true)
    apply!(sys, find_pressure!)
    apply!(sys, internal_force!)

	N_of_particles = length(sys.particles)
	@show(N_of_particles)
	@show(m)
	

	step_final = Int64(round(t_end/dt))
	times = Float64[] #time instants
	#Ss = Float64[] # Entropy values
	#Ekin = Float64[] # Kinetic energy values
	#thermalize!(sys)
	E0 = energy(sys)[1]
	ls = Float64[]
	#Eg = Float64[] # Gravitational energy values
	#Ewall = Float64[] # Wall energy values
	#Eint = Float64[] # Internal energy values
	#Etot = Float64[] # Internal energy values

	for k = 0 : step_final
        verlet_step!(sys)
        save_results!(out, sys, k, E0)
    	if k % round(step_final/100) == 0 # store a number of entropy values
			#distr = velocity_histogram(sys, N = 100)
			#S = entropy_2D_MB(distr)
			#push!(Ss, S)
			#@show(S)

			push!(times, k*dt)
			left_number = left(sys)
			push!(ls, left_number)
			@show left_number

        	#(E_tot, E_kin, E_int, E_gra, E_wal) = energy(sys)
			#push!(Ekin, E_kin)
			#push!(Eg, E_gra)
			#push!(Eint, E_int)
			#push!(Ewall, E_wal)
            #push!(Etot, E_tot)
        	println()
		end
	end

	# Plotting the velocity distribution in comparison with Maxwell-Boltzmann
	#T = plot_velocity_distr(sys, m, path*"/energy_distribution_middle.pdf")

	# Plotting the entropy in time
	#Sred_eq_E = [(1+log(Ekin[k]/(m*length(sys.particles)))) for k in 1:length(Ss)]
	#Sred_eq_T= (1+log(kB*T/m))*ones(Float64, length(Ss))
	#p = plot(times, [Ss Sred_eq_T Sred_eq_E], label = ["entropy" "S_eq(T)" "S_eq(E)"],legend=:bottomright)
	#savefig(p, path*"/entropy_middle.pdf")
	df = DataFrame(time_steps = times, left = ls)
	CSV.write(path*"/results.csv", df)
	p = plot(times, ls, label = "left_number",legend=:bottomright)

	if revert
		#revert velocities
		println("--------------------")
		println("Reverting velocities")
		println("--------------------")
		for p in sys.particles
			p.v = -p.v
		end
	#	Ss_rev = Float64[]
 #       Ekin_rev = Float64[] # Kinetic energy values
 #       Eg_rev = Float64[] # Gravitational energy values
 #       Ewall_rev = Float64[] # Wall energy values
 #       Eint_rev = Float64[] # Internal energy values
 #       Etot_rev = Float64[] # Internal energy values		for k = step_final:-1:0

		for k = step_final:-1:0
			verlet_step!(sys)
        	save_results!(out, sys, k, E0)
			if k % round(step_final/100) == 0 # store a number of entropy values
	#			distr = velocity_histogram(sys, v_max = sqrt(2*norm(g)*water_column_height), N = 100)
	#			S = entropy_2D_MB(distr)
	#			push!(Ss_rev, S)
	#			@show(S)

    #    		(E_tot, E_kin, E_int, E_gra, E_wal) = energy(sys)
				#push!(Ekin_rev, E_kin)
				#push!(Eg_rev, E_gra)
				#push!(Eint_rev, E_int)
				#push!(Ewall_rev, E_wal)
				#push!(Etot_rev, E_tot)

				println()
			end
		end
		#plot_velocity_distr(sys, m, path*"/energy_distribution_final.pdf")

		# Plotting the entropy in time
		#p = plot(times, [Ss Ss_rev Sred_eq_T Sred_eq_E], label = ["entropy forward" "entropy backward" "S_eq(T)" "S_eq(E)"], legend=:bottomright)
		#savefig(p, path*"/entropy_final.pdf")
		#df = DataFrame(time_steps = times, S_Boltzmann = Ss, S_eq_T = Sred_eq_T, S_eq_E = Sred_eq_E)
		#CSV.write(path*"/entropy_final.csv", df)
        # Plotting the energies in time
        #p = plot(times, [Ekin Ekin_rev], label = ["Ekin" "Ekin_rev"], legend=:bottomright)
		#savefig(p, path*"/energy_kin.pdf")
  #      p = plot(times, [Ewall Ewall_rev], label = ["Ewall" "Ewall_rev"], legend=:bottomright)
		#savefig(p, path*"/energy_wall.pdf")
  #      p = plot(times, [Eg Eg_rev ], label = ["Eg" "Eg_rev"], legend=:bottomright)
		#savefig(p, path*"/energy_g.pdf")
  #      p = plot(times, [Eint Eint_rev ], label = ["Eint" "Eint_r"], legend=:bottomright)
		#savefig(p, path*"/energy_int.pdf")
  #      p = plot(times, [Etot Etot_rev], label = ["Etot" "Etot_r"], legend=:bottomright)
		#savefig(p, path*"/energy_tot.pdf")
		#df = DataFrame(time_steps = times, E_kinetic = Ekin, E_wall = Ewall, E_graviational = Eg, E_internal = Eint, E_total = Etot, E_kinetic_rev = Ekin_rev, E_wall_rev = Ewall_rev, E_graviational_rev = Eg_rev, E_internal_rev = Eint_rev, E_total_rev = Etot_rev)
		#CSV.write(path*"/energy.csv", df)	
	end

	save_pvd_file(out)

end ## function main

function plot_left(left_file::String; fit=false)
    df = DataFrame(CSV.File(left_file))
    times = df[:, "time_steps"]
    ls = df[:, "left"]
	if fit
		model(t, p) = p[1] * exp.(-p[2] * t) .+ p[3]
		p0 = [ls[1]/2, 10.0, ls[1]/2]
		fitting = curve_fit(model, times, ls, p0)
		param = fitting.param
		@show param
		curve = param[1] * exp.(-param[2] .* times) .+ param[3]
		plot()
    	p = plot!(times, [ls, curve], legend=:topright, labels=["number of left particles"  "exponential fit"], linewidth=5, thickness_scaling=1)
	else
		plot()
    	p = plot!(times, ls, legend=:topright, label="number of left particles", linewidth=5, thickness_scaling=1)
	end
	savefig(p, path*"/left.pdf")
end


function plot_energy(energy_file::String)
    df = DataFrame(CSV.File(energy_file))
    times = df[:, "time_steps"]
    e_pot = df[:, "E_graviational"]
	Delta_e_pot = e_pot[1]-e_pot[end]
	print("Delta e pot = ", Delta_e_pot)
    e_tot = df[:, "E_total"]
	e_tot0 = e_tot[1]
	e_tot = (e_tot .- e_tot0)./Delta_e_pot
    p = plot(times, e_tot, legend=:topright, label=L"\frac{E_{tot}-E_{tot}(0)}{E_g(end)-E_g(0)}")
	savefig(p, "./energy_tot_scaled.pdf")
end

end ## module

