#=

# Symplecticity & reversibility

```@raw html
	<img src='../assets/fixpa.png' alt='missing' width="50%" height="50%" /><br>
```

Simulation of a water column collapsing under its own weight onto dry bottom.
Here we use a symplectic scheme and get a reversible simulation. 
At the end of the simulation, the velocities are reverted and the simulation goes back to its initial conditions.
Despite the reversibility, Boltzmann entropy grows and attains its maximum value just before the velocities are reverted.
=#

module collapse_symplectic_viscous

using Printf
using SmoothedParticles
using Parameters
using Plots
using DataFrames # to store the csv file
using CSV# to store the csv file
include("utils/FixPA.jl")
include("utils/entropy.jl")
using .FixPA
using .entropy


#using ReadVTK  #not implemented
#using VTKDataIO

#=
Declare constant parameters
=#

##physical
const dr = 1.0e-2          # average particle distance (decrease to refine, increase to speed up)
const h = 3.0*dr           # kernel radius
const rho0 = 1000.   	   # fluid density
const m = rho0*dr^2        # particle mass
const g = -9.8*VECY        # gravitational acceleration
const mu = 8.4e-4      # dynamic viscosity of water
const cv = 4184.0          # specific heat capacity for water
const T0 = 293.15          # initial temperature

##geometrical
const water_column_width = 1.0
const water_column_height = 2.0
const box_height = 3.0
const box_width = 4.0
const wall_width = 2.5*dr


##artificial
const c = 50.0             #numerical speed of sound
const dr_wall = 0.95*dr
const E_wall = 10*norm(g)*water_column_height
const eps = 1e-16

##temporal
const dt = 0.1*h/c
const t_end = 1.0
const dt_frame = t_end/100

## output
const OUT_DIR = "results/collapse_fixpa_viscous"

##particle types
const FLUID = 0.
const WALL = 1.

@with_kw mutable struct Particle <: AbstractParticle
	x::RealVector #position
	v::RealVector = VEC0 #velocity
	a::RealVector = VEC0 #acceleration
	P::Float64 = 0. #pressure
	S::Float64 = 0. #entropy
	T::Float64 = T0 #temperature
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
    # Barotropic part of the pressure
	p.P = c^2*(p.rho - p.rho0)
end

function find_temperature!(p::Particle)
	p.T = T0 * exp(p.S / (m * cv))
end

@inbounds function internal_force!(p::Particle, q::Particle, r::Float64)
	if p.type == FLUID && q.type == FLUID
        x_pq = p.x - q.x
		ker = m*rDwendland2(h,r)
		p.a += -ker*(p.P/p.rho^2 + q.P/q.rho^2)*(p.x - q.x)
		#p.a += +2*ker*mu/rho0^2*(p.v - q.v)
    	p.a += 8.0*ker*mu/(p.rho*q.rho)*dot(p.v-q.v, x_pq)/(r*r + 0.01*h*h)*x_pq # viscous
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
		p.v = rev_add(p.v, 0.5*dt*(p.a + g))
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

function energy(sys::ParticleSystem)
	fluid_particles = Iterators.filter(p -> p.type == FLUID, sys.particles)

	E_kin = sum(p -> 0.5 * m * dot(p.v, p.v), fluid_particles)
    # The internal energy is the sum of "cold" compressive energy and thermal energy.
    # This formula for E_comp is the integral of P/ρ^2 dρ, making it consistent with the pressure P = c^2(ρ - ρ₀).
	E_comp = sum(p -> m * c^2 * (p.rho0/p.rho + log(p.rho/p.rho0) - 1.0), p for p in fluid_particles if p.rho > 0 && p.rho0 > 0)
    # Thermal energy is defined relative to the initial temperature T0.
	E_therm = sum(p -> m * cv * (p.T - T0), fluid_particles)
	E_int = E_comp + E_therm
	E_gra = sum(p -> -m * dot(g, p.x), fluid_particles)
	E_wal = 0.5 * sum(p -> SmoothedParticles.sum(sys, LJ_potential, p), fluid_particles)
	
	E_tot = E_kin + E_int + E_gra + E_wal
	return (E_tot, E_kin, E_int, E_gra, E_wal)
end

function entropy_production!(p::Particle, q::Particle, r::Float64)
	if p.type == FLUID && q.type == FLUID
		ker = rDwendland2(h,r)
		x_pq = p.x - q.x
		v_pq = p.v - q.v
    	p.S += - 4.0*m*m*ker*mu/(p.T*p.rho*q.rho)*dot(v_pq, x_pq)^2/(r*r + 0.01*h*h)*dt #viscous
	end
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
	apply!(sys, entropy_production!)
    apply!(sys, find_temperature!)
    apply!(sys, find_pressure!)
    apply!(sys, reset_a!)
    apply!(sys, internal_force!)
    apply!(sys, accelerate!)
end

function save_results!(out::SmoothedParticles.DataStorage, sys::ParticleSystem, k::Int64)
    if (k %  Int64(round(dt_frame/dt)) == 0)
        @printf("t = %.6e\n", k*dt)
        (E_tot, E_kin, E_int, E_gra, E_wal) = energy(sys)
		@show E_tot
		@show E_kin
		@show E_int
		@show E_gra
		@show E_wal
        println("# of part. = ", length(sys.particles))
        println()
        save_frame!(out, sys, :v, :a, :P, :rho, :rho0, :S, :T)
    end
end

function main(;revert = true) #if revert=true, velocities are inverted at the end of the simulation and the simulation then goes backward
	sys = make_system()
    if !isdir(OUT_DIR)
        mkdir(OUT_DIR)
    end
	out = new_pvd_file(OUT_DIR)

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
	Etot_arr = Float64[]
	Ekin_arr = Float64[]
	Eint_arr = Float64[]
	Egra_arr = Float64[]
	Ewal_arr = Float64[]
	Ss_arr = Float64[] # Entropy values

	for k = 0 : step_final
        verlet_step!(sys)
        save_results!(out, sys, k)
    	if k % round(step_final/100) == 0 # store a number of entropy values
			distr = velocity_histogram(sys, N = 100)
			S = sum(p->p.S, sys.particles)
			(Etot, Ekin, Eint, Egra, Ewal) = energy(sys)

			push!(times, k*dt)
			push!(Ss_arr, S)
			push!(Etot_arr, Etot)
			push!(Ekin_arr, Ekin)
			push!(Eint_arr, Eint)
			push!(Egra_arr, Egra)
			push!(Ewal_arr, Ewal)

			@show(S)
			@show(Etot, Ekin, Eint, Egra, Ewal)
        	println()
		end
	end

	gr() # Set the backend for Plots.jl to ensure savefig works
	p1 = plot(times, Ss_arr, label = "entropy", legend=:bottomright)
	savefig(p1, joinpath(OUT_DIR, "entropy.pdf"))
	p2 = plot(times, [Etot_arr, Ekin_arr, Eint_arr, Egra_arr, Ewal_arr], label=["Total" "Kinetic" "Internal" "Gravity" "Wall"], legend=:left)
	savefig(p2, joinpath(OUT_DIR, "energy.pdf"))
	df = DataFrame(time=times, E_total=Etot_arr, E_kinetic=Ekin_arr, E_internal=Eint_arr, E_gravity=Egra_arr, E_wall=Ewal_arr, Entropy=Ss_arr)
	CSV.write(joinpath(OUT_DIR, "diagnostics.csv"), df)

	save_pvd_file(out)

end ## function main

end ## module
