#=

# 8: Vertical heat convection

```@raw html
	<img src='../assets/fixpa.png' alt='missing' width="50%" height="50%" /><br>
```
Simulation of heat convection
=#

module heat 

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

#=
Declare constant parameters
=#

##physical
const dr = 2.0e-2          #average particle distance (decrease to make finer simulation)
const h = 3.0*dr           #size of kernel support
const g = -9.8*VECY*10.0  #gravitational acceleration
const mu = 8.4e-2          #dynamic viscosity of water
const gamma = 1.6

const folder_name = "heat"
#const cv = 4184.0
const cv = 1.0
const p0 = 1.0e01
const rho0 = 10.0
const c0 = 10.0*sqrt(p0*gamma/rho0)
const m = rho0*dr^2        #particle mass
const kB = 1.380649E-23
const T0 = c0^2/(gamma*(gamma-1.0)*cv)
@show T0

const m0 = rho0*dr*dr
const S0 = m0
@show rho0
@show m0
@show S0
@show c0
const Tdown = 700.0
const Tup= 680.0
@show Tdown
@show Tup
const lambda = 1.0e-2 #heat exchange coefficient at the boundary
const bc_width = h
const lambda_F = 1.0e-2 #heat conductivity (times temperature squared)

##geometrical
const box_height = 0.5
const box_width = 2.0
const wall_width = 2.5*dr

##artificial
const dr_wall = 0.95*dr
const E_wall = 2*norm(g)
const eps = 1e-6

##temporal
const dt = 0.01*h/c0
@show dt
const t_end = 5.0
@show t_end
const dt_frame = t_end/500
@show dt_frame

##particle types
const FLUID = 0.
const WALL = 1.
const EMPTY = 2.


mutable struct Particle <: AbstractParticle
	x::RealVector
    m::Float64
    S::Float64 
    v::RealVector 
    a::RealVector 
    rho::Float64 
    rho0::Float64 
    s::Float64 
    P::Float64 
    T::Float64 
    q::RealVector #heat
	type::Float64 #particle type
	Particle(x::RealVector, type::Float64) = begin
		return new(x, m0, S0, VEC0, VEC0, 0.0, 0.0, 0.0, 0.0, 0.0, VEC0, type)
	end
end

#=
Define geometry and make particles
=#

function make_system()
	grid = Grid(dr, :square)
	box = Rectangle(0., 0., box_width, box_height)
	wall = BoundaryLayer(box, grid, wall_width)
	sys = ParticleSystem(Particle, box + wall, h)

	#wallL = Specification(wallL, x -> x[2] <= wall_width)
	#wallR = Specification(wallR, x -> x[2] <= wall_width)

	generate_particles!(sys, grid, box, x -> Particle(x, FLUID))
	#generate_particles!(sys, grid, boxR, x -> Particle(x, FLUID))
	generate_particles!(sys, grid, wall, x -> Particle(x, WALL))

	return sys
end

#=
Define particle interactions
=#

@inbounds function internal_force!(p::Particle, q::Particle, r::Float64)
	if p.type == FLUID && q.type == FLUID
		ker = q.m*rDwendland2(h,r)
        x_pq = p.x - q.x
		p.a += -ker*(p.P/p.rho^2 + q.P/q.rho^2)*(p.x - q.x)
    	#p.a += -ker*(p.T/p.rho*(q.S/q.m - p.s/p.rho))*(p.x - q.x)
    	#p.a += -ker*(q.T/q.rho*(p.S/p.m - q.s/q.rho))*(p.x - q.x)
		#p.a += +2*ker*mu/rho0^2*(p.v - q.v)
    	p.a += 8.0*ker*mu/(p.rho*q.rho)*dot(p.v-q.v, x_pq)/(r*r + 0.01*h*h)*x_pq
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
		p.a = VEC0
		p.x += dt*p.v
		#reset rho, s and a
		p.rho = 0.
		p.s = 0.
        p.q = VEC0
		p.a = VEC0
	end
end

function accelerate!(p::Particle)
	if p.type == FLUID
		#p.v = rev_add(p.v, 0.5*dt*p.a)
		p.v = p.v + 0.5*dt*(p.a + g)
	end
end

@inbounds function find_rho!(p::Particle, q::Particle, r::Float64)
	if p.type == FLUID && q.type == FLUID
        p.rho += q.m*wendland2(h,r)
        #p.s   += q.S*wendland2(h,r)
    end
end

@inbounds function apply_rho0!(p::Particle)
	if p.type == FLUID
		p.rho += p.rho0
	end
end

@inbounds function find_s!(p::Particle)
	if p.type == FLUID
        p.s   = p.S*p.rho/p.m
    end
end

@inbounds function find_rho0!(p::Particle, q::Particle, r::Float64)
    if p.type == FLUID && q.type == FLUID
		p.rho0 += m*wendland2(h,r)
	end
end

@inbounds function set_rho0!(p::Particle)
    if p.type == FLUID
		p.rho0 = rho0*c0^2/(gamma*(gamma-1.0)*cv*T0)
	end
end

function eint(rho::Float64, s::Float64)::Float64
	return rho*c0^2/(gamma*(gamma-1.0))*(rho/rho0)^(gamma-1.0)*exp(s/(cv*rho))+rho0*c0^2/gamma - p0
end

function find_P!(p::Particle)
	if p.type == FLUID
    	p.T = c0^2/(gamma*(gamma-1.0))*(p.rho/rho0)^(gamma-1.0)*exp(p.s/(p.rho*cv))/cv
		#p.P = p.rho*c0^2 *(p.rho/rho0)^(gamma-1.0) * exp(p.s/(cv*p.rho))/gamma - rho0*c0^2/gamma + p0 - (p.rho/rho0)^(gamma-1.0)*exp(p.s/(cv*p.rho))*p.s/(cv*gamma*(gamma-1.0)) + p.rho*p.S*p.T/p.m
		p.P = (gamma-1.0)*eint(p.rho, p.s)- (rho0*c0^2-gamma*p0)
	end
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

function energy(sys::ParticleSystem)
	(E_kin, E_int, E_gra, E_wal, E_tot) = (0., 0., 0., 0., 0.)
	for p in sys.particles
		if p.type == FLUID
			E_kin += 0.5*m*dot(p.v, p.v)
			E_int += eint(p.rho, p.s)/p.rho*m
			#E_int +=  0.5*m*c^2*(p.rho - p.rho0)^2/rho0^2
			E_gra += -m*dot(g, p.x)
			E_wal += SmoothedParticles.sum(sys, LJ_potential, p)
		end
	end
	#E_tot = E_kin + E_int + E_gra + E_wal
	E_tot = E_kin +E_int + E_wal + E_gra
	return (E_tot, E_kin, E_int, E_gra, E_wal)
end

function bc!(sys::ParticleSystem)
	for p in sys.particles
		if p.type == FLUID && p.x[2] < bc_width # bottom
			p.S += p.m * lambda * cv * (Tdown-p.T) * dt
		end
		#if p.type == FLUID && p.x[2] > box_height - bc_width #top
		#	p.S += p.m * lambda * cv * (Tup -p.T) * dt
		#end
	end
end

#=
Put everything into a time loop
=#

function find_heat!(p::Particle, q::Particle, r::Float64)
	if p.type == FLUID && q.type == FLUID
		ker = rDwendland2(h,r)
		x_pq = p.x - q.x
        p.q += lambda_F * q.m/q.rho*(p.T-q.T)*ker*x_pq 
    end
end

function entropy_production!(p::Particle, q::Particle, r::Float64)
	if p.type == FLUID && q.type == FLUID
		ker = rDwendland2(h,r)
		x_pq = p.x - q.x
		v_pq = p.v - q.v
    	p.S += - 4.0*p.m*q.m*ker*mu/(p.T*p.rho*q.rho)*dot(v_pq, x_pq)^2/(r*r + 0.01*h*h)*dt #viscous
        p.S += p.m*q.m/(p.rho*p.T*q.rho)*dot(p.q+q.q,x_pq)*ker* dt #Fourier
	end
end

function verlet_step!(sys::ParticleSystem)
    apply!(sys, accelerate!)
    apply!(sys, move!)
    create_cell_list!(sys)
    #apply!(sys, reset_rho!)
    #apply!(sys, find_rho!, self = true)
    #apply!(sys, find_pressure!)
    apply!(sys, reset_a!)
    apply!(sys, find_rho!, self=true)
	apply!(sys, apply_rho0!)
    apply!(sys, find_s!)
    apply!(sys, find_P!)
	apply!(sys, find_heat!)
	apply!(sys, entropy_production!)
    apply!(sys, internal_force!)
    apply!(sys, accelerate!)
end

function save_results!(out::SmoothedParticles.DataStorage, sys::ParticleSystem, k::Int64, E0::Float64)
    if (k %  Int64(round(dt_frame/dt)) == 0)
         save_frame!(out, sys, :v, :a, :type, :P, :s, :T, :rho)
    end
end

function ensure_uint64_vtk_header!(filepath::String)
	old_header = codeunits("<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">")
	new_header = codeunits("<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\" compressor=\"vtkZLibDataCompressor\">")
	header_type_marker = codeunits("header_type=\"UInt64\"")
	data = read(filepath)
	if findfirst(header_type_marker, data) !== nothing
		return
	end
	match = findfirst(old_header, data)
	if match === nothing
		return
	end
	patched = Vector{UInt8}(undef, length(data) - length(old_header) + length(new_header))
	patched[1:first(match)-1] = data[1:first(match)-1]
	patched[first(match):first(match)+length(new_header)-1] = new_header
	patched[first(match)+length(new_header):end] = data[last(match)+1:end]
	write(filepath, patched)
end

function fix_vtp_headers!(folder::String)
	for file in readdir(folder; join=true)
		endswith(file, ".vtp") || continue
		ensure_uint64_vtk_header!(file)
	end
end

function main(;heating = true) #if heating=true, the bottom edge is heated to Tdown and the upper edge cooled to Tup
	sys = make_system()
	out = new_pvd_file(folder_name)
    #initialization
    create_cell_list!(sys)
    #apply!(sys, find_rho0!, self = true)
    #for p in sys.particles
    #    p.S = S0*p.x[2]
    #end
	apply!(sys, set_rho0!)
    apply!(sys, find_rho!, self = true)
	apply!(sys, apply_rho0!)
    apply!(sys, find_s!)
    apply!(sys, find_P!)
    apply!(sys, find_heat!)
    apply!(sys, internal_force!)

	N_of_particles = length(sys.particles)
	@show(N_of_particles)
	@show(m)

	step_final = Int64(round(t_end/dt))
	times = Float64[] #time instants
	#thermalize!(sys)
	E0 = energy(sys)[1]
	initial_T = average_T(sys)
	@show initial_T
	Ts = Float64[] # Entropy values
	Ekin = Float64[] # Kinetic energy values
	Ewall = Float64[] # Wall energy values
	Eint = Float64[] # Internal energy values
	Eg = Float64[] # gravitational energy values
	Etot = Float64[] # Internal energy values

	for k = 0 : step_final
        verlet_step!(sys)
		if heating
			bc!(sys)
		end
        save_results!(out, sys, k, E0)
    	if k % round(step_final/100) == 0 # store a number of entropy values
       		@printf("t = %.6e\n", k*dt)
			#distr = velocity_histogram(sys, N = 100)
			#S = entropy_2D_MB(distr)
			#push!(Ss, S)
			#@show(S)

			push!(times, k*dt)

			#energy
			(E_tot, E_kin, E_int, E_g, E_wal) = energy(sys)
			@show E_tot
			push!(Etot, E_tot)
			@show E_kin
			push!(Ekin, E_kin)
			@show E_int
			push!(Eint, E_int)
			@show E_g
			push!(Eg, E_g)
			@show E_wal	
			push!(Ewall, E_wal)
			E_err = E_tot - E0
			@show E_err	
			T = average_T(sys)
			push!(Ts, T)
			@show T
			println("# of part. = ", length(sys.particles))
			println()
		end
	end

	# Plotting the energies in time
	p = plot(times, Etot, label = "E_tot",legend=:bottomright)
	savefig(p, folder_name*"/Etot.pdf")
	p = plot(times, Ekin, label = "E_kin",legend=:bottomright)
	savefig(p, folder_name*"/Ekin.pdf")
	p = plot(times, Eint, label = "E_int",legend=:bottomright)
	savefig(p, folder_name*"/Eint.pdf")
	p = plot(times, Eg, label = "E_g",legend=:bottomright)
	savefig(p, folder_name*"/Eg.pdf")
	p = plot(times, Ewall, label = "E_wall",legend=:bottomright)
	savefig(p, folder_name*"/Ewall.pdf")
	p = plot(times, Ts, label = "T",legend=:bottomright)
	savefig(p, folder_name*"/T.pdf")

	df = DataFrame(time_steps = times, E_total = Etot, E_kinetic = Ekin, E_internal = Eint, E_walls = Ewall, E_g = Eg, temperature = Ts)
	CSV.write(folder_name*"/results.csv", df)

	final_T = average_T(sys)
	@show initial_T
	@show final_T

	save_pvd_file(out)
	fix_vtp_headers!(folder_name)

end ## function main

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

function average_T(sys::ParticleSystem)::Float64
    T = 0.0
    n = 0
    for p in sys.particles
        if p.type == FLUID
            T += p.T
            n += 1
        end
    end
    return T/n
end

end ## module
