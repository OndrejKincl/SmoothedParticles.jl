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
using Cubature
using Plots
using CSV
using DataFrames
#include("utils/ICR.jl")

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
const gamma = 1.0
#csL = sqrt(gamma * pL / rhoL) #speed of sound left
#csR = sqrt(gamma * pR / rhoR) #speed of sound right
#@show csL
#@show csR
#const c = max(csL, csR)	#numerical speed of sound
const c = 1.0
#const mu = 1.0e-3		#dynamic viscosity of water
#const nu = 1.0e-3		#pressure stabilization
#const P0 = 1.2          #anti-clump term

#geometry parameters
const chan_l = 1.0           #length of the channel
const chan_w = 0.01     #width of the channel
const dr = chan_w / 40 		     #average particle distance for the reference density
@show dr
const m = rho0*dr^2		#particle mass
const drL = dr * sqrt(rho0 / rhoL) 		     #average particle distance (decrease to make finer simulation)
const drR = dr * sqrt(rho0 / rhoR)		     #average particle distance (decrease to make finer simulation)
@show drL
@show drR
const dxR = drR^2/drL # drR^2 = dxR * drL 
@show dxR
const h = 2.5*max(drL, drR)		     #size of kernel support
const wall_w = 2.5*dr        #width of the wall
const LRboundary = 0.3 #boundary position between L and R
const dr_wall = 0.95*dr
const E_wall = 1.0
const eps = 1e-6


#temporal parameters
const dt = 0.2*h/c      #time step
@show dt
const t_end = 0.02      #end of simulation
const dt_frame = dt    #how often data is saved
const steps = t_end/dt
@show steps
const number_of_profile_frames = 10
const dt_profile = steps/number_of_profile_frames*dt    #how often data is saved
@show dt_frame

const bins = 10

#particle types
const FLUID = 0.
const WALL = 1.

#=
Declare variables to be stored in a Particle
=#

#function LJ_potential(p::Particle, q::Particle, r::Float64)::Float64
#	if q.type == WALL && p.type == FLUID && r < dr_wall
#		s2 = (dr_wall^2 + eps^2)/(r^2 + eps^2)
#		return m*E_wall*(0.25*s2^2 - 0.5*s2 + 0.25)
#	else
#		return 0.0
#	end
#end

mutable struct Particle <: AbstractParticle
    x::RealVector #position
    v::RealVector #velocity
    a::RealVector #acceleration
    rho::Float64 #density
    rho0::Float64 #initial density
    Drho::Float64 #rate of density
    P::Float64 #pressure
    type::Float64 #particle type
    Particle(x,type) = begin
        return new(x, VEC0, VEC0, rho0, rho0, 0., 0., type)
    end
    Particle(x, rho, type) = begin
        return new(x, VEC0, VEC0, rho, rho, 0., 0., type)
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

	wall = Specification(wall, x -> (x[2] > chan_w || x[2] < 0.0 || x[1] < 0.0))

	generate_particles!(sys, gridL, boxL, x -> Particle(x, rhoL, FLUID))
	generate_particles!(sys, gridR, boxR, x -> Particle(x, rhoR, FLUID))
	#ICR.renormalize!(sys, dr)

	#for i in 1:(Int64(round(chan_w/drL))-1)
	#	y = i * drL
	#	for j in 1:Int64(round((chan_l-LRboundary)/dxR))
	#		x = LRboundary + j*dxR
	#		push!(sys.particles, Particle(RealVector((x, y, 0.0)), FLUID))	
	#	end
	#end
	#ICR.renormalize!(sys, dr)

	generate_particles!(sys, grid, wall, x -> Particle(x, WALL))

    return sys
end

#Define interactions between particles

@inbounds function balance_of_mass!(p::Particle, q::Particle, r::Float64)
	ker = m*rDwendland2(h,r)
	#p.Drho += ker*(dot(p.x-q.x, p.v-q.v) + 2*nu*(p.rho-q.rho))
	p.Drho += ker*(dot(p.x-q.x, p.v-q.v))
end

function reset_rho!(p::Particle)
    p.rho = 0.0
end

@inbounds function find_rho0!(p::Particle, q::Particle, r::Float64)
    if p.type == FLUID && q.type == FLUID
		p.rho0 += m*wendland2(h,r)
	end
end

@inbounds function find_rho!(p::Particle, q::Particle, r::Float64)
    if p.type == FLUID && q.type == FLUID
		p.rho += m*wendland2(h,r)
	end
end

function find_pressure!(p::Particle)
	p.rho += p.Drho*dt
	p.Drho = 0.0
	#p.P = rho0*c^2*((p.rho/rho0)^7 - 1.0)/7
    #p.P = c^2 * (p.rho-p.rho0)/gamma
    p.P = c^2 * p.rho/gamma
end

@inbounds function internal_force!(p::Particle, q::Particle, r::Float64)
	if p.type == FLUID && q.type == FLUID
		ker = m*rDwendland2(h,r)
		p.a += -ker*(p.P/p.rho^2 + q.P/p.rho^2)*(p.x - q.x)
		#p.a += +2*ker*mu/rho0^2*(p.v - q.v)
	elseif p.type == FLUID && q.type == WALL && r < dr_wall
		s2 = (dr_wall^2 + eps^2)/(r^2 + eps^2)
		p.a += -E_wall/(r^2 + eps^2)*(s2 - s2^2)*(p.x - q.x)
	end	
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

function interpolate_field2D(sys::ParticleSystem, property:: Function, kernel::Function, h::Float64)::Function
	result = function (x,y)
		value = 0.0
		for	p in sys.particles
			value += kernel(h, sqrt((p.x[1]-x)^2 + (p.x[2]-y)^2)) * property(p)
		end	
		return value
	end
	return result
end

function calculate_profile(field:: Function, xmin:: Float64, xmax:: Float64, ymin:: Float64, ymax:: Float64, bins:: Int64)::Array{Float64,1}
	profile = zeros(bins)	
	dx = xmax/bins
	for i in 1:bins
		#println("Integrating over ", ((i-1)*dx, ymin), ", ", (i*dx, ymax))
		if i % Int64(round(bins/10)) == 0 
			println(Int64(round(i/bins*100)), "%")
		end
		(val,err) = hcubature(field, ((i-1)*dx, ymin), (i*dx, ymax);	reltol=1e-8, abstol=0, maxevals=0)	
		profile[i] = val / (dx * (ymax-ymin))
	end
	return profile
end

function export_profile(time:: Float64, profile:: Array{Float64, 1}, csv_file::IO)::String
	line = string(time, ", ")
	for i in 1:(length(profile)-1)
		value = profile[i]
		if isnan(value)
			value = 0.0
		end
		line = string(line, value, ", ")
	end
	line = string(line, string(profile[length(profile)]), "\n")
	write(csv_file, line)
	return line
end

function density_field(sys::ParticleSystem, x::RealVector)::Float64
	normalization = SmoothedParticles.sum(sys, (p,r) -> wendland2(h,r), x)
	return SmoothedParticles.sum(sys, (p,r) -> wendland2(h,r)*p.rho, x)/normalization
end

function pressure_field(sys::ParticleSystem, x::RealVector)::Float64
	normalization = SmoothedParticles.sum(sys, (p,r) -> wendland2(h,r), x)
	return SmoothedParticles.sum(sys, (p,r) -> wendland2(h,r)*p.P, x)/normalization
end



function  main(find_density_profile = false, find_pressure_profile = false)
    sys = make_system()
    out = new_pvd_file(folder_name)
    if find_density_profile
		csv_density = open(string(folder_name,"/density.csv"), "w")
	end
	if find_pressure_profile
    	csv_pressure = open(string(folder_name,"/pressure.csv"), "w")
	end

    apply!(sys, find_rho0!, self = true)

    #a modified Verlet scheme
    for k = 0 : Int64(round(t_end/dt))
        t = k*dt
        apply!(sys, move!)

        create_cell_list!(sys)
	apply!(sys, balance_of_mass!)
    	#apply!(sys, reset_rho!)
    	#apply!(sys, find_rho!, self = true)
        apply!(sys, find_pressure!)
        apply!(sys, internal_force!)
        apply!(sys, accelerate!)
        #save data at selected frames
        if (k %  Int64(round(dt_frame/dt)) == 0)
            @show t
            save_frame!(out, sys, :v, :P, :rho, :type)
        end

        if (k %  Int64(round(dt_profile/dt)) == 0) || (k==Int64(round(t_end/dt)))
		if find_density_profile
			println("Calculating density profile.")
			#density_field = interpolate_field2D(sys, p -> p.rho, wendland2, h)
			density_profile = calculate_profile(x->density_field(sys, RealVector((x[1],x[2],0.0))), 0.0, chan_l, 0.0, chan_w, bins)
			export_profile(t, density_profile, csv_density)
			plot(density_profile)
			gui()
		end
		if find_pressure_profile
			println("Calculating pressure profile.")
			pressure_profile = calculate_profile(x->pressure_field(sys, RealVector((x[1],x[2],0.0))), 0.0, chan_l, 0.0, chan_w, bins)
			export_profile(t, pressure_profile, csv_pressure)
			plot(pressure_profile)
			gui()
		end
	end
	apply!(sys, accelerate!)
    end
    save_pvd_file(out)

    if find_density_profile
		close(csv_density)
	end
	if find_pressure_profile
    	close(csv_pressure)
	end
end

function animate(csv_file_name::String, gif_file_name::String, field_name::String, field_label::String)
	println("Reading from ", string(folder_name, "/", csv_file_name))
	data = Matrix(CSV.read(string(folder_name, "/", csv_file_name), DataFrame; header=false))
	line_length = length(data[1, :])
	rows = length(data[:,1])
	xs = [i*chan_l/line_length for i in 1:line_length]
	#p1 = plot(xs, data[1,:], label = "density", xlabel = "x", ylabel = "rho")
	#savefig(p1, string(folder_name, "/", "plot.pdf"))

	anim = @animate for i in 1:rows
    		plot(xs, data[i,:], label = field_name, xlabel = "x", ylabel = field_label)
	end
	gif(anim, string(folder_name, "/", gif_file_name), fps = 15)
	gui()
end

function animate_density()
	animate("density.csv", "density.gif", "density", "rho")
end

function animate_pressure()
	animate("pressure.csv", "pressure.gif", "pressure", "P")
end

end
