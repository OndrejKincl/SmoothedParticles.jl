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
const dr = chan_l/1000	     #average particle distance for the reference density
@show dr
const m = rho0*dr		#particle mass
#const drL = dr * sqrt(rho0 / rhoL) 		     #average particle distance (decrease to make finer simulation)
#const drR = dr * sqrt(rho0 / rhoR)		     #average particle distance (decrease to make finer simulation)
const drL = dr * rho0 / rhoL 		     #average particle distance (decrease to make finer simulation)
const drR = dr * rho0 / rhoR		     #average particle distance (decrease to make finer simulation)
@show drL
@show drR
const dxR = drR^2/drL # drR^2 = dxR * drL 
@show dxR
const h = 2.5*max(drL, drR)		     #size of kernel support
const wall_w = 2.2*dr        #width of the wall
const LRboundary = 0.3 #boundary position between L and R
const dr_wall = 0.95*drL
const E_wall = 10.0
const eps = 1e-6


#temporal parameters
const dt = 0.002*h/c      #time step
@show dt
const t_end = 0.2      #end of simulation
const dt_frame = 20*dt    #how often data is saved
const steps = t_end/dt
@show steps
#const number_of_profile_frames = 100
#const dt_profile = steps/number_of_profile_frames*dt    #how often data is saved
const dt_profile = dt_frame * 2
@show dt_frame

const bins = 10

#particle types
const FLUID = 0.
const WALL = 1.

@fastmath function wendland1(h::Float64, r::Float64)::Float64
	x = r/h
	if (x > 1.0)
		return 0.0
	end
	return 1.5*((1.0 - x)^4)*(1.0 + 4.0*x)/h
end

@fastmath function Dwendland1(h::Float64, r::Float64)::Float64
	x = r/h
	if (x > 1.0)
		return 0.0
	end
	return -30.0*x*((1.0 - x)^3)/h^2
end

@fastmath function rDwendland1(h::Float64, r::Float64)::Float64
	x = r/h
	if (x > 1.0)
		return 0.0
	end
	return -30.0*((1.0 - x)^3)/h^3
end

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
	#1D version
	line = Rectangle(-0.2*chan_l, 0., 1.4*chan_l, h)
	sys = ParticleSystem(Particle, line, h)
	i1 = Int64(round(LRboundary/drL))
	i2 = Int64(round(LRboundary/drR))
	i3 = Int64(round((chan_l-LRboundary)/drR))
	for i in 1:i1
		push!(sys.particles, Particle(RealVector(((i-1)*drL, 0.0, 0.0)), rhoL, FLUID))	
	end
	for i in i2:i3
		push!(sys.particles, Particle(RealVector((i*drR, 0.0, 0.0)), rhoR, FLUID))	
	end
	sys.particles[1].type = WALL
	sys.particles[2].type = WALL
	sys.particles[3].type = WALL
	N = length(sys.particles)
	sys.particles[N-2].type = WALL
	sys.particles[N-1].type = WALL
	sys.particles[N].type = WALL

    return sys
end

#Define interactions between particles

@inbounds function balance_of_mass!(p::Particle, q::Particle, r::Float64)
	ker = m*rDwendland1(h,r)
	#p.Drho += ker*(dot(p.x-q.x, p.v-q.v) + 2*nu*(p.rho-q.rho))
	p.Drho += ker*(dot(p.x-q.x, p.v-q.v))
end

function reset_rho!(p::Particle)
    p.rho = 0.0
end

@inbounds function find_rho0!(p::Particle, q::Particle, r::Float64)
    if p.type == FLUID && q.type == FLUID
		p.rho0 += m*wendland1(h,r)
	end
end

@inbounds function find_rho!(p::Particle, q::Particle, r::Float64)
    if p.type == FLUID && q.type == FLUID
		p.rho += m*wendland1(h,r)
	end
end

function find_pressure!(p::Particle)
	p.rho += p.Drho*dt
	p.Drho = 0.0
	#p.P = rho0*c^2*((p.rho/rho0)^gamma - 1.0)/gamma #Tait equation
    #p.P = c^2 * (p.rho-p.rho0)/gamma
    p.P = c^2 * p.rho/gamma
end

@inbounds function internal_force!(p::Particle, q::Particle, r::Float64)
	if p.type == FLUID && q.type == FLUID
		ker = m*rDwendland1(h,r)
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

function calculate_profile_1D(field:: Function, xmin:: Float64, xmax:: Float64, bins:: Int64)::Array{Float64,1}
	profile = zeros(bins)	
	dx = xmax/bins
	for i in 1:bins
		#println("Integrating over ", ((i-1)*dx, ymin), ", ", (i*dx, ymax))
		if i % Int64(round(bins/10)) == 0 
			println(Int64(round(i/bins*100)), "%")
		end
		(val,err) = hquadrature(field, (i-1)*dx, i*dx;	reltol=1e-8, abstol=0, maxevals=0)	
		profile[i] = val / dx 
	end
	return profile
end



function export_profile(time:: Float64, profile:: Array{Float64,1}, csv_file::IO)::String
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
	normalization = SmoothedParticles.sum(sys, (p,r) -> wendland1(h,r), x)
	return SmoothedParticles.sum(sys, (p,r) -> wendland1(h,r)*p.rho, x)/normalization
end

function pressure_field(sys::ParticleSystem, x::RealVector)::Float64
	normalization = SmoothedParticles.sum(sys, (p,r) -> wendland1(h,r), x)
	return SmoothedParticles.sum(sys, (p,r) -> wendland1(h,r)*p.P, x)/normalization
end

function velocity_field(sys::ParticleSystem, x::RealVector)::Float64
	normalization = SmoothedParticles.sum(sys, (p,r) -> wendland1(h,r), x)
	return SmoothedParticles.sum(sys, (p,r) -> wendland1(h,r)*p.v[1], x)/normalization
end


function  main(find_density_profile = false, find_pressure_profile = false, find_velocity_profile = false)
    sys = make_system()
    out = new_pvd_file(folder_name)
    if find_density_profile
		csv_density = open(string(folder_name,"/density.csv"), "w")
	end
	if find_pressure_profile
    	csv_pressure = open(string(folder_name,"/pressure.csv"), "w")
	end
	if find_velocity_profile
    	csv_velocity= open(string(folder_name,"/velocity.csv"), "w")
	end

    apply!(sys, find_rho0!, self = true)

	xs = [i*drR for i in 1:Int64(round(chan_l/drR))]

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
			println("N = ", length(sys.particles))
            save_frame!(out, sys, :v, :P, :rho, :type)
        end

        if (k %  Int64(round(dt_profile/dt)) == 0) || (k==Int64(round(t_end/dt)))
			if find_density_profile
				println("Calculating density profile.")
				#density_field = interpolate_field2D(sys, p -> p.rho, wendland1, h)
				#density_profile = calculate_profile(x->density_field(sys, RealVector((x[1],x[2],0.0))), 0.0, chan_l, 0.0, chan_w, bins)
				#density_profile = calculate_profile(x->density_field(sys, RealVector((x[1],x[2],0.0))), 0.0, chan_l, 0.0, h, bins)
				#export_profile(t, density_profile, csv_density)
				#plot(density_profile)
				ys = zeros(length(xs))
				for i in 1:length(ys)
					ys[i] = density_field(sys, RealVector((xs[i],0.0,0.0))) 
				end
				export_profile(t, ys, csv_density)
				plot(xs, ys, labels = string("t=", k*dt), ylims = (0., 1.))
				gui()
			end
			if find_pressure_profile
				println("Calculating pressure profile.")
				#pressure_profile = calculate_profile(x->pressure_field(sys, RealVector((x[1],x[2],0.0))), 0.0, chan_l, 0.0, chan_w, bins)
				#pressure_profile = calculate_profile(x->pressure_field(sys, RealVector((x[1],x[2],0.0))), 0.0, chan_l, 0.0, h, bins)
				#export_profile(t, pressure_profile, csv_pressure)
				#plot(pressure_profile)
				ys = zeros(length(xs))
				for i in 1:length(ys)
					ys[i] = pressure_field(sys, RealVector((xs[i],0.0,0.0))) 
				end
				export_profile(t, ys, csv_pressure)
				plot(xs, ys, labels = string("t=", k*dt), ylims = (0., 1.))
				gui()
			end
			if find_velocity_profile
				println("Calculating velocity profile.")
				ys = zeros(length(xs))
				for i in 1:length(ys)
					ys[i] = velocity_field(sys, RealVector((xs[i],0.0,0.0))) 
				end
				export_profile(t, ys, csv_velocity)
				plot(xs, ys, labels = string("t=", k*dt), ylims = (0., 2.))
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
	if find_velocity_profile
    	close(csv_velocity)
	end
end

function animate(csv_file_name::String, gif_file_name::String, field_name::String, field_label::String)
	println("Reading from ", string(folder_name, "/", csv_file_name))
	data = Matrix(CSV.read(string(folder_name, "/", csv_file_name), DataFrame; header=false))
	line_length = length(data[1,:])
	rows = length(data[:,1])
	xs = [i*chan_l/(line_length-1) for i in 2:line_length] #assuming equidistant
	#p1 = plot(xs, data[1,:], label = "density", xlabel = "x", ylabel = "rho")
	#savefig(p1, string(folder_name, "/", "plot.pdf"))

	anim = @animate for i in 1:rows
    		plot(xs, data[i,2:line_length], label = string("t=",data[i,1]), xlabel = "x", ylabel = field_name, ylims = (0., 1.))
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

function animate_velocity()
	animate("velocity.csv", "velocity.gif", "velocity", "v")
end

function animate_all()
	animate_density()
	animate_pressure()
	animate_velocity()
end

end
