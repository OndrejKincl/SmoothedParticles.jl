# lid-driven cavity flow benchmark using SHTC
module SHTC_ldc

using SmoothedParticles
using LinearAlgebra
using StaticArrays
using Parameters
#using Plots
#using DifferentialEquations
#include("SPH_AVG.jl")
#include("make_graphs.jl")


#Declare const parameters
#------------------------
#all dims in SI

const llid = 1.0                #length of the lid
const vlid = 1.0				#flow speed of the lid
const rho0 = 1.0                #density
const c_l = 20.0			#longitudinal speed of sound
const c_s = 20.0			#shear speed of sound
const Re = 100.0
const tau = 6*vlid*llid/(Re*c_s^2)

const dr = llid/100 	        #average particle distance
const h = 2.4*dr		        #size of kernel support
const m = rho0*dr^2             #particle mass
const acf = 1e-3				#anti-clump factor
#const spf = 30					#number of steps per one Shepard filter
#const fst = 1e-2					#filter strength
const wwall = 1.5h

#temporal parameters
const dt = 0.05*h/c_l
const t_end = 1.0        			  #end of simulation
const dt_frame = max(t_end/200, dt) 	  #how often save data

#particle types
const FLUID = 0.
const WALL = 1.
const LID = 2.

#Declare fields
#--------------


@with_kw mutable struct Particle <: AbstractParticle
	x::RealVector
    v::RealVector = VEC0
    rho::Float64 = rho0
    type::Float64
    A::RealMatrix = MAT1
    stress::RealMatrix = MAT0
end

function compute_fluxes(sys::ParticleSystem, res = 100)
    s = range(0.,1.,length=res)
    v1 = zeros(res)
    v2 = zeros(res)
    fluxes = open("results/SHTC_LDC_Re"*string(Re)*"/fluxes.txt", "w")
    for i in 1:res
		#x-velocity along y-centerline
		x = RealVector(0.5, s[i], 0.)
		gamma = SmoothedParticles.sum(sys, (p,r) -> m/p.rho*wendland2(h,r), x)
        v1[i] = SmoothedParticles.sum(sys, (p,r) -> m/p.rho*p.v[1]*wendland2(h,r), x)/gamma
		#y-velocity along x-centerline
		x = RealVector(s[i], 0.5, 0.)
		gamma = SmoothedParticles.sum(sys, (p,r) -> m/p.rho*wendland2(h,r), x)
        v2[i] = SmoothedParticles.sum(sys, (p,r) -> m/p.rho*p.v[2]*wendland2(h,r), x)/gamma
		#save results into csv
        write(fluxes, string(s[i])*" "*string(v1[i])*" "*string(v2[i])*"\n")
    end
    close(fluxes)
end

function main()
#Define geometry and create particles
#------------------------------------
    grid = Grid(dr, :hexagonal)
	box = Rectangle(0., 0., llid, llid)
	walls = BoundaryLayer(box, grid, wwall)
	lid   = Specification(walls, x -> x[2] > llid)
	walls = Specification(walls, x -> x[2] <= llid)
    sys = ParticleSystem(Particle, walls + lid + box, h)
	generate_particles!(sys, grid, box, x -> Particle(x=x, type=FLUID))
	generate_particles!(sys, grid, lid, x -> Particle(x=x, type=LID, v = vlid*VECX))
	generate_particles!(sys, grid, walls, x -> Particle(x=x, type=WALL))

#utility functions
#-----------------

	@inbounds function deviatoric(G::RealMatrix)::RealMatrix
		return G - 1/3*(G[1,1] + G[2,2] + G[3,3])*MAT1
	end

#Define interactions between particles
#-------------------------------------
	@inbounds function update_rho!(p::Particle, q::Particle, r::Float64)
		if p.type == FLUID
            p.rho += dt*m*rDwendland2(h,r)*SmoothedParticles.dot(p.x-q.x, p.v-q.v)
        end
    end

	@inbounds function convect_A!(p::Particle, q::Particle, r::Float64)
		if p.type != LID
			p.A += dt*m/p.rho*rDwendland2(h,r)*p.A*((p.v - q.v)*(p.x - q.x)')
		end
	end

	function relax_f(A)
		return -3/tau*A*deviatoric(A'*A)
	end

    @inbounds function relax_A!(p::Particle)
		#RK4 scheme
		A_new = p.A
		K = relax_f(p.A)
		A_new += dt*K/6
		K = relax_f(p.A + dt*K/2)
		A_new += dt*K/3
		K = relax_f(p.A + dt*K/2)
		A_new += dt*K/3
		K = relax_f(p.A + dt*K)
		A_new += dt*K/6
		p.A = A_new
	end

    @inbounds function find_stress!(p::Particle)
        finger = p.A'*p.A
        p.stress = c_l^2*(p.rho - rho0/(1.0 + acf))*MAT1 + c_s^2*p.rho*finger*deviatoric(finger)
    end

	@inbounds function update_v!(p::Particle, q::Particle, r::Float64)
		if p.type == FLUID
			p.v += -dt*m*rDwendland2(h,r)*(p.stress/p.rho^2 + q.stress/q.rho^2)*(p.x - q.x)
		end
	end

	@inbounds function move!(p::Particle)
		if p.type == FLUID
			p.x += p.v*dt
		end
	end

	@inbounds function find_gamma!(p::Particle, q::Particle, r::Float64)
		p.gamma += m/q.rho*wendland2(h,r)
	end

	@inbounds function interpolate_A!(p::Particle, q::Particle, r::Float64)
		p.stress += m/q.rho*wendland2(h,r)*p.A
	end

	@inbounds function filter_A!(p::Particle, q::Particle, r::Float64)
		p.stress += m/q.rho*wendland2(h,r)*p.A
	end

#Time iteration
#--------------
	out = new_pvd_file("results/SHTC_LDC_Re"*string(Re))
	for k = 0 : Int64(round(t_end/dt))
		if (k % Int64(round(dt_frame/dt)) == 0)
			println("t = ", k*dt)
            println("N = ", length(sys.particles))
            println()
			save_frame!(out, sys, :v, :A, :rho, :type)
		end
		create_cell_list!(sys)
        apply!(sys, find_stress!)
        apply!(sys, update_v!)
        apply!(sys, update_rho!)
		apply!(sys, convect_A!)
		apply!(sys, relax_A!)
		apply!(sys, move!)
		#=shepard filter
		if k % spf == 0
			SPH_AVG.average!(sys, (p,q,r) -> m/q.rho*wendland2(h,r), :A, :stress, fst)
		end
		=#
	end
	save_pvd_file(out)
    compute_fluxes(sys)
	#make_graphs.main([100.0])
end

end