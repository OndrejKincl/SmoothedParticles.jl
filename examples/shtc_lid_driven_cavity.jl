module LDC
using SmoothedParticles
using Parameters
import LinearAlgebra
#using Plots
#using CSV
#using DataFrames
#include("make_graphs.jl")

const llid = 1.0                #length of the lid
const vlid = 1.0				#flow speed of the lid
const rho0 = 1.0                #density
const c_l = 20.0			#longitudinal speed of sound
const c_s = 40.0			#shear speed of sound
const Re = 100.0
const c_0 = sqrt(c_l^2 + 4/3*c_s^2)
const dr = llid/60	        #average particle distance
const tau = 6*rho0*vlid*llid/(Re*c_s^2)

const name = "ldc_c40"

#artificial_viscosities
const nu_v = 1e-6
const nu_A = 1e-6

const h = 2.4dr
const vol = dr*dr
const m = rho0*vol
const wwall = h
const acf = 1e-3

const dt = min(0.1*h/c_l, tau/3)
const t_end = 2.0
const t_acc = 0.4
const dt_plot = max(t_end/100, dt)
const dt_frame = max(t_end/100, dt)

#particle types
const FLUID = 0.
const WALL = 1.
const LID = 2.


@with_kw mutable struct Particle <: AbstractParticle
	x::RealVector         #position
    v::RealVector = VEC0  #velocity
    f::RealVector = VEC0  #force
    A::RealMatrix = MAT1  #distortion
    H::RealMatrix = MAT0  #density matrix
    T::RealMatrix = MAT0  #stress tensor
    DA::RealMatrix = MAT0
    type::Float64
    rho::Float64 = 0.
    rho0::Float64 = 0.
end

@inbounds function outer(x::RealVector, y::RealVector)::RealMatrix
    return RealMatrix(
        x[1]*y[1], x[2]*y[1], 0., 
        x[1]*y[2], x[2]*y[2], 0.,
        0., 0., 0.
    )
end

@inbounds function dev2(G::RealMatrix)::RealMatrix
    return G - 1/3*(G[1,1] + G[2,2] + G[3,3])*MAT1
end

function det2(A::RealMatrix)::Float64
    return (A[1,1]*A[2,2] - A[2,1]*A[1,2])*A[3,3]
end

function inv2(H::RealMatrix)::RealMatrix
    idet = 1.0/(H[1,1]*H[2,2] - H[2,1]*H[1,2])
    return RealMatrix(
         idet*H[2,2], -idet*H[2,1],     0.,
        -idet*H[1,2],  idet*H[1,1],     0.,
              0.,     0.,     0.,
    )
end

function compute_fluxes(sys::ParticleSystem, res = 100)
    s = range(0.,1.,length=res)
    v1 = zeros(res)
    v2 = zeros(res)
    fluxes = open("SHTC_LDC_Re"*string(Re)*"/fluxes.txt", "w")
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

#CREATE INITIAL STATE

function make_geometry()
    grid = Grid(dr, :hexagonal)
	box = Rectangle(0., 0., llid, llid)
	walls = BoundaryLayer(box, grid, wwall)
	lid   = Specification(walls, x -> x[2] > llid)
	walls = Specification(walls, x -> x[2] <= llid)
    sys = ParticleSystem(Particle, walls + lid + box, h)
	generate_particles!(sys, grid, box, x -> Particle(x=x, type=FLUID))
	generate_particles!(sys, grid, lid, x -> Particle(x=x, type=LID))
	generate_particles!(sys, grid, walls, x -> Particle(x=x, type=WALL))
    create_cell_list!(sys)
    apply!(sys, find_rho!, self = true)
    for p in sys.particles
        p.rho0 = p.rho
    end
    force_computation!(sys)
    return sys
end

function force_computation!(sys::ParticleSystem, t = 0.)
    apply!(sys, reset!)
    apply!(sys, p -> enforce_dbc!(p,t))
    apply!(sys, find_rho!, self = true)
    apply!(sys, find_H!)
    apply!(sys, find_T!)
    apply!(sys, find_f!)
    apply!(sys, artificial_viscosity!)
    apply!(sys, update_A!)
end

#PHYSICS
#-------------------------------------

@inbounds function reset!(p::Particle)
    p.H = MAT0
    p.f = VEC0
    p.rho = rho0 - p.rho0
    p.DA = MAT0
end

@inbounds function find_rho!(p::Particle, q::Particle, r::Float64)
    p.rho += m*wendland2(h,r)
end

@inbounds function find_H!(p::Particle, q::Particle, r::Float64)
    ker = wendland2(h,r)
    x_pq = p.x - q.x
    v_pq = p.v - q.v
    p.H += -ker*outer(x_pq, x_pq)
    p.DA += ker*outer(v_pq, x_pq)
end

function relax_func(A)
    return -3/tau*abs(det2(A))^(5/3)*A*dev2(A'*A)
end

@inbounds function relax_A!(p::Particle)
    #RK4 scheme
    A_new = p.A
    K = relax_func(p.A)
    A_new += dt*K/6
    K = relax_func(p.A + dt*K/2)
    A_new += dt*K/3
    K = relax_func(p.A + dt*K/2)
    A_new += dt*K/3
    K = relax_func(p.A + dt*K)
    A_new += dt*K/6
    p.A = A_new
end

@inbounds function find_T!(p::Particle)
    Hi = inv2(p.H)
    p.A += dt*p.A*p.DA*Hi
    G = p.A'*p.A
    p.T = m*c_s^2*G*dev2(G)*Hi
end

@inbounds function find_f!(p::Particle, q::Particle, r::Float64)
    x_pq = p.x-q.x
    P_p = c_l^2*(p.rho-rho0)/(1.0+acf)
    P_q = c_l^2*(q.rho-rho0)/(1.0+acf)
    p.f += -m^2*(P_p/p.rho^2 + P_q/q.rho^2)*rDwendland2(h,r)*x_pq
    p.f += -(p.T + q.T)*wendland2(h,r)*x_pq
end


@inbounds function update_v!(p::Particle)
    if p.type == FLUID
        p.v += 0.5*dt*p.f/m
    end
end

@inbounds function update_x!(p::Particle)
    if p.type == FLUID
        p.x += dt*p.v
    end
end

@inbounds function enforce_dbc!(p::Particle, t::Float64)
    if p.type == LID
        p.v = (t > t_acc ? 1.0 : t/t_acc)*vlid*VECX
    end
end

@inbounds function artificial_viscosity!(p::Particle, q::Particle, r::Float64)
    rDker = rDwendland2(h,r)
    p.f += 2*m*vol*rDker*nu_v*(p.v-q.v)
    p.DA += 2*vol*rDker*nu_A*(p.A-q.A)
end

function update_A!(p::Particle)
    p.A += dt*p.DA
    p.DA = MAT0
end

#Compute energy

function particle_energy(p::Particle)
    d = abs(p.rho/rho0)
    G = p.A'*p.A
    G0 = dev2(G)
    E_kinet = 0.5*m*dot(p.v, p.v)
    E_shear = 0.25*m*c_s^2*LinearAlgebra.norm(G0,2)^2
    E_press = m*c_l^2*(d - 1.0 - log(d))
    return E_kinet + E_shear + E_press
end 

#Time iteration
#--------------


function main()
    sys = make_geometry()
    out = new_pvd_file(name)
    @time for k = 0 : Int64(round(t_end/dt))
        t = k*dt
        if (k % Int64(round(dt_plot/dt)) == 0)
            println("t = ", k*dt)
            println("N = ", length(sys.particles))
            E = sum(p -> particle_energy(p), sys.particles)
            println("E = ", E)
            println()
        end
        if (k % Int64(round(dt_frame/dt)) == 0)
            save_frame!(out, sys, :v, :A, :type, :rho, :rho0, :f)
        end
        apply!(sys, update_v!)
        apply!(sys, update_x!)
        create_cell_list!(sys)
        force_computation!(sys, t)
        apply!(sys, update_v!)
        apply!(sys, relax_A!)
    end
    save_pvd_file(out)
    compute_fluxes(sys)
    #make_graphs.main([100.0])
end

end

