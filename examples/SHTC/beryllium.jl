module beryllium
using SmoothedParticles
using Parameters
import LinearAlgebra
import StaticArrays
#using Plots
using CSV
using DataFrames

#DECLARE CONSTANTS
#-----------------

const L = 0.06
const W = 0.01
const c_s = 9046.59
const c_0 = c_s
const c_p = 4c_0
const c = sqrt(c_0^2 + 4/3*c_s^2)

const rho0 = 1845.0

const dr = W/20 #W/40
const h = 3.0dr
const m0 = rho0*dr*dr

const dt = 0.05*dr/c
const t_end = 3e-5
const dt_plot = max(t_end/200, dt)

function init_velocity(x::RealVector)::RealVector
    A = 4.3369e-5
    omega = 2.3597e5
    alpha = 78.834
    a1 = 56.6368
    a2 = 57.6455
    s = alpha*(x[1] + L/2)
    v = A*omega*(a1*(sinh(s) + sin(s)) - a2*(cosh(s) + cos(s)))
    return v*VECY
end

#STRUCTURAL KERNELS
#------------------

@fastmath function wendland2h(h::Float64, r::Float64)::Float64
    x = r/h
    return x < 1.0 ? 14.0*(1.0 - x)^3*(14.0*x^2 - 3.0*x - 1.0)/(pi*h^2) : 0.
end

@fastmath function rDwendland2h(h::Float64, r::Float64)::Float64
    x = r/h
    return x < 1.0 ? 140.0*(1.0 - x)^2*(4.0 - 7.0*x)/(pi*h^4) : 0.
end

#DECLARE VARIABLES
#-----------------

@with_kw mutable struct Particle <: AbstractParticle
	m::Float64  = m0 #mass
    x::RealVector         #position
    v::RealVector = init_velocity(x)  #velocity

    P::Float64 = 0.0 #pressure
    f::RealVector = VEC0  #force
    A::RealMatrix = MAT1  #distortion
    T::RealMatrix = MAT0  #stress
    L::RealMatrix = MAT0
    
    J::Float64 = 0.
    K::Float64 = 0.

    J0::Float64 = 0.
    K0::Float64 = 0.
end


#ALGEBRAIC TOOLS
#---------------

@inbounds function outer(x::RealVector, y::RealVector)::RealMatrix
    return RealMatrix(
        x[1]*y[1], x[2]*y[1], 0., 
        x[1]*y[2], x[2]*y[2], 0.,
        0., 0., 0.
    )
end

@inbounds function det(A::RealMatrix)::Float64
    return A[1,1]*A[2,2] - A[1,2]*A[2,1]
end

@inbounds function inv(A::RealMatrix)::RealMatrix
    idet = 1.0/det(A)
    return RealMatrix(
        +idet*A[2,2], -idet*A[2,1], 0., 
        -idet*A[1,2], +idet*A[1,1], 0.,
        0., 0., 1.
    )
end

@inbounds function dev(G::RealMatrix)::RealMatrix
    tr = 1/3*(G[1,1] + G[2,2] + G[3,3])
    return G - tr*MAT1
end


#CREATE INITIAL STATE
#--------------------

function make_geometry()
    grid = Grid(dr, :hexagonal)
    rod = Rectangle(-L/2, -W/2, L/2, W/2)
    dom = BoundaryLayer(rod, grid, W)
    sys = ParticleSystem(Particle, dom, h)
    generate_particles!(sys, grid, rod, x -> Particle(x=x))
    create_cell_list!(sys)
    apply!(sys, find_J!)
    for p in sys.particles
        p.J0 = 1.0 - p.J
        p.K0 = -p.K
    end
    apply!(sys, reset!)
    apply!(sys, find_J!)
    apply!(sys, find_T!)
    apply!(sys, find_f!)   
    return sys
end


#DECLARE PHYSICS
#---------------

function update_v!(p::Particle)
    p.v += 0.5*dt*p.f/p.m
end

function update_x!(p::Particle)
    p.x += 0.5*dt*p.v
end

function find_L!(p::Particle, q::Particle, r::Float64)
    ker = q.m/rho0*rDwendland2(h,r)
    x_pq = p.x-q.x
    v_pq = p.v-q.v
    p.T += ker*outer(x_pq, x_pq)
    p.L += ker*outer(v_pq, x_pq)
end

function update_A!(p::Particle)
    p.L = p.L*inv(p.T)
    p.A = p.A*(MAT1 - 0.5*dt*p.L)*StaticArrays.inv(MAT1 + 0.5*dt*p.L)
end

function find_J!(p::Particle, q::Particle, r::Float64)
    x_pq = p.x-q.x
    p.T += q.m/rho0*rDwendland2(h,r)*outer(x_pq, x_pq)
    p.J += q.m/rho0*wendland2(h,r)
    p.K += q.m/rho0*wendland2h(h,r)
end

function find_T!(p::Particle)
    G = transpose(p.A)*p.A
    p.P = 0.5*rho0*c_0^2*((1.0-1.0/p.J)/p.J^2 + log(p.J)/p.J)
    p.T = p.P/rho0*MAT1 - c_s^2*G*dev(G)*inv(p.T)
end

function find_f!(p::Particle, q::Particle, r::Float64)
    ker = q.m/rho0*rDwendland2(h,r)
    kerh = q.m/rho0*rDwendland2h(h,r)
    x_pq = p.x-q.x
    #stress
    p.f += -p.m*ker*(p.T*x_pq)
    p.f += -p.m*ker*(q.T*x_pq)
    #anti-clumping force
    p.f += -p.m*kerh*c_p^2*(p.K + q.K)*x_pq
end

function reset!(p::Particle)
    p.f = VEC0
    p.L = MAT0
    p.L = MAT0
    p.T = MAT0
    p.J = p.J0
    p.K = p.K0
end

#ENERGY COMPUTATION
#------------------

function pE_kinetic(p::Particle)::Float64
    return 0.5*p.m*dot(p.v, p.v)
end 

function pE_bulk(p::Particle)::Float64
    return 0.25*p.m*c_0^2*((1.0 - 1.0/p.J)^2 + log(p.J)^2)
end

function pE_shear(p::Particle)::Float64
    G = transpose(p.A)*p.A
    return 0.25*p.m*c_s^2*LinearAlgebra.norm(dev(G),2)^2
end

function pE_penalty(p::Particle)::Float64
    return 0.5*p.m*c_p^2*p.K^2
end


#TIME ITERATION
#--------------

function main()
    sys = make_geometry()
    centre = argmin(p -> LinearAlgebra.norm(p.x), sys.particles)
    out = new_pvd_file("results/beryllium")
    csv_data = open("results/beryllium/data.csv", "w")
    @time for k = 0 : Int64(round(t_end/dt))
        t = k*dt
        if (k % Int64(round(dt_plot/dt)) == 0)
            println("t = ", t)
            println("N = ", length(sys.particles)) 
            y = centre.x[2]
            E_kinetic = sum(p -> pE_kinetic(p), sys.particles)
            E_bulk = sum(p -> pE_bulk(p), sys.particles)
            E_shear = sum(p -> pE_shear(p), sys.particles)
            E_penalty = sum(p -> pE_penalty(p), sys.particles)
            E_total = E_kinetic + E_bulk + E_shear + E_penalty
            @show E_total
            @show y
            println()
            write(csv_data, vec2string([t, y, E_total, E_kinetic, E_bulk, E_shear, E_penalty]))
            save_frame!(out, sys, :v, :A, :P)
        end
        apply!(sys, update_v!)
        apply!(sys, update_x!)
        create_cell_list!(sys)
        apply!(sys, reset!)
        apply!(sys, find_L!)
        apply!(sys, update_A!)
        apply!(sys, update_x!)
        create_cell_list!(sys)
        apply!(sys, reset!)
        apply!(sys, find_J!)
        apply!(sys, find_T!)
        apply!(sys, find_f!)        
        apply!(sys, update_v!)
    end
    save_pvd_file(out)
    close(csv_data)
end

function vec2string(a::Vector)::String
    out = ""
    for i in 1:length(a)-1
        out = out*string(a[i])*","
    end
    if length(a) > 0
        out = out*string(a[end])
    end
    out = out*"\n"
end

end #module
