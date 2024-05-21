module taco #TAylor-COuette flow
using SmoothedParticles
using Parameters
import LinearAlgebra
import StaticArrays
include("tools.jl")

#DECLARE CONSTANTS
#-----------------

const R1 = 1.0
const R2 = 2.0
const omega = 1.0
const Re = 20.0
const c_s = 30.0
const c_0 = 15.0
const rho0 = 1.0

const tau = 6*omega*R2*(R2-R1)/(Re*c_s^2)
const dr = (R2-R1)/20
const h = 3.0*dr
const wwall = 1.5h
const c_p = 0.01c_0
const c = sqrt(c_0^2 + 4/3*c_s^2)

const m0 = rho0*dr*dr

const dt = 0.05*dr/c
const t_end = 2.0
const dt_plot = max(t_end/20, dt)

@show dr
@show dt
@show c_0
@show c_s
@show c_p
@show h/2

function vexact(x::RealVector)::RealVector
    r = norm(x)
    return R2/r*(r/R1 - R1/r)/(R2/R1 - R1/R2)*RealVector(-omega*x[2], omega*x[1], 0.)
end

const inv = StaticArrays.inv
const det = StaticArrays.det

#PARTICLE FLAGS
const FLUID = 0.
const INNER = 1.
const OUTER = 2.

#DECLARE VARIABLES
#-----------------

@with_kw mutable struct Particle <: AbstractParticle
	m::Float64  = m0                  #mass
    x::RealVector                     #position
    x0::RealVector = x
    v::RealVector = VEC0  #velocity

    P::Float64 = 0.0      #pressure
    f::RealVector = VEC0  #force
    A::RealMatrix = MAT1  #distortion
    T::RealMatrix = MAT0  #stress
    L::RealMatrix = MAT0
    
    rho::Float64 = 0.
    lambda::Float64 = 0.

    C_rho::Float64 = 0.
    C_lambda::Float64 = 0.

    type::Float64
    err::Float64 = 0.
end

#CREATE INITIAL STATE
#--------------------

function make_geometry()
    grid = Grid(dr, :vogel)
    fluid = Circle(0., 0., R2) - Circle(0., 0., R1)
    walls = BoundaryLayer(fluid, grid, wwall)
    inner = Specification(walls, x -> norm(x) < 0.5(R1+R2))
    outer = Specification(walls, x -> norm(x) > 0.5(R1+R2))
    domain = BoundaryLayer(fluid, grid, 10wwall)
    sys = ParticleSystem(Particle, domain, h)
    generate_particles!(sys, grid, fluid, x -> Particle(x=x, type=FLUID))
    generate_particles!(sys, grid, inner, x -> Particle(x=x, type=INNER))
    generate_particles!(sys, grid, outer, x -> Particle(x=x, type=OUTER))
    create_cell_list!(sys)
    apply!(sys, find_rho!, self = true)
    for p in sys.particles
        p.C_rho = rho0 - p.rho
        p.C_lambda = -p.lambda
    end
    apply!(sys, reset!)
    apply!(sys, find_rho!, self = true)
    apply!(sys, find_T!)
    apply!(sys, find_f!)   
    return sys
end


#DECLARE PHYSICS
#---------------

function update_v!(p::Particle)
    if p.type == FLUID
        p.v += 0.5*dt*p.f/p.m
    else
        p.v = vexact(p.x)
    end
end

function update_x!(p::Particle, t::Float64)
    if p.type == FLUID
        p.x += 0.5*dt*p.v
    elseif p.type == OUTER
        p.x = RealVector(
            p.x0[1]*cos(omega*t) - p.x0[2]*sin(omega*t), 
            p.x0[1]*sin(omega*t) + p.x0[2]*cos(omega*t), 
            0.
        )
    end
end

function find_L!(p::Particle, q::Particle, r::Float64)
    ker = q.m*rDwendland2(h,r)
    x_pq = p.x-q.x
    v_pq = p.v-q.v
    p.T += ker*outer(x_pq, x_pq)
    p.L += ker*outer(v_pq, x_pq)
end

function update_A!(p::Particle)
    p.L = p.L*subinv(p.T)
    p.A = p.A*(MAT1 - 0.5*dt*p.L)*inv(MAT1 + 0.5*dt*p.L)
end

function find_rho!(p::Particle, q::Particle, r::Float64)
    x_pq = p.x-q.x
    p.T += q.m*rDwendland2(h,r)*outer(x_pq, x_pq)
    p.rho += q.m*wendland2(h,r)
    p.lambda += q.m*wendland2h(h,r)
end

function find_T!(p::Particle)
    G = transpose(p.A)*p.A
    p.P = c_0^2*(p.rho - rho0)*rho0/p.rho
    p.T = -p.P/p.rho^2*MAT1 + c_s^2*G*dev(G)*subinv(p.T)
end

function find_f!(p::Particle, q::Particle, r::Float64)
    ker = q.m*rDwendland2(h,r)
    kerh = q.m*rDwendland2h(h,r)
    x_pq = p.x-q.x
    #stress
    p.f += p.m*ker*(p.T + q.T)*x_pq
    #anti-clumping force
    p.f += -p.m*kerh*(c_p/rho0)^2*(p.lambda + q.lambda)*x_pq
end

function reset!(p::Particle)
    p.f = VEC0
    p.L = MAT0
    p.T = MAT0
    p.rho = p.C_rho
    p.lambda = p.C_lambda
end

function find_err!(p::Particle)
    p.err = norm(p.v - vexact(p.x))
end

function relax_f(A)
    return -3/tau*A*dev(transpose(A)*A)
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

#ENERGIES
#--------

function pE_kinetic(p::Particle)::Float64
    return 0.5*p.m*dot(p.v, p.v)
end 

function pE_bulk(p::Particle)::Float64
    return 0.5*p.m*c_0^2*(p.rho - rho0)^2/p.rho^2
end

function pE_shear(p::Particle)::Float64
    G = transpose(p.A)*p.A
    return 0.25*p.m*c_s^2*LinearAlgebra.norm(dev(G),2)^2
end

function pE_penalty(p::Particle)::Float64
    return 0.5*p.m*c_p^2*(p.lambda/rho0)^2
end


#TIME ITERATION
#--------------

function main()
    sys = make_geometry()
    out = new_pvd_file("results/taco")
    csv_data = open("results/taco/data.csv", "w")
    write(csv_data, string("t,E_total,E_kinetic,E_bulk,E_shear,E_penalty\n"))
    @time for k = 0 : Int64(round(t_end/dt))
        t = k*dt
        if (k % Int64(round(dt_plot/dt)) == 0)
            @show t
            N = length(sys.particles)
            @show N
            E_kinetic = sum(p -> pE_kinetic(p), sys.particles)
            E_bulk = sum(p -> pE_bulk(p), sys.particles)
            E_shear = sum(p -> pE_shear(p), sys.particles)
            E_penalty = sum(p -> pE_penalty(p), sys.particles)
            E_total = E_kinetic + E_bulk + E_shear + E_penalty
            apply!(sys, find_err!)
            err = sum(p -> p.err, sys.particles)/N
            @show E_total
            @show err
            println()
            write(csv_data, vec2string([t, err, E_total, E_kinetic, E_bulk, E_shear, E_penalty]))
            save_frame!(out, sys, :v, :A, :P, :type, :err)
        end
        apply!(sys, update_v!)
        apply!(sys, p -> update_x!(p,t + 0.5dt))
        create_cell_list!(sys)
        apply!(sys, reset!)
        apply!(sys, find_L!)
        apply!(sys, update_A!)
        apply!(sys, relax_A!)
        apply!(sys, p -> update_x!(p,t + dt))
        create_cell_list!(sys)
        apply!(sys, reset!)
        apply!(sys, find_rho!, self=true)
        apply!(sys, find_T!)
        apply!(sys, find_f!)        
        apply!(sys, update_v!)
    end
    save_pvd_file(out)
    close(csv_data)


    #save velocity distribution
    csv_data = open("results/taco/velocity.csv", "w")
    write(csv_data, string("s,v_ana,v_num\n"))
    s = range(R1, R2, length=50)
    for i in 1:length(s)
        x = RealVector(s[i], 0., 0.)
        gamma = SmoothedParticles.sum(sys, (p,r) -> m0/p.rho*wendland2(h,r), x)
        v_num = SmoothedParticles.sum(sys, (p,r) -> m0/p.rho*p.v[2]*wendland2(h,r), x)/gamma
        v_ana = vexact(x)[2]
        write(csv_data, vec2string([s[i], v_ana, v_num]))
    end
    close(csv_data)
end

end #module

