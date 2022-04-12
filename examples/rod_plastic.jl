module Rod
using SmoothedParticles
using Parameters
using Plots
using CSV
using DataFrames
import LinearAlgebra

const L = 5.0
const W = 0.5
const D = 1.0
const pull_D = 0.1
const pull_force = 40.0 #N/m^3
const pull_T = 0.4

const c_0 = 40.0
const c_s = 100.0
const c_l = sqrt(c_0^2 + 4/3*c_s^2)  #total sound speed
const rho0 = 1.0
const v_char = pull_force/rho0*pull_T
const tau = pull_T/200
const nu = 1e-3
const s_yield = 80

const dr = W/12
const h = 2.5dr
const vol = dr*dr
const m = rho0*vol

const dt = min(0.2*h/c_l, tau/3)
const t_end = 10.0
const dt_plot = max(t_end/400, dt)
const dt_frame = max(t_end/200, dt)

@with_kw mutable struct Particle <: AbstractParticle
	x::RealVector         #position
    v::RealVector = VEC0  #velocity
    f::RealVector = VEC0  #force
    X::RealVector = x     #Lag. position
    A::RealMatrix = MAT1  #distortion
    H::RealMatrix = MAT0  #density matrix
    T::RealMatrix = MAT0  #stress tensor
    DA::RealMatrix = MAT0 #tmp variable for A increments
    s::Float64 = 0.       #measure the amount of stress : if greater than s_Y, material undergoes permanent deformation
end

#ALGEBRAIC TOOLS

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
        0., 0., 0.
    )
end

@inbounds function trans(A::RealMatrix)::RealMatrix
    return RealMatrix(
        A[1,1], A[1,2], 0., 
        A[2,1], A[2,2], 0.,
        0.,  0., 0.
    )
end

@inbounds function dev(G::RealMatrix)::RealMatrix
    λ = 1/3*(G[1,1] + G[2,2] + 1.0)
    return RealMatrix(
        G[1,1] - λ,G[2,1], 0.0,
        G[1,2], G[2,2] - λ, 0.0,
        0.0, 0.0, 1.0 - λ
    )
end

#CREATE INITIAL STATE

function make_geometry()
    grid = Grid(dr, :hexagonal)
    rod = Rectangle(0., .0, L, W)
    dom = Rectangle(-D, -3*D, L+D, W+3*D)
    sys = ParticleSystem(Particle, dom, h)
    generate_particles!(sys, grid, rod, x -> Particle(x=x))
    create_cell_list!(sys)
    force_computation!(sys, 0.)
    return sys
end

function force_computation!(sys::ParticleSystem, t::Float64)
    apply!(sys, find_H!)
    apply!(sys, find_T!)
    apply!(sys, find_f!)
    if t < pull_T
        apply!(sys, pull!)
    end
end

#PHYSICS
#-------------------------------------

@inbounds function find_H!(p::Particle, q::Particle, r::Float64)
    ker = wendland2(h,r)
    x_pq = p.x-q.x
    v_pq = p.v-q.v
    p.H  += -ker*outer(x_pq, x_pq)
    p.DA +=  ker*outer(v_pq, x_pq)
end

@inbounds function find_T!(p::Particle)
    Hi = inv(p.H)
    p.A += dt*p.A*p.DA*Hi
    p.DA = MAT0
    At = trans(p.A)
    G = At*p.A
    P = c_0^2*(det(p.A)-1.0)
    p.T = m*(P*inv(At) + c_s^2*p.A*dev(G))*Hi
end

@inbounds function find_f!(p::Particle, q::Particle, r::Float64)
    ker = wendland2(h,r)
    #rDker = rDwendland2(h,r)
    x_pq = p.x-q.x
    #main force
    p.f += -ker*(trans(p.A)*(p.T*x_pq))
    p.f += -ker*(trans(q.A)*(q.T*x_pq))
end

function relax_f(A::RealMatrix)::RealMatrix
    return -3/tau*A*dev(A'*A)
end

@inbounds function relax_A!(p::Particle)
    if p.s > s_yield
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
end

@inbounds function artificial_viscosity!(p::Particle, q::Particle, r::Float64)
    rDker = rDwendland2(h,r)
    p.f += 2*m*vol*rDker*nu*(p.v-q.v)
end
    
@inbounds function update_v!(p::Particle)
    p.v += 0.5*dt*p.f/m
    if p.X[1] < h
        p.v = VEC0
    end
end

@inbounds function update_x!(p::Particle)
    p.x += dt*p.v
    #reset vars
    p.H = MAT0
    p.f = VEC0
    p.A += dt*p.DA #for applying artificial viscosity on A 
    p.DA = MAT0
end

@inbounds function find_s!(p::Particle)
    G = p.A'*p.A
    sigma = c_s^2*G*dev(G)
    p.s = sqrt(1.5)*LinearAlgebra.norm(dev(sigma))
end

function pull!(p::Particle)
    if p.X[1] > L - pull_D
        p.f += RealVector(0., vol*pull_force, 0.)
    end
    #=
    if p.X[1] < pull_D
        p.f -= RealVector(0., vol*pull_force, 0.)
    end
    =#
end

#Compute energy

function particle_energy(p::Particle)
    d = abs(det(p.A))
    G = trans(p.A)*p.A
    G0 = dev(G)
    E_kinet = 0.5*m*dot(p.v, p.v)
    E_shear = 0.25*m*c_s^2*LinearAlgebra.norm(G0,2)^2
    E_press = m*c_0^2*(d - 1.0 - log(d))
    return E_kinet + E_shear + E_press
end 

#Time iteration
#--------------


function main()
    sys = make_geometry()
    p_sel = argmax(p -> abs(p.x[1]) + abs(p.x[2]), sys.particles)
    out = new_pvd_file("results/rod_plastic")
    csv_data = open("rod_plastic/rod.csv", "w")
    @time for k = 0 : Int64(round(t_end/dt))
        t = k*dt
        if (k % Int64(round(dt_plot/dt)) == 0)
            println("t = ", k*dt)
            println("N = ", length(sys.particles))
            E = sum(p -> particle_energy(p), sys.particles)
            println("E = ", E)
            println("h = ", p_sel.x[2])
            println()
            write(csv_data, string(k*dt,",",p_sel.x[2],",",E,"\n"))
        end
        if (k % Int64(round(dt_frame/dt)) == 0)
            save_frame!(out, sys, :v, :A, :f, :T, :s)
        end
        apply!(sys, update_v!)
        apply!(sys, update_x!)
        apply!(sys, find_s!)
        apply!(sys, artificial_viscosity!)
        apply!(sys, relax_A!)
        create_cell_list!(sys)
        force_computation!(sys, t)
        apply!(sys, update_v!)
    end
    save_pvd_file(out)
    close(csv_data)
    #plot result
    data = CSV.read("results/rod_plastic/rod.csv", DataFrame; header=false)
    p1 = plot(data[:,1], data[:,2], label = "", xlabel = "time", ylabel = "amplitude")
    p2 = plot(data[:,1], data[:,3], label = "", xlabel = "time", ylabel = "energy")
    savefig(p1, "results/rod_plastic/amplitude.pdf")
    savefig(p2, "results/rod_plastic/energy.pdf")
end

end

