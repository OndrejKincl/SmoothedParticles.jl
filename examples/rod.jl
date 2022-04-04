module rod
using SmoothedParticles
using Parameters
using Plots
using CSV
using DataFrames
import LinearAlgebra

#CONSTANT PARAMETERS
#-------------------------------

const L = 5.0   #rod length
const W = 0.5   #rod width
const r_free = 1.0   #how much free space we want around the rod

const pull_force = 1.0 #pulling force [N]
const pull_time = 0.5  #for how long we pull

const c_l = 20.0   #longitudinal sound speed
const c_s = 200.0  #shear sound speed
const c_0 = sqrt(c_l^2 + 4/3*c_s^2)  #total sound speed
const rho0 = 1.0   #density
const nu = 1.0e-4    #artificial viscosity (surpresses noise but is not neccessary)

const dr = W/16    #discretization step
const h = 2.5dr    #support radius
const vol = dr^2   #particle volume
const m = rho0*vol #particle mass

const dt = 0.1h/c_0 #time step
const t_end = 5.0  #total simulation time
const dt_plot = max(t_end/400, dt) #how often save txt data (cheap)
const dt_frame = max(t_end/100, dt) #how often save pvd data (expensive)

#ALGEBRAIC TOOLS
#----------------------------------------

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

#DEFINE VARIABLES
#------------------------------

@with_kw mutable struct Particle <: AbstractParticle
	x::RealVector         #position
    v::RealVector = VEC0  #velocity
    f::RealVector = VEC0  #force
    X::RealVector = x     #Lag. position
    A::RealMatrix = MAT0  #distortion
    H::RealMatrix = MAT0  #correction matrix
    B::RealMatrix = MAT0  #derivative of energy wrt A
    e::Float64 = 0.       #fronorm squared of eta
end

#CREATE INITIAL STATE
#----------------------------

function make_geometry()
    grid = Grid(dr, :hexagonal)
    rod = Rectangle(0., .0, L, W)
    dom = Rectangle(-r_free, -r_free, L + r_free, W + r_free)
    sys = ParticleSystem(Particle, dom, h)
    generate_particles!(sys, grid, rod, x -> Particle(x=x))
    create_cell_list!(sys)
    force_computation!(sys, 0.)
    return sys
end

function force_computation!(sys::ParticleSystem, t::Float64)
    apply!(sys, find_A!)
    apply!(sys, find_B!)
    apply!(sys, find_f!)
    if t < pull_time
        apply!(sys, pull!)
    end
end


#PHYSICS
#-------------------------------------

function find_A!(p::Particle, q::Particle, r::Float64)
    ker = wendland2(h,r)
    x_pq = p.x - q.x
    X_pq = p.X - q.X
    p.A += -ker*outer(X_pq, x_pq)
    p.H += -ker*outer(x_pq, x_pq)
end

function find_B!(p::Particle)
    Hi = inv(p.H)
    p.A = p.A*Hi
    At = trans(p.A)
    G = At*p.A
    P = c_l^2*(det(p.A)-1.0)
    p.B = m*(P*inv(At) + c_s^2*p.A*dev(G))*Hi
end

function find_f!(p::Particle, q::Particle, r::Float64)
    ker = wendland2(h,r)
    rDker = rDwendland2(h,r)
    x_pq = p.x - q.x
    X_pq = p.X - q.X
    #force
    p.f += -ker*(trans(p.A)*(p.B*x_pq))
    p.f += -ker*(trans(q.A)*(q.B*x_pq))
    #"eta" correction (remove this -> energy will not be conserved!)
    k_pq = +trans(p.B)*(X_pq - p.A*x_pq)
    k_qp = -trans(q.B)*(X_pq - q.A*x_pq)
    p.f += rDker*dot(x_pq, k_pq)*x_pq + ker*k_pq
    p.f -= rDker*dot(x_pq, k_qp)*x_pq + ker*k_qp
    #artificial_viscosity
    p.f += 2*m*vol*rDker*nu*(p.v - q.v)
end

function pull!(p::Particle)
    if p.X[1] > L-h
        p.f += RealVector(0., (vol*pull_force)/(h*W), 0.)
    end
end

function update_v!(p::Particle)
    p.v += 0.5*dt*p.f/m
    #dirichlet bc
    if p.X[1] < h
        p.v = VEC0
    end
end

function update_x!(p::Particle)
    p.x += dt*p.v
    #reset vars
    p.H = MAT0
    p.A = MAT0
    p.f = VEC0
    p.e = 0.
end

function find_e!(p::Particle, q::Particle, r::Float64)
    eta = inv(p.A)*(p.X - q.X) - (p.x - q.x)
    p.e += dot(eta, eta)
end

function particle_energy(p::Particle)
    d = abs(det(p.A))
    G = trans(p.A)*p.A
    G0 = dev(G)
    E_kinet = 0.5*m*dot(p.v, p.v)
    E_shear = 0.25*m*c_s^2*LinearAlgebra.norm(G0,2)^2
    E_press = m*c_l^2*(d - 1.0 - log(d))
    return E_kinet + E_shear + E_press
end 

#TIME ITERATION
#--------------

function main()
    sys = make_geometry()
    out = new_pvd_file("results/rod")
    csv_data = open("results/rod/rod.csv", "w")

    #select top-right corner
    p_sel = argmax(p -> abs(p.x[1]) + abs(p.x[2]), sys.particles) 
    
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
            apply!(sys, find_e!)
            save_frame!(out, sys, :v, :A, :e)
        end
        #verlet scheme:
        apply!(sys, update_v!)
        apply!(sys, update_x!)
        create_cell_list!(sys)
        force_computation!(sys, t)
        apply!(sys, update_v!)
    end
    save_pvd_file(out)
    close(csv_data)
    #plot result
    data = CSV.read("results/rod/rod.csv", DataFrame; header=false)
    p1 = plot(data[:,1], data[:,2], label = "rod-test", xlabel = "time", ylabel = "amplitude")
    p2 = plot(data[:,1], data[:,3], label = "rod-test", xlabel = "time", ylabel = "energy")
    savefig(p1, "results/rod/amplitude.pdf")
    savefig(p2, "results/rod/energy.pdf")
end

end