module twist
using SmoothedParticles
using Parameters
import LinearAlgebra
import StaticArrays
#using Plots
#using CSV
#using DataFrames

#DECLARE CONSTANTS
#-----------------

const H = 6.0
const W = 1.0
const omega = 105.0


const rho0 = 1100.0
const Y = 17e6 #Young modulus
const nu = 0.495 #Poisson ratio

const c_s = sqrt(0.5/rho0*Y/(1.0 + nu))
const c_0 = sqrt(nu*Y/(rho0*(1.0 + nu)*(1.0 - 2*nu)))
const c_p = c_0
const c = sqrt(c_0^2 + 4/3*c_s^2)


const dr = W/24
const h = 3.0dr
const m0 = rho0*dr*dr*dr

const dt = 0.2*dr/c
const t_end = 0.01 #0.2
const dt_plot = max(t_end/200, dt)

function init_velocity(x::RealVector)::RealVector
    return (x[3] > 0.)*omega*sin(0.5*pi*x[3]/H)*RealVector(x[2], -x[1], 0.)
end

#STRUCTURAL KERNELS
#------------------

@fastmath function wendland3h(h::Float64, r::Float64)::Float64
    x = r/h
    return x < 1.0 ? 21.0*(1.0 - x)^3*(14.0*x^2 - 3.0*x - 1.0)/(pi*h^3) : 0.
end

@fastmath function rDwendland3h(h::Float64, r::Float64)::Float64
    x = r/h
    return x < 1.0 ? 210.0*(1.0 - x)^2*(4.0 - 7.0*x)/(pi*h^5) : 0.
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

const inv = StaticArrays.inv
const det = StaticArrays.det

@inbounds function outer(x::RealVector, y::RealVector)::RealMatrix
    return RealMatrix(
        x[1]*y[1], x[2]*y[1], x[3]*y[1], 
        x[1]*y[2], x[2]*y[2], x[3]*y[2],
        x[1]*y[3], x[2]*y[3], x[3]*y[3]
    )
end

@inbounds function dev(G::RealMatrix)::RealMatrix
    tr = 1/3*(G[1,1] + G[2,2] + G[3,3])
    return G - tr*MAT1
end

@inbounds function trace(G::RealMatrix)::Float64
    return G[1,1] + G[2,2] + G[3,3]
end


#CREATE INITIAL STATE
#--------------------

function make_geometry()
    grid = Grid(dr, :bodycentered)
    column = Box(-0.5W, -0.5W, -h, 0.5W, 0.5W, H + 0.1dr)
    dom = column + BoundaryLayer(column, grid, W)
    sys = ParticleSystem(Particle, dom, h)
    generate_particles!(sys, grid, column, x -> Particle(x=x))
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
    if p.x[3] > 0.
        p.v += 0.5*dt*p.f/p.m
    end
end

function update_x!(p::Particle)
    p.x += 0.5*dt*p.v
end

function find_L!(p::Particle, q::Particle, r::Float64)
    ker = q.m/rho0*rDwendland3(h,r)
    x_pq = p.x-q.x
    v_pq = p.v-q.v
    p.T += ker*outer(x_pq, x_pq)
    p.L += ker*outer(v_pq, x_pq)
end

function update_A!(p::Particle)
    p.L = p.L*inv(p.T)
    p.A = p.A*(MAT1 - 0.5*dt*p.L)*inv(MAT1 + 0.5*dt*p.L)
end

function find_J!(p::Particle, q::Particle, r::Float64)
    x_pq = p.x-q.x
    p.T += q.m/rho0*rDwendland3(h,r)*outer(x_pq, x_pq)
    p.J += q.m/rho0*wendland3(h,r)
    p.K += q.m/rho0*wendland3h(h,r)
end

function find_T!(p::Particle)
    F = StaticArrays.inv(p.A)
    B = F*transpose(F)
    detF = 1.0/p.J
    p.P = -rho0*c_0^2*detF^2*(detF - 1.0)
    p.T = -p.P/rho0*MAT1 - c_s^2*(B - MAT1)*inv(p.T)
end

function find_f!(p::Particle, q::Particle, r::Float64)
    ker = q.m/rho0*rDwendland3(h,r)
    kerh = q.m/rho0*rDwendland3h(h,r)
    x_pq = p.x-q.x
    #stress
    p.f += p.m*ker*(p.T*x_pq)
    p.f += p.m*ker*(q.T*x_pq)
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
    F = StaticArrays.inv(p.A)
    B = F*transpose(F)
    detF = StaticArrays.det(F)
    return 0.5*p.m*c_s^2*(trace(B) - 3.0 - 2.0*log(detF))
end

function pE_shear(p::Particle)::Float64
    detF = 1.0/p.J
    return 0.5*p.m*c_0^2*(detF - 1.0)^2
end

function pE_penalty(p::Particle)::Float64
    return 0.5*p.m*c_p^2*p.K^2
end


#TIME ITERATION
#--------------

function find_minimizer(sys::ParticleSystem, f::Function)::Particle
    p = sys.particles[1]
    for q in sys.particles
        if f(q) < f(p)
            p = q
        end
    end
    return p
end

function main()
    sys = make_geometry()
    p_top = find_minimizer(sys, p -> p.x[1]^2 + p.x[2]^2 - p.x[3])
    out = new_pvd_file("results/twist")
    csv_data = open("results/twist/data.csv", "w")
    @time for k = 0 : Int64(round(t_end/dt))
        t = k*dt
        if (k % Int64(round(dt_plot/dt)) == 0)
            println("t = ", t)
            println("N = ", length(sys.particles)) 
            z = p_top.x[3]/H
            E_kinetic = sum(p -> pE_kinetic(p), sys.particles)
            E_bulk = sum(p -> pE_bulk(p), sys.particles)
            E_shear = sum(p -> pE_shear(p), sys.particles)
            E_penalty = sum(p -> pE_penalty(p), sys.particles)
            E_total = E_kinetic + E_bulk + E_shear + E_penalty
            @show E_total
            @show z
            println()
            write(csv_data, vec2string([t, z, E_total, E_kinetic, E_bulk, E_shear, E_penalty]))
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

#=
function plot_result()
    data = CSV.read("twist4/data.csv", DataFrame; header=false)
    t = data[:,1]
    z = data[:,2]
    E0 = data[1,3]
    E_total = data[:,3]./E0
    E_kinetic = data[:,4]./E0
    E_bulk = data[:,5]./E0
    E_bulk = data[:,6]./E0
    E_penalty = data[:,7]./E0
    #data_ref = Matrix(CSV.read("reference.csv", DataFrame; delim = ',', header=false))
    #t_ref = data_ref[:, 1]
    #z_ref = data_ref[:, 2]
    p1 = plot(
        t, z, 
        label = "SPH",
        xlabel = "time [s]", ylabel = "H [1]",
        legend = :topright,
        linewidth = 2,
        color = :orange,
        size = (600,350)
    )
    #=
    scatter!(p1, t_ref, z_ref, label = "REF",
        markershape = :diamond, ms = 4, color = :blue
    )
    =#
    p2 = plot(
        t, [E_total E_kinetic E_bulk E_shear E_penalty],
        label = ["total" "kinetic" "internal-bulk" "internal-shear" "structural penalty"], 
        xlabel = "time [s]", ylabel = "energy [1]",
        linewidth = 2.0, legend = :outertopleft
    )
    p3 = plot(
        t, E_total .- 1.0,
        xlabel = "time [s]", ylabel = "energy error [1]",
        linewidth = 2.0, legend = :none
    )
    savefig(p1, "twist4/amplitude.pdf")
    savefig(p2, "twist4/energy.pdf")
    savefig(p3, "twist4/energy_error.pdf")
    println("DE = ", maximum(E_total) - minimum(E_total))
end
=#

end #module
