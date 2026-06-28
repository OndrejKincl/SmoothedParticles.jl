module irreversible
using SmoothedParticles
using Parameters
using Printf
using Match
using Random

const N = 30
const L = 1e-2
const dr = L/N      # average particle distance
const h = 3.0*dr        # size of kernel support
const rho0 =  145.2352  # fluid density
const m = rho0*dr^2     # particle mass
const u1 = 25.0         # numerical speed of first sound
const u2 = 20.37        # [m/s]
const g = VEC0 #-9.8*VECY         #gravitational acceleration
const dt = 0.001*h/u1
const T0 = 1.65
const s0 = 335.0
const xn0 = 28.09/rho0
const xs0 = 1.0 - xn0
const C = T0*s0*s0*xs0/(xn0*u2*u2)
const x_prime = 0.66*T0/C
const t_end = 1e-3
const dt_frame = max(dt, t_end/50)

const FLUID = 0.0
const ADIABATIC = 1.0
const COOLER = 2.0

const alpha = 1.0e-3/(rho0*u1*u1)
const beta = 1.0e-2*rho0*s0/T0
const mu =  1.0e-2       #dynamic viscosity of helium


@with_kw mutable struct Particle <: AbstractParticle
    x::RealVector
    v::RealVector = VEC0
    vs::RealVector = VEC0
    vn::RealVector = VEC0
    a::RealVector = VEC0
    as::RealVector = VEC0
    #R::RealMatrix = MAT0
    
    P::Float64 = 0.
    rho::Float64 = rho0
    Ds::Float64 = 0.
    Drho::Float64 = 0.
    
    T::Float64 = 0.
    s::Float64 = 0.
    xn::Float64 = xn0
    xs::Float64 = xs0
    type::Float64 = FLUID
end


function move!(p::Particle)
    p.x += 0.5*dt*p.v
end

function find_P_and_T!(p::Particle)
    vns = (p.v - p.vs)/p.xn
    p.P = u1*u1*(p.rho - rho0)
    p.T = T0*p.s/C - 0.5*x_prime*dot(vns, vns)
    p.xn = xn0 + p.s*x_prime
    p.xs = xs0 - p.s*x_prime
end

#variant with entropy production
function find_Drho_and_Ds!(p::Particle, q::Particle, r::Float64)
    ker = rDwendland2(h,r)
    x_pq = p.x - q.x
    p.Drho +=  2.0*alpha*m*ker*(p.P - q.P)*p.rho/q.rho
    p.Ds   +=  m*ker*beta/(p.rho*q.rho)*(p.T - q.T)*(1.0 + (T0 + q.T)/(T0 + p.T))
    p.Ds   +=  -4.0*m*ker*mu/((T0 + p.T)*p.rho*q.rho)*dot(v_n(p) - v_n(q), x_pq)^2/(r*r + 0.01*h*h) 
    p.Ds   +=  -m*ker*alpha/((T0 + p.T)*p.rho*q.rho)*(p.P - q.P)^2
end

#variant without entropy production 
#=

function find_Drho_and_Ds!(p::Particle, q::Particle, r::Float64)
    ker = rDwendland2(h,r)
    p.Ds   +=   2.0*m*ker*beta/(p.rho*q.rho)*(p.T - q.T)
end

=#

function v_ns(p::Particle)::RealVector
    return (p.v - p.vs)/p.xn
end

function v_n(p::Particle)::RealVector
    return (p.v - p.xs*p.vs)/p.xn
end

function update_rho_and_s!(p::Particle)
    p.rho += dt*p.Drho 	
    if p.type == FLUID
        p.s += dt*p.Ds
    end
    p.Drho = 0.0
    p.Ds = 0.0
end

function reset_a_and_as!(p::Particle)
    p.a  = g
    p.as = g
end

function find_vn!(p::Particle)
    p.vn = v_n(p)    
end

function find_a_and_as!(p::Particle, q::Particle, r::Float64)
    ker = rDwendland2(h,r)
    x_pq = p.x - q.x
    p.a +=   8.0*m*ker*mu/(p.rho*q.rho)*dot(v_n(p)-v_n(q), x_pq)/(r*r + 0.01*h*h)*x_pq
end

function accelerate!(p::Particle)
    p.v += 0.5*dt*p.a
    p.vs += 0.5*dt*p.as
end

function find_energy(sys::ParticleSystem)::NTuple{4,Float64}
    kinetic = 0.
    bulk = 0.
    heat = 0.
    potential = 0.
    for p in sys.particles
    	v_ns = (p.v - p.vs)/p.xn
        kinetic += m*(0.5*dot(p.v,p.v) + 0.5*p.xs*p.xn*dot(v_ns, v_ns))
        bulk += m*u1*u1*( log(abs(p.rho/rho0)) + rho0/p.rho - 1.0) 
        heat += m*(T0*p.s + 0.5*T0*p.s*p.s/C)
        potential += -m*dot(g, p.x)
    end
    total = kinetic + bulk + heat + potential
    return (kinetic, heat, bulk, total)
end

function save_energy(file, sys::ParticleSystem, t::Float64, e_char::Float64)
    if t == 0.
        write(file, "t,kinetic,bulk,heat,total\n")
    end
    (kinetic, heat, bulk, total) = find_energy(sys)
    @show total
    write(file, string(t/t_end, ",", kinetic/e_char, ",", bulk/e_char ,",", heat/e_char, ",", total/e_char, "\n"))
end

function make_geometry()::ParticleSystem
    grid = Grid(dr, :hexagonal)
    box = Rectangle(-L/2, -L/2, L/2, L/2)
    dom = Rectangle(-5*L, -5*L, 5*L, 5*L)
    sys = ParticleSystem(Particle, dom, h)
    generate_particles!(sys, grid, box, x -> Particle(x=x, type=FLUID))
    Random.seed!(135)
    for p in sys.particles
        p.v = RealVector(randn(), randn(), 0.0)
        p.vs = RealVector(randn(), randn(), 0.0)
        p.s = randn()
        p.rho += 5.0 - 2.5*rand()
    end
    create_cell_list!(sys)
    apply!(sys, find_P_and_T!)
    apply!(sys, reset_a_and_as!)
    apply!(sys, find_a_and_as!)
    return sys
end

function vec2string(a::AbstractVector)::String
    out = ""
    for i in 1:length(a)-1
        out = out*string(a[i])*","
    end
    if length(a) > 0
        out = out*string(a[end])
    end
    out = out*"\n"
end

function main()
    sys = make_geometry()
    out = new_pvd_file("results/irreversible")
    energy_data = open("results/irreversible/energy_data.csv", "w")
    e_char = find_energy(sys)[end]
    for k in 0:round(Int, t_end/dt)
        t = k*dt
		if (k %  Int64(round(dt_frame/dt)) == 0)
			@printf("t = %.6e\n", t)
            println("=========================")
			println("N = ", length(sys.particles))
			E = find_energy(sys)
			println("E_err = ", E[end] - e_char)
            apply!(sys, find_vn!)
			save_frame!(out, sys, :P, :type, :v, :vs, :vn, :T, :s, :rho)
            save_energy(energy_data, sys, t, e_char)
            println()
		end
        apply!(sys, accelerate!)
        apply!(sys, move!)
        create_cell_list!(sys)
        apply!(sys, find_Drho_and_Ds!)
        apply!(sys, update_rho_and_s!)
        apply!(sys, move!)
        create_cell_list!(sys)
        apply!(sys, find_P_and_T!)
        apply!(sys, reset_a_and_as!)
        apply!(sys, find_a_and_as!)
        apply!(sys, accelerate!)
    end
    save_pvd_file(out)
    close(energy_data)
end

end