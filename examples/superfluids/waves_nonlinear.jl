module waves_nonlinear
using SmoothedParticles
using Parameters
using Printf

const L = 1e-2
const T0 = 1.9
const rho0 = 145.4684
const s0 = 725.5
const xn0 = 61.03/rho0
const xs0 = 1.0 - xn0 
const u1 = 40.0
const u2 = 18.83
const t_end = sqrt(2)*L/u1
const dt_frame = t_end/100
const C = T0*s0*s0*xs0/(xn0*u2*u2)
const x_prime = 1.17*T0/C

const WALL = 1.0
const FLUID = 0.0

@inbounds function invert2x2(A::RealMatrix)::RealMatrix
    idet = 1.0/(A[1]*A[5] - A[2]*A[4])
    return RealMatrix(
        idet*A[5], -idet*A[2],  0., 
       -idet*A[4],  idet*A[1],  0.,
               0.,         0.,  1.
    )
end

@inbounds function outer2x2(x::RealVector, y::RealVector)::RealMatrix
    return RealMatrix(
        x[1]*y[1], x[2]*y[1], 0.,
        x[1]*y[2], x[2]*y[2], 0.,
               0.,        0., 1.
    )
end

@with_kw mutable struct Particle <: AbstractParticle
    x::RealVector
    v::RealVector = VEC0
    vs::RealVector = VEC0
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
    err::Float64 = 0.
    s_exact::Float64 = 0.
end

struct GlobalParameters
    dr::Float64
    m::Float64
    h::Float64
    dt::Float64
    A::Float64
    GlobalParameters(N::Int64, A::Float64, dt_hat::Float64) = begin
        dr = L/N
        m = dr*dr*rho0
        h = 6.0*dr
        dt = dt_hat*h/u1
        return new(dr,m,h,dt,A)
    end
end

function s_init(x::RealVector, gp::GlobalParameters)::Float64
    return gp.A*s0*sin(pi*x[1]/L)*sin(pi*x[2]/L)
end

function s_exact(x::RealVector, t::Float64, gp::GlobalParameters)::Float64
    return gp.A*s0*sin(pi*x[1]/L)*sin(pi*x[2]/L)*cos(sqrt(2)*pi*u2*t/L)
end

function make_system(gp::GlobalParameters)::ParticleSystem
    dom = Rectangle(-L/2, -L/2, L/2, L/2)
    grid = Grid(gp.dr, :vogel)
    wall = BoundaryLayer(dom, grid, 2*gp.h)
    sys = ParticleSystem(Particle, dom + wall, gp.h)
    generate_particles!(sys, grid, dom, x -> Particle(x=x, s=s_init(x,gp)))
    #generate_particles!(sys, grid, wall, x -> Particle(x=x, s=s_init(x), type=WALL))
    create_cell_list!(sys)
    #compute_R!(sys, gp)
    apply!(sys, find_P_and_T!)
    apply_with_gp!(sys, find_a_and_as!, gp)
    return sys
end

function move!(p::Particle, gp::GlobalParameters)
    if p.type == FLUID
        p.x += 0.5*gp.dt*p.v
    end
end

function find_P_and_T!(p::Particle)
    vns = v_ns(p)
    p.P = u1*u1*(p.rho - rho0)
    p.T = T0*p.s/C - 0.5*x_prime*dot(vns, vns)
    p.xn = xn0 + p.s*x_prime
    p.xs = xs0 - p.s*x_prime
end

function compute_R!(sys::ParticleSystem, gp::GlobalParameters)
    apply!(sys, reset_R!)
    apply_with_gp!(sys, find_R!, gp)
    apply!(sys, invert_R!)
end

function reset_R!(p::Particle)
    p.R = MAT0
end

function find_R!(p::Particle, q::Particle, r::Float64, gp::GlobalParameters)
    x_pq = p.x - q.x
    p.R += -gp.m/p.rho*rDwendland2(gp.h,r)*outer2x2(x_pq, x_pq)
end

function invert_R!(p::Particle)
    p.R = invert2x2(p.R)
    #p.R = MAT1
end

function v_ns(p::Particle)::RealVector
    return (p.v - p.vs)/p.xn
end

function v_n(p::Particle)::RealVector
    return (p.v - p.xs*p.vs)/p.xn
end

function find_Drho_and_Ds!(p::Particle, q::Particle, r::Float64, gp::GlobalParameters)
    ker = rDwendland2(gp.h,r)
    x_pq = p.x - q.x
    gradw_p = ker*x_pq
    j_p = p.rho*(s0 + p.s)*p.xs*v_ns(p)
    j_q = q.rho*(s0 + q.s)*q.xs*v_ns(q)
    #j_p = p.rho*s0*xs0/xn0*(p.v - p.vs)
    #j_q = q.rho*s0*xs0/xn0*(q.v - q.vs)
    p.Drho += gp.m*dot(gradw_p, p.v - q.v)
    p.Ds   += -gp.m*(dot(j_p, gradw_p)/(p.rho^2) + dot(j_q, gradw_p)/(q.rho^2))
end

function update_rho_and_s!(p::Particle, gp::GlobalParameters)
    p.rho += gp.dt*p.Drho 
    p.s += gp.dt*p.Ds
    p.Drho = 0.0
    p.Ds = 0.0
end

function reset_a_and_as!(p::Particle)
    p.a  = VEC0
    p.as = VEC0
end

function find_a_and_as!(p::Particle, q::Particle, r::Float64, gp::GlobalParameters)
    ker = rDwendland2(gp.h,r)
    x_pq = p.x - q.x
    gradw_p = ker*x_pq
    vns_p = v_ns(p)
    vns_q = v_ns(q)
    pressure = -gp.m*(p.P/p.rho^2*gradw_p + q.P/q.rho^2*gradw_p)
    p.a += pressure
    p.a += - gp.m*p.xn*p.xs/p.rho*dot(gradw_p, vns_p)*vns_p
    p.a += - gp.m*q.xn*q.xs/q.rho*dot(gradw_p, vns_q)*vns_q
    p.as += pressure
    p.as += -gp.m/p.rho*p.xn*dot(v_n(p) - v_n(q), vns_p)*gradw_p
    p.as += -gp.m/p.rho*(s0 + p.s)*(p.T - q.T)*gradw_p
end

function accelerate!(p::Particle, gp::GlobalParameters)
    if p.type == FLUID
        p.v += 0.5*gp.dt*p.a
    end
    p.vs += 0.5*gp.dt*p.as
end

function find_energy(sys::ParticleSystem, gp::GlobalParameters)::NTuple{4,Float64}
    kinetic = 0.
    bulk = 0.
    heat = 0.
    for p in sys.particles
    	vns = v_ns(p)
        kinetic += gp.m*(0.5*dot(p.v,p.v) + 0.5*p.xs*p.xn*dot(vns, vns))
        bulk += gp.m*u1*u1*(log(abs(p.rho/rho0)) + rho0/p.rho - 1.0) 
        heat += gp.m*(T0*p.s + 0.5*T0*p.s*p.s/C)
    end
    total = kinetic + bulk + heat
    return (kinetic, heat, bulk, total)
end

function save_energy(file, sys::ParticleSystem, gp::GlobalParameters, t::Float64, e_char::Float64)
    if t == 0.
        write(file, "t,kinetic,bulk,heat,total\n")
    end
    (kinetic, heat, bulk, total) = find_energy(sys, gp)
    energy_err = (total - e_char)/e_char
    @show energy_err
    write(file, string(t/t_end, ",", kinetic/e_char, ",", bulk/e_char ,",", heat/e_char, ",", total/e_char, "\n"))
end

function save_error(file, sys::ParticleSystem, t::Float64, gp::GlobalParameters)
    if t == 0.
        write(file, "t,error\n")
    end
    l2_error = 0.
    for p in sys.particles
        p.err = (p.s - s_exact(p.x, t, gp))/(gp.A*s0)
        l2_error += p.err*p.err
    end
    l2_error = sqrt(l2_error/length(sys.particles))
    @show l2_error
    write(file, string(t/t_end, ",", l2_error, "\n"))
end

function save_midpoint(file, sys::ParticleSystem, t::Float64, gp::GlobalParameters)
    if t == 0.
        write(file, "t,p,s_computed,s_exact\n")
    end
    x = RealVector(L/4, L/4, 0.)
    p = SmoothedParticles.sum(sys, (p,r) -> gp.m*p.P/rho0*wendland2(gp.h,r), x)
    computed = SmoothedParticles.sum(sys, (p,r) -> gp.m*p.Ds/rho0*wendland2(gp.h,r), x)
    exact = s_exact(x, t, gp)
    write(file, string(t,",",p,",",computed,",",exact,"\n"))
    return
end

function apply_with_gp!(sys::ParticleSystem, fun!::Function, gp::GlobalParameters)
    if hasmethod(fun!, (Particle, Particle, Float64, GlobalParameters))
        apply!(sys, (p::Particle,q::Particle,r::Float64) -> fun!(p,q,r,gp))
    elseif hasmethod(fun!, (Particle, GlobalParameters))
        apply!(sys, (p::Particle) -> fun!(p,gp))
    else
        throw("invalid function type")
    end
end


function main(;N::Int64 = 200, A::Float64 = 0.1, dt_hat = 0.1)
    @show x_prime
    gp = GlobalParameters(N, A, dt_hat)
    sys = make_system(gp)
    name = string("_t", round(Int64, 100*dt_hat),"_A", round(Int64, 100*A))
    #out = new_pvd_file("results/waves_nonlinear")
    energy_data = open("results/waves_nonlinear/energy_data"*name*".csv", "w")
    error_data = open("results/waves_nonlinear/error_data"*name*".csv", "w") 
    midpoint_data = open("results/waves_nonlinear/midpoint_data"*name*".csv", "w")
    e_char = find_energy(sys, gp)[end]
    @time for k = 0 : Int64(round(t_end/gp.dt))
	    if (k %  Int64(round(dt_frame/gp.dt)) == 0)
            t = k*gp.dt
            println("N = ", length(sys.particles))
            @printf("t = %.6e\n", t)
            save_energy(energy_data, sys, gp, t, e_char)
            save_error(error_data, sys, t, gp)
            save_midpoint(midpoint_data, sys, t, gp)
            #save_frame!(out, sys, :v, :vs, :P, :T, :s, :type, :s_exact, :err)
	    end
	    apply_with_gp!(sys, accelerate!, gp)
	    apply_with_gp!(sys, move!, gp)
	    create_cell_list!(sys)
        #compute_R!(sys, gp)
        apply_with_gp!(sys, find_Drho_and_Ds!, gp)
	    apply_with_gp!(sys, update_rho_and_s!, gp)
	    apply_with_gp!(sys, move!, gp)
	    create_cell_list!(sys)
        #compute_R!(sys, gp)
        apply!(sys, find_P_and_T!)
        apply!(sys, reset_a_and_as!)
	    apply_with_gp!(sys, find_a_and_as!, gp)
	    apply_with_gp!(sys, accelerate!, gp)
	end
	#save_pvd_file(out)
    close(energy_data)
end ##function main()


end ##module
