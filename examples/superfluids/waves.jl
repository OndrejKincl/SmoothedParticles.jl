module waves
using SmoothedParticles
using Parameters
using Printf

const L = 1e-2
const T0 = 1.9
const rho0 = 145.4684
const s0 = 725.5
const xn0 = 61.03/rho0
const xs0 = 1.0 - xn0 
const u1 = 229.0 #20.0 #229.0
const u2 = 18.83
const Ds_max = 1e-2*s0
const t_end = sqrt(2)*L/u2
const dt_frame = t_end/100

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
    vn::RealVector = VEC0
    vs::RealVector = VEC0
    an::RealVector = VEC0
    as::RealVector = VEC0
    #iH::RealMatrix = MAT0
    P::Float64 = 0.
    Drho::Float64 = 0.
    DT::Float64 = 0.
    Ds::Float64 = 0.
    err::Float64 = 0.
end

struct GlobalParameters
    dr::Float64
    m::Float64
    h::Float64
    dt::Float64
    GlobalParameters(N::Int64) = begin
        dr = L/N
        m = dr*dr*rho0
        h = 2.4*dr
        dt = 0.1*h/u1
        return new(dr,m,h,dt)
    end
end

function Ds_init(x::RealVector)::Float64
    return Ds_max*sin(pi*x[1]/L)*sin(pi*x[2]/L)
end

function Ds_exact(x::RealVector, t::Float64)::Float64
    return Ds_max*sin(pi*x[1]/L)*sin(pi*x[2]/L)*cos(sqrt(2)*pi*u2*t/L)
end

function make_system(gp::GlobalParameters)::ParticleSystem
    dom = Rectangle(-L/2, -L/2, L/2, L/2)
    grid = Grid(gp.dr, :vogel)
    extension = Rectangle(-gp.h-L/2, -gp.h-L/2, L/2+gp.h, L/2+gp.h)
    sys = ParticleSystem(Particle, extension, gp.h)
    generate_particles!(sys, grid, dom, x -> Particle(x=x, Ds=Ds_init(x)))
    create_cell_list!(sys)
    #apply_with_gp!(sys, find_H!, gp)
    #apply!(sys, invert_H!)
    apply!(sys, find_P_and_T!)
    apply_with_gp!(sys, find_an_and_as!, gp)
    apply_with_gp!(sys, accelerate!, gp)
    return sys
end

function move!(p::Particle, gp::GlobalParameters)
    p.x += 0.5*gp.dt*(xn0*p.vn + xs0*p.vs)
end

function find_P_and_T!(p::Particle)
    p.P = u1*u1*p.Drho
    p.DT = xn0/xs0*u2*u2*p.Ds/(s0*s0)
end

function reset_H!(p::Particle)
    p.iH = MAT0
end

function find_H!(p::Particle, q::Particle, r::Float64, gp::GlobalParameters)
    x_pq = p.x - q.x
    p.iH += -gp.m/rho0*rDwendland2(gp.h,r)*outer2x2(x_pq, x_pq)
end

function invert_H!(p::Particle)
    p.iH = invert2x2(p.iH)
    #p.iH = MAT1
end

function update_rho_and_s!(p::Particle, q::Particle, r::Float64, gp::GlobalParameters)
    ker = rDwendland2(gp.h,r)
    x_pq = p.x - q.x
    #gradw_p = ker*p.iH*x_pq
    #gradw_q = ker*q.iH*x_pq
    gradw_p = ker*x_pq
    gradw_q = ker*x_pq
    div_v = -gp.m/rho0*dot(xn0*(p.vn - q.vn) + xs0*(p.vs - q.vs), gradw_p)
    div_vns = gp.m/rho0*(dot(p.vn - p.vs, gradw_p) + dot(q.vn - q.vs, gradw_q))
    p.Drho += -gp.dt*rho0*div_v
    p.Ds   += -gp.dt*xs0*s0*div_vns
end

function reset_an_and_as!(p::Particle)
    p.an = VEC0
    p.as = VEC0
end

function find_an_and_as!(p::Particle, q::Particle, r::Float64, gp::GlobalParameters)
    ker = rDwendland2(gp.h,r)
    x_pq = p.x - q.x
    #gradw_p = ker*p.iH*x_pq
    #gradw_q = ker*q.iH*x_pq
    gradw_p = ker*x_pq
    gradw_q = ker*x_pq
    gradP = gp.m/rho0*(p.P*gradw_p + q.P*gradw_q)
    gradT = -gp.m/rho0*(p.DT - q.DT)*gradw_p
    p.an += -1.0/rho0*gradP - xs0/xn0*s0*gradT
    p.as += -1.0/rho0*gradP + s0*gradT
end

function accelerate!(p::Particle, gp::GlobalParameters)
    p.vn += 0.5*gp.dt*p.an
    p.vs += 0.5*gp.dt*p.as
end

function find_energy(sys::ParticleSystem, gp::GlobalParameters)::NTuple{4,Float64}
    kinetic = 0.
    bulk = 0.
    heat = 0.
    for p in sys.particles
        kinetic += gp.m*(0.5*xn0*dot(p.vn,p.vn) + 0.5*xs0*dot(p.vs,p.vs))
        bulk += 0.5*gp.m*(p.Drho/rho0*u1)^2
        heat += 0.5*gp.m*xn0/xs0*(p.Ds/s0*u2)^2
    end
    total = kinetic + bulk + heat
    return (kinetic, heat, bulk, total)
end

function save_energy(file, sys::ParticleSystem, gp::GlobalParameters, t::Float64, e_char::Float64)
    if t == 0.
        write(file, "t,kinetic,bulk,heat,total\n")
    end
    (kinetic, heat, bulk, total) = find_energy(sys, gp)
    @show total
    write(file, string(t/t_end, ",", kinetic/e_char, ",", bulk/e_char ,",", heat/e_char, ",", total/e_char, "\n"))
end

function save_error(file, sys::ParticleSystem, t::Float64)
    if t == 0.
        write(file, "t,error\n")
    end
    l2_error = 0.
    for p in sys.particles
        p.err = (p.Ds - Ds_exact(p.x, t))/Ds_max
        l2_error += p.err*p.err
    end
    l2_error = sqrt(l2_error/length(sys.particles))
    @show l2_error
    write(file, string(t/t_end, ",", l2_error, "\n"))
end

function save_midpoint(file, sys::ParticleSystem, t::Float64, gp::GlobalParameters)
    if t == 0.
        write(file, "t,Ds_computed,Ds_exact\n")
    end
    x = RealVector(L/4, L/4, 0.)
    computed = SmoothedParticles.sum(sys, (p,r) -> gp.m*p.Ds/rho0*wendland2(gp.h,r), x)
    exact = Ds_exact(x, t)
    write(file, string(t,",",computed,",",exact,"\n"))
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


function main(N::Int64)
    gp = GlobalParameters(N)
	sys = make_system(gp)
	out = new_pvd_file("results/waves")
    energy_data = open("results/waves/energy_data"*string(N)*".csv", "w")
    error_data = open("results/waves/error_data"*string(N)*".csv", "w") 
    midpoint_data = open("results/waves/midpoint_data"*string(N)*".csv", "w") 
    e_char = find_energy(sys, gp)[end]
	@time for k = 0 : Int64(round(t_end/gp.dt))
		if (k %  Int64(round(dt_frame/gp.dt)) == 0)
            t = k*gp.dt
			println("N = ", length(sys.particles))
			@printf("t = %.6e\n", t)
            save_energy(energy_data, sys, gp, t, e_char)
            save_error(error_data, sys, t)
            save_midpoint(midpoint_data, sys, t, gp)
			save_frame!(out, sys, :vn, :vs, :P, :DT, :err, :Ds)
		end
		apply_with_gp!(sys, accelerate!, gp)
		apply_with_gp!(sys, move!, gp)
		create_cell_list!(sys)
        #apply!(sys, reset_H!)
        #apply_with_gp!(sys, find_H!, gp)
        #apply!(sys, invert_H!)
		apply_with_gp!(sys, update_rho_and_s!, gp)
		apply_with_gp!(sys, move!, gp)
		create_cell_list!(sys)
        #apply!(sys, reset_H!)
        #apply_with_gp!(sys, find_H!, gp)
        #apply!(sys, invert_H!)
        apply!(sys, find_P_and_T!)
        apply!(sys, reset_an_and_as!)
		apply_with_gp!(sys, find_an_and_as!, gp)
		apply_with_gp!(sys, accelerate!, gp)
	end
	save_pvd_file(out)
    close(energy_data)
    close(error_data)
    close(midpoint_data)
end ##function main()


end ##module
