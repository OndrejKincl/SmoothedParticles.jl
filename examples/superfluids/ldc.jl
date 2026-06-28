module ldc
using SmoothedParticles
using Parameters
using Printf

const Re = 100.0
const dr = L/N
const T0 = 3.0

const s0 = 1.0
const xn0 = 1.0
const xs0 = 0.0 
const u1 = 20.0
const u2 = 0.0
const h = 2.4*dr
const m = rho0*dr*dr
const dt = 0.1*h/u1
const t_end = 5.0
const dt_frame = max(dt, t_end/100)
const wwall = h

const FLUID = 0.0
const WALL = 1.0
const LID = 2.0

@inbounds function invert2x2(A::RealMatrix, tol = 1e-8)::RealMatrix
    idet = 1.0/(A[1]*A[5] - A[2]*A[4])
    if abs(idet) < tol
        return MAT1
    end
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
    iH::RealMatrix = MAT0
    P::Float64 = 0.
    Drho::Float64 = 0.
    DT::Float64 = 0.
    Ds::Float64 = 0.
    type::Float64 = FLUID
end

function make_system()::ParticleSystem
    grid = Grid(dr, :hexagonal)
    box = Rectangle(0., 0., L, L)
    dom = Rectangle(-0.3L, -0.3L, 1.3L, 1.3L)
    walls = BoundaryLayer(box, grid, wwall)
    lid   = Specification(walls, x -> x[2] >= L)
    wall = walls - lid
    fluid = box - walls
    sys = ParticleSystem(Particle, dom, h)
    generate_particles!(sys, grid, fluid, x -> Particle(x=x, type=FLUID))
    generate_particles!(sys, grid, lid, x -> Particle(x=x, type=LID))
    generate_particles!(sys, grid, wall, x -> Particle(x=x, type=WALL))
    create_cell_list!(sys)
    apply!(sys, find_H!)
    apply!(sys, invert_H!)
    apply!(sys, find_P_and_T!)
    apply!(sys, find_an_and_as!)
    apply!(sys, accelerate!)
    return sys
end

function move!(p::Particle)
    if p.type == FLUID
        p.x += 0.5*dt*(xn0*p.vn + xs0*p.vs)
    end
end

function find_P_and_T!(p::Particle)
    p.P = u1*u1*p.Drho
    p.DT = 0. #xn0/xs0*u2*u2*p.Ds/(s0*s0)
end

function reset_H!(p::Particle)
    p.iH = MAT0
end

function find_H!(p::Particle, q::Particle, r::Float64)
    x_pq = p.x - q.x
    p.iH += -m/rho0*rDwendland2(h,r)*outer2x2(x_pq, x_pq)
end

function invert_H!(p::Particle)
    #p.iH = invert2x2(p.iH)
    p.iH = MAT1
end

function update_rho_and_s!(p::Particle, q::Particle, r::Float64)
    ker = rDwendland2(h,r)
    x_pq = p.x - q.x
    gradw_p = ker*p.iH*x_pq
    gradw_q = ker*q.iH*x_pq
    div_v = -m/rho0*dot(xn0*(p.vn - q.vn) + xs0*(p.vs - q.vs), gradw_p)
    div_vns = m/rho0*(dot(p.vn - p.vs, gradw_p) + dot(q.vn - q.vs, gradw_q))
    p.Drho += -dt*rho0*div_v
    p.Ds   += -dt*xs0*s0*div_vns
end

function reset_an_and_as!(p::Particle)
    p.an = VEC0
    p.as = VEC0
end

function find_an_and_as!(p::Particle, q::Particle, r::Float64)
    ker = rDwendland2(h,r)
    x_pq = p.x - q.x
    gradw_p = ker*p.iH*x_pq
    gradw_q = ker*q.iH*x_pq
    gradP = m/rho0*(p.P*gradw_p + q.P*gradw_q)
    gradT = -m/rho0*(p.DT - q.DT)*gradw_p
    p.an += -1.0/rho0*gradP - xs0/xn0*s0*gradT
    p.as += -1.0/rho0*gradP + s0*gradT
    #viscosity
    if p.type == FLUID && q.type != LID
        p.an += 2*m*ker*nu/rho0*(p.vn - q.vn)
    end
    if p.type == FLUID && q.type == LID
        p.an += 2*m*ker*nu/rho0*(p.vn - vlid*VECX)
    end
end

function accelerate!(p::Particle)
    if p.type == FLUID
        p.vn += 0.5*dt*p.an
        p.vs += 0.5*dt*p.as
    end
end

function find_energy(sys::ParticleSystem)::NTuple{4,Float64}
    kinetic = 0.
    bulk = 0.
    heat = 0.
    for p in sys.particles
        kinetic += m*(0.5*xn0*dot(p.vn,p.vn) + 0.5*xs0*dot(p.vs,p.vs))
        bulk += 0.5*m*(p.Drho/rho0*u1)^2
        #heat += 0.5*m*xn0/xs0*(p.Ds/s0*u2)^2
    end
    total = kinetic + bulk + heat
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

function compute_fluxes(sys::ParticleSystem, res = 100)
    s = range(0.,1.,length=res)
    fluxes = open("results/ldc/Re100/fluxes.csv", "w")
    write(fluxes, "s,v1,v2\n")
    for i in 1:res
		#x-velocity along y-centerline
		x = RealVector(0.5, s[i], 0.)
		gamma = SmoothedParticles.sum(sys, (p,r) -> Float64(p.type==FLUID)*m*wendland2(h,r), x)
        v1 = SmoothedParticles.sum(sys, (p,r) -> Float64(p.type==FLUID)*m*p.vn[1]*wendland2(h,r), x)/gamma
		#y-velocity along x-centerline
		x = RealVector(s[i], 0.5, 0.)
		gamma = SmoothedParticles.sum(sys, (p,r) -> Float64(p.type==FLUID)*m*wendland2(h,r), x)
        v2 = SmoothedParticles.sum(sys, (p,r) -> Float64(p.type==FLUID)*m*p.vn[2]*wendland2(h,r), x)/gamma
		#save results into csv
        write(fluxes, string(s[i],",",v1,",",v2,"\n"))
    end
    close(fluxes)
end

function main()
	sys = make_system()
	out = new_pvd_file("results/ldc/Re100")
    energy_data = open("results/ldc/Re100/energy.csv", "w")
    e_char = find_energy(sys)[end]
	for k = 0 : Int64(round(t_end/dt))
		if (k %  Int64(round(dt_frame/dt)) == 0)
            t = k*dt
			println("N = ", length(sys.particles))
			@printf("t = %.6e\n", t)
            save_energy(energy_data, sys, t, e_char)
			save_frame!(out, sys, :vn, :vs, :P, :DT, :type)
		end
		apply!(sys, accelerate!)
		apply!(sys, move!)
		create_cell_list!(sys)
        apply!(sys, reset_H!)
        apply!(sys, find_H!)
        apply!(sys, invert_H!)
		apply!(sys, update_rho_and_s!)
		apply!(sys, move!)
		create_cell_list!(sys)
        apply!(sys, reset_H!)
        apply!(sys, find_H!)
        apply!(sys, invert_H!)
        apply!(sys, find_P_and_T!)
        apply!(sys, reset_an_and_as!)
		apply!(sys, find_an_and_as!)
		apply!(sys, accelerate!)
	end
	save_pvd_file(out)
    close(energy_data)
    compute_fluxes(sys)
end ##function main()


end ##module