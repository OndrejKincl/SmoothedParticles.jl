module ldc
using SmoothedParticles
using Parameters
using Printf

const Re = 400
const u = 20.0
const P0 = 5.0
const t_end = 20.0
const dt_frame = t_end/100
const alpha = 1e-6

const FLUID = 0.0
const WALL = 1.0
const LID = 2.0

@with_kw mutable struct Particle <: AbstractParticle
    x::RealVector
    v::RealVector = VEC0
    a::RealVector = VEC0
    P::Float64 = 0.
    rho::Float64 = 1.0
    Drho::Float64 = 0.
    type::Float64 = FLUID
end

struct SimVars
    dr::Float64
    m::Float64
    h::Float64
    dt::Float64
    N::Int64
    eta::Float64
    SimVars(N::Int64) = begin
        dr = 1.0/N
        m = dr*dr
        h = 3.0*dr
        dt = 0.1*h/u
        eta = 0.1*h
        return new(dr,m,h,dt,N,eta)
    end
end

function apply2!(sys::ParticleSystem, fun!::Function, sv::SimVars)
    apply!(sys, (p::Particle,q::Particle,r::Float64) -> fun!(p,q,r,sv))
end

function apply1!(sys::ParticleSystem, fun!::Function, sv::SimVars)
    apply!(sys, (p::Particle) -> fun!(p,sv))
end

function make_system(sv::SimVars)::ParticleSystem
    grid = Grid(sv.dr, :hexagonal)
    box = Rectangle(0., 0., 1.0, 1.0)
    walls = BoundaryLayer(box, grid, 1.5*sv.h)
    lid   = Specification(walls, x -> x[2] >= 1.0)
    wall = walls - lid
    fluid = box - walls
    sys = ParticleSystem(Particle, walls, sv.h)
    generate_particles!(sys, grid, fluid, x -> Particle(x=x, type=FLUID))
    generate_particles!(sys, grid, lid, x -> Particle(x=x, type=LID))
    generate_particles!(sys, grid, wall, x -> Particle(x=x, type=WALL))
    create_cell_list!(sys)
    apply1!(sys, find_P!, sv)
    apply2!(sys, find_a!, sv)
    apply1!(sys, update_v!, sv)
    return sys
end

function move!(p::Particle, sv::SimVars)
    if p.type == FLUID
        p.x += 0.5*sv.dt*p.v
    end
end

function find_P!(p::Particle, sv::SimVars)
    p.P = u*u*(p.rho-1.0) + P0
end

function find_Drho!(p::Particle, q::Particle, r::Float64, sv::SimVars)
    ker = rDwendland2(sv.h,r)
    x_pq = p.x - q.x
    v_pq = p.v - q.v
    p.Drho += sv.m*ker*dot(v_pq, x_pq)
    if p.type == FLUID && q.type == FLUID
        p.Drho += 2.0*alpha*sv.m*ker*p.rho/q.rho*(p.P - q.P)
    end
end

function update_rho!(p::Particle, sv::SimVars)
    p.rho += sv.dt*p.Drho
    p.Drho = 0.0
end

function reset_a!(p::Particle, sv::SimVars)
    p.a = VEC0
end

function find_a!(p::Particle, q::Particle, r::Float64, sv::SimVars)
    ker = rDwendland2(sv.h,r)
    x_pq = p.x - q.x
    p.a += -sv.m*ker*(p.P/p.rho^2 + q.P/q.rho^2)*x_pq
    v_pq = p.v - q.v
    if q.type == LID
        s = abs(p.x[2] - q.x[2])/(sv.eta + abs(p.x[2] - 1.0))
        v_pq = s*(p.v - VECX)
    end
    p.a += 8/(Re*p.rho*q.rho)*sv.m*ker*dot(v_pq, x_pq)/(r^2 + sv.eta^2)*x_pq
end

function update_v!(p::Particle, sv::SimVars)
    if p.type == FLUID
        p.v += 0.5*sv.dt*p.a
    end
end

function compute_fluxes(sys::ParticleSystem, sv::SimVars, res = 100)
    s = range(0.,1.,length=res)
    fluxes = open(string("results/ldc/Re", Re, "/fluxes", sv.N, ".csv"), "w")
    write(fluxes, "s,v1,v2\n")
    for i in 1:res
		#x-velocity along y-centerline
		x = RealVector(0.5, s[i], 0.)
		gamma = SmoothedParticles.sum(sys, (p,r) -> Float64(p.type==FLUID)*sv.m*wendland2(sv.h,r), x)
        v1 = SmoothedParticles.sum(sys, (p,r) -> Float64(p.type==FLUID)*sv.m*p.v[1]*wendland2(sv.h,r), x)/gamma
		#y-velocity along x-centerline
		x = RealVector(s[i], 0.5, 0.)
		gamma = SmoothedParticles.sum(sys, (p,r) -> Float64(p.type==FLUID)*sv.m*wendland2(sv.h,r), x)
        v2 = SmoothedParticles.sum(sys, (p,r) -> Float64(p.type==FLUID)*sv.m*p.v[2]*wendland2(sv.h,r), x)/gamma
		#save results into csv
        write(fluxes, string(s[i],",",v1,",",v2,"\n"))
    end
    close(fluxes)
end

function main(N::Int64)
    sv = SimVars(N)
	sys = make_system(sv)
	out = new_pvd_file(string("results/ldc/Re",Re))
	for k = 0 : Int64(round(t_end/sv.dt))
		if (k %  Int64(round(dt_frame/sv.dt)) == 0)
            t = k*sv.dt
			println("N = ", length(sys.particles))
			@printf("t = %.6e\n", t)
			save_frame!(out, sys, :v, :P, :type)
		end
		apply1!(sys, update_v!, sv)
		apply1!(sys, move!, sv)
		create_cell_list!(sys)
        apply2!(sys, find_Drho!, sv)
		apply1!(sys, update_rho!, sv)
		apply1!(sys, move!, sv)
		create_cell_list!(sys)
        apply1!(sys, find_P!, sv)
        apply1!(sys, reset_a!, sv)
		apply2!(sys, find_a!, sv)
		apply1!(sys, update_v!, sv)
	end
	save_pvd_file(out)
    compute_fluxes(sys, sv)
end ##function main()


end ##module