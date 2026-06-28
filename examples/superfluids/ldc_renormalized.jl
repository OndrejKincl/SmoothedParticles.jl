module ldc
using SmoothedParticles
using Parameters
using Printf

const L = 1.0
const Re = 400
const path = "results/ldc_renormalized/Re"*string(Re)

const rho0 = 1.0
const vlid = 1.0
const mu = rho0*vlid*L/Re

const u1 = 20.0
const alpha = 0.0#1.0e-4
const t_end = 1.0
const dt_frame = t_end/100

const FLUID = 0.0
const WALL = 1.0
const LID = 2.0

const P0 = 10.0

@inbounds function invert2x2(A::RealMatrix, tol = 0.01)::RealMatrix
    det = (A[1]*A[5] - A[2]*A[4])
	if abs(det) < tol
		return MAT1
	end
	idet = 1.0/det
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
    #vs::RealVector = VEC0
    an::RealVector = VEC0
    #as::RealVector = VEC0
    R::RealMatrix = MAT0
    P::Float64 = 0.
    rho::Float64 = rho0
    Drho::Float64 = 0.
    type::Float64 = FLUID
end

struct GlobalParameters
    dr::Float64
    m::Float64
    h::Float64
    dt::Float64
    N::Int64
    GlobalParameters(N::Int64) = begin
        dr = L/N
        m = dr*dr*rho0
        h = 6.0*dr
        dt = 0.1*h/u1
        return new(dr,m,h,dt,N)
    end
end

function apply_!(sys::ParticleSystem, fun!::Function, gp::GlobalParameters)
    if hasmethod(fun!, (Particle, Particle, Float64, GlobalParameters))
        apply!(sys, (p::Particle,q::Particle,r::Float64) -> fun!(p,q,r,gp))
    elseif hasmethod(fun!, (Particle, GlobalParameters))
        apply!(sys, (p::Particle) -> fun!(p,gp))
    else
        throw("invalid function type")
    end
end

function make_system(gp::GlobalParameters)::ParticleSystem
    grid = Grid(gp.dr, :hexagonal)
    box = Rectangle(0., 0., L, L)
    dom = Rectangle(-0.3L, -0.3L, 1.3L, 1.3L)
    walls = BoundaryLayer(box, grid, 2*gp.h)
    lid   = Specification(walls, x -> x[2] >= L)
    wall = walls - lid
    fluid = box - walls
    sys = ParticleSystem(Particle, dom, gp.h)
    generate_particles!(sys, grid, fluid, x -> Particle(x=x, type=FLUID))
    generate_particles!(sys, grid, lid, x -> Particle(x=x, type=LID))
    generate_particles!(sys, grid, wall, x -> Particle(x=x, type=WALL))
    create_cell_list!(sys)
    compute_R!(sys, gp)
    apply!(sys, find_P!)
    apply_!(sys, find_an!, gp)
    apply_!(sys, accelerate!, gp)
    return sys
end

function move!(p::Particle, gp::GlobalParameters)
    if p.type == FLUID
        p.x += 0.5*gp.dt*p.vn
    end
end

function find_P!(p::Particle)
    p.P = u1*u1*(p.rho - rho0) + P0
end

function find_Drho!(p::Particle, q::Particle, r::Float64, gp::GlobalParameters)
    ker = rDwendland2(gp.h,r)
    x_pq = p.x - q.x
    v_pq = p.vn - q.vn
    P_pq = p.P - q.P
    p.Drho += gp.m*ker*dot(p.R*x_pq, v_pq)
    if p.type == FLUID && q.type == FLUID
        p.Drho += gp.m*ker*2.0*alpha*P_pq*p.rho/q.rho
    end
end

function update_rho!(p::Particle, gp::GlobalParameters)
    p.rho += gp.dt*p.Drho
    p.Drho = 0.
end

function reset_an!(p::Particle)
    p.an = VEC0
end

function find_an!(p::Particle, q::Particle, r::Float64, gp::GlobalParameters)
    ker = rDwendland2(gp.h,r)
    x_pq = p.x - q.x
    v_pq = (q.type != LID) ? (p.vn - q.vn) : (p.vn - vlid*VECX)
    p.an += gp.m/p.rho^2*(p.P - q.P)*ker*x_pq
    p.an += 8.0*gp.m*ker*mu/(p.rho*q.rho)*dot(v_pq, x_pq)/(r^2 + 0.01*gp.h^2)*x_pq
end

function accelerate!(p::Particle, gp::GlobalParameters)
    if p.type == FLUID
        p.vn += 0.5*gp.dt*p.an
    end
end

function compute_R!(sys::ParticleSystem, gp::GlobalParameters)
    apply!(sys, reset_R!)
    apply_!(sys, find_R!, gp)
    apply!(sys, invert_R!)
end

function compute_fluxes(sys::ParticleSystem, gp::GlobalParameters, res = 100)
    s = range(0.,1.,length=res)
    fluxes = open(string(path*"/fluxes", gp.N, ".csv"), "w")
    write(fluxes, "s,v1,v2\n")
    for i in 1:res
		#x-velocity along y-centerline
		x = RealVector(0.5, s[i], 0.)
		gamma = SmoothedParticles.sum(sys, (p,r) -> Float64(p.type==FLUID)*gp.m*wendland2(gp.h,r), x)
        v1 = SmoothedParticles.sum(sys, (p,r) -> Float64(p.type==FLUID)*gp.m*p.vn[1]*wendland2(gp.h,r), x)/gamma
		#y-velocity along x-centerline
		x = RealVector(s[i], 0.5, 0.)
		gamma = SmoothedParticles.sum(sys, (p,r) -> Float64(p.type==FLUID)*gp.m*wendland2(gp.h,r), x)
        v2 = SmoothedParticles.sum(sys, (p,r) -> Float64(p.type==FLUID)*gp.m*p.vn[2]*wendland2(gp.h,r), x)/gamma
		#save results into csv
        write(fluxes, string(s[i],",",v1,",",v2,"\n"))
    end
    close(fluxes)
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
end

function main(N::Int64)
    gp = GlobalParameters(N)
    sys = make_system(gp)
    out = new_pvd_file(path)
    for k = 0 : Int64(round(t_end/gp.dt))
        if (k %  Int64(round(dt_frame/gp.dt)) == 0)
            t = k*gp.dt
            println("N = ", length(sys.particles))
            @printf("t = %.6e\n", t)
            save_frame!(out, sys, :vn, :P, :type)
        end
        apply_!(sys, accelerate!, gp)
        apply_!(sys, move!, gp)
        create_cell_list!(sys)
        compute_R!(sys, gp)
        apply_!(sys, find_Drho!, gp)
        apply_!(sys, update_rho!, gp)
        apply_!(sys, move!, gp)
        create_cell_list!(sys)
        compute_R!(sys, gp)
        apply!(sys, find_P!)
        apply!(sys, reset_an!)
        apply_!(sys, find_an!, gp)
        apply_!(sys, accelerate!, gp)
    end
    save_pvd_file(out)
    compute_fluxes(sys, gp)
end ##function main()


end ##module
