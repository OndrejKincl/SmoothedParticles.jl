# 5: Surface tension simulation in 3D

WARNING: Takes very long time to compute (hours on cluster for 0.1 second)

````julia
module drop

using Printf
include("../src/SPHLib.jl")
using .SPHLib
using Parameters
````

Declare constant parameters

````julia
##geometrical
const dr = 3.7e-5          #average particle distance (decrease to make finer simulation)
const h = 3.0*dr           #size of kernel support
const rad = 1e-3
const deskw= 0.9h

##physical
const rho0 = 1000.   	   #fluid density
const m = rho0*dr^3        #particle mass
const mu = 0.1         #dynamic viscosity of water
const beta = 72e-3      #surface tension
const vol = dr^3
const g = -9.8*VECZ
const c = 10.0*max(sqrt(beta/rho0/dr), sqrt(4*norm(g)*rad))            #numerical speed of sound

##temporal
const dt = 0.3*dr/c
const t_end = 2e-5
#const t_end = 0.02
const dt_frame = max(t_end/50,dt)

#artificial
const s0 = dr*dr/100

const FLUID = 0.
const SOLID = 1.

@with_kw mutable struct Particle <: AbstractParticle
	x::RealVector #position
	v::RealVector = VEC0 #velocity
	a::RealVector = VEC0 #acceleration
	P::Float64 = 0. #pressure
	rho::Float64 = 0. #density
    rho0::Float64 = 0.
    n::RealVector = VEC0 #normal
    type::Float64
end
````

Define geometry and make particles

````julia
function make_system()
	grid = Grid(dr, :cubic)
    drop = Ball(0., 0., rad + h, rad)
    desk = Box(-2rad, -2rad, -deskw, 2rad, 2rad, 0.)
    dom = Box(-2rad, -2rad, -2deskw, 2rad, 2rad, 2.2rad)
    sys = ParticleSystem(Particle, dom, h)
	generate_particles!(sys, grid, drop, x -> Particle(x=x, type=FLUID))
    generate_particles!(sys, grid, desk, x -> Particle(x=x, type=SOLID))
	return sys
end
````

Define particle interactions

````julia
@inbounds function find_n!(p::Particle, q::Particle, r::Float64)
    p.n += 2*vol*vol*rDwendland3(h,r)*(p.x - q.x)
end

function reset_n!(p::Particle)
    p.n = VEC0
end

function normalize_n!(p::Particle)
    s = norm(p.n)
    p.n /= (s + s0)
end

@inbounds function find_rho!(p::Particle, q::Particle, r::Float64)
    p.rho += m*wendland3(h,r)
end

@inbounds function find_rho0!(p::Particle, q::Particle, r::Float64)
    p.rho0 += m*wendland3(h,r)
end

function find_pressure!(p::Particle)
	p.P = c^2*(p.rho - p.rho0)
end

@inbounds function internal_force!(p::Particle, q::Particle, r::Float64)
		ker = m*rDwendland3(h,r)
        #pressure
		p.a += -ker*(p.P/rho0^2 + q.P/rho0^2)*(p.x - q.x)
		#viscosity
        p.a += 2*ker*mu/rho0^2*(p.v - q.v)
        #surface tension
        p.a -= 2*beta/rho0^2*(
            (m*DDwendland3(h,r)-ker)*dot(p.x-q.x, p.n-q.n)*(p.x-q.x)/(r^2 + s0)
            +ker*(p.n-q.n)
        )
end

function reset_a!(p::Particle)
    p.a = VEC0
end

function reset_rho!(p::Particle)
    p.rho = 0.0
end

function move!(p::Particle)
	p.x += (p.type==FLUID)*dt*p.v
end

function accelerate!(p::Particle)
	p.v += (p.type==FLUID)*0.5*dt*(p.a+g)
end

function energy(p::Particle)::Float64
	kinetic = 0.5*m*dot(p.v, p.v)
	internal =  0.5*m*c^2*(p.rho - p.rho0)^2/rho0^2
    s = norm(p.n)
    tensile = beta*(s - s0*log(s/s0+1))
    potential = m*dot(p.x, -g)
    return kinetic + internal + tensile + potential
end
````

Put everything into a time loop

````julia
function verlet_step!(sys)
    apply!(sys, accelerate!)
    apply!(sys, move!)
    create_cell_list!(sys)
    apply!(sys, reset_rho!)
    apply!(sys, find_rho!, self = true)
    apply!(sys, reset_n!)
    apply!(sys, find_n!, self = true)
    apply!(sys, normalize_n!)
    apply!(sys, find_pressure!)
    apply!(sys, reset_a!)
    apply!(sys, internal_force!)
    apply!(sys, accelerate!)
end

function save_results!(out, sys, k)
    @printf("t = %.6e\n", k*dt)
    apply!(sys, reset_n!)
    apply!(sys, find_n!, self = true)
    save_frame!(out, sys, :v, :a, :P, :rho, :rho0, :type, :n)
end


function main()
    E0 = 0.0
	sys = make_system()
	out = new_pvd_file("results/drop")
    #initialization
    create_cell_list!(sys)
    apply!(sys, find_rho0!, self = true)
    apply!(sys, find_rho!, self = true)
    apply!(sys, find_pressure!)
    apply!(sys, find_n!)
    apply!(sys, normalize_n!)
    apply!(sys, internal_force!)
	for k in 0 : Int64(round(t_end/dt))
        verlet_step!(sys)
        if (k %  Int64(round(dt_frame/dt)) == 0)
            save_results!(out, sys, k)
            E = sum(energy, sys.particles)
            if k == 0
                E0 = E
            end
            E = E - E0
            @show E
            println("# of part. = ", length(sys.particles))
            println()
        end
	end
	save_pvd_file(out)
end ## function main


end ## module
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

