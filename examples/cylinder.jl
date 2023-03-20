#=

# Flow around cylinder

```@raw html
	<img src='../assets/cylinder.png' width="50%" height="50%" alt='missing' /><br>
```
 
```@raw html
A simulation of flow around cylinder.
All parameters of this benchmark can be found
 <a href="http://www.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark1_re20.html">here.</a>
```

The initial setup was created using the packing technique from a paper *Particle packing algorithm for SPH schemes* by Colagrosi et al.
=#

module cylinder

using Printf
using SmoothedParticles

const folder_name = "results/cylinder"

#=
Declare constants
=#

#geometry parameters
const chan_l = 2.2
const chan_w = 0.41          #width of the channel
const cyl1 = 0.2             #x coordinate of the cylinder
const cyl2 = 0.005           #y coordinate of the cylinder
const cyl_r = 0.05           #radius of the cylinder
const dr = pi*cyl_r/20 		 #average particle distance (decrease to make finer simulation)
const h = 2.4*dr
const bc_width = 6*dr
const x2_min = -chan_w/2 - 6*dr
const x2_max =  chan_w/2 + 6*dr


#physical parameters
const U_max = 0.3       #maximum inflow velocity
const rho0 = 1.0		#referential fluid density
const m0 = rho0*dr^2	#particle mass
const c = 20.0*U_max	#numerical speed of sound
const mu = 1.0e-3		#dynamic viscosity
const nu = 0.1*h*c      #pressure stabilization


#temporal parameters
const dt = 0.1*h/c                     #time step
const t_end = 10.0                     #end of simulation
const dt_frame = max(dt, t_end/200)    #how often data is saved
const t_acc = 1.0                      #time to accelerate to full speed
const t_measure = t_end/2              #time from which we start measuring drag and lift

#particle types
const FLUID = 0.
const INFLOW = 1.
const WALL = 2.
const OBSTACLE = 3.

#=
Declare variables to be stored in a Particle
=#

function check_symmetry(sys::ParticleSystem)::Bool
    for p in sys.particles
        has_mirror = false
        _x = RealVector(p.x[1], -p.x[2], 0.0)
        _v = RealVector(p.v[1], -p.v[2], 0.0)
        _a = RealVector(p.a[1], -p.a[2], 0.0)
        for q in sys.particles
            if (q.x == _x) && (q.v == _v) && (q.a == _a)
                has_mirror = true
            end
        end
        if !has_mirror
            return false
        end
    end
    return true
end

mutable struct Particle <: AbstractParticle
    x::RealVector #position
    v::RealVector #velocity
    a::RealVector #acceleration
    rho::Float64 #density
    Drho::Float64 #rate of density
    P::Float64 #pressure
    m::Float64 #mass
    type::Float64 #particle type
    Particle(x, type=FLUID) = begin
        return new(x, VEC0, VEC0,  rho0, 0., 0., m0, type)
    end
end

function make_system()
    domain = Rectangle(-bc_width, x2_min, chan_l, x2_max)
    sys = ParticleSystem(Particle, domain, h)
    import_particles!(sys, "init/cylinder.vtp", x -> Particle(x))
    @show check_symmetry(sys)
    return sys
end

#Inflow function

function set_inflow_speed!(p::Particle, t::Float64)
    if p.type == INFLOW
        s = min(1.0, t/t_acc)
        v1 = s*U_max*(1.0 - (2.0*p.x[2]/chan_w)^2)
        p.v = v1*VECX
    end
end

#Define interactions between particles


@inbounds function balance_of_mass!(p::Particle, q::Particle, r::Float64)
	ker = q.m*rDwendland2(h,r)
	p.Drho += ker*(dot(p.x-q.x, p.v-q.v))
    if p.type == FLUID && q.type == FLUID
        p.Drho += 2*nu/p.rho*(p.rho - q.rho)
    end
end

function find_pressure!(p::Particle)
    if p.x[1] >= -bc_width + h
	     p.rho += p.Drho*dt
    end
	p.Drho = 0.0
	p.P = c^2*(p.rho - rho0)
end

@inbounds function internal_force!(p::Particle, q::Particle, r::Float64)
	ker = q.m*rDwendland2(h,r)
    x_pq = p.x - q.x
	p.a += -ker*(p.P/p.rho^2 + q.P/q.rho^2)*x_pq
    p.a += 8.0*ker*mu/(p.rho*q.rho)*dot(p.v - q.v, x_pq)/(r*r + 0.01*h*h)*x_pq
end

function move!(p::Particle)
	p.a = VEC0
	if p.type == FLUID || p.type == INFLOW
		p.x += dt*p.v
	end
end

function gravity(p::Particle)
    #f = (RealVector(cyl1, cyl2, 0.0) - p.x)
    f = RealVector(cyl1 - p.x[1], -p.x[2], 0.0)
    absf2 = (cyl1 - p.x[1])^2 + p.x[2]^2
    return 0.3*U_max^2*f/absf2
end

function accelerate!(p::Particle)
	if p.type == FLUID
		p.v += 0.5*dt*(p.a + gravity(p))
	end
end

function add_new_particles!(sys::ParticleSystem)
    new_particles = Particle[]
    for p in sys.particles
        if p.type == INFLOW && p.x[1] >= 0
            p.type = FLUID
            x = p.x - bc_width*VECX
            newp = Particle(x, INFLOW)
            push!(new_particles, newp)
        end
    end
    append!(sys.particles, new_particles)
end

function calculate_force(obstacle::Vector{Particle})::RealVector
    F = sum(p -> p.m*p.a, obstacle)
    L_char = 0.1
    U_mean = 2/3*U_max
    C = 2.0*F/(L_char*U_mean^2)
    return C
end

function  main()
    sys = make_system()
	out = new_pvd_file(folder_name)
    save_frame!(out, sys, :v, :P, :rho, :type)
    C_SPH = VEC0
    C_ref = RealVector(5.57953523384, 0.010618948146, 0.)
    nsteps = Int64(round(t_end/dt))
    nsamples = 0
    obstacle = filter(p -> p.type==OBSTACLE, sys.particles)
    #a modified Verlet scheme
	for k = 1 : nsteps 
        t = k*dt
        apply!(sys, accelerate!)
        apply!(sys, move!)
        add_new_particles!(sys)
        apply!(sys, p -> set_inflow_speed!(p,t))
        create_cell_list!(sys)
		apply!(sys, balance_of_mass!)
        apply!(sys, find_pressure!)
        apply!(sys, internal_force!)
        apply!(sys, accelerate!)
        
        if t > t_measure
            nsamples += 1
            C_SPH += calculate_force(obstacle)
        end

        #save data at selected frames
        
        if (k %  Int64(round(dt_frame/dt)) == 0)
            @show t
            @show check_symmetry(sys)
            println("N = ", length(sys.particles))
            println("C_drag = ", C_SPH[1]/nsamples)
            println("ref value = ", C_ref[1]) 
            println("C_lift = ", C_SPH[2]/nsamples)
            println("ref value = ", C_ref[2]) 
            save_frame!(out, sys, :v, :P, :rho, :type)
        end
	end
	save_pvd_file(out)
    println()
    C_SPH = C_SPH/nsamples
    relative_error = norm(C_SPH - C_ref)/norm(C_ref)
    @show C_SPH
    @show C_ref
    println("relative error = ",100*relative_error,"%")
end

end