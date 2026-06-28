module fountain
using SmoothedParticles
using Parameters
using Printf
using Match

#=
Declare constant parameters
=#


const N = 8
const r_cap = 1.55e-3/2
const dr = r_cap/N       # average particle distance
const h = 2.8*dr        # size of kernel support
const rho0 = 145.2352   # fluid density
const m = rho0*dr^2     # particle mass
const c1 = 40.0         # numerical speed of first sound
const c2 = 20.37        # [m/s]
const g = -9.8*VECY         #gravitational acceleration
const mu = 1e-4       #dynamic viscosity of helium

const Xn = 28.09/rho0
const Xs = 1.0 - Xn

const T = 1.65          # [K]
const S = 335.0        # [J/(kg*K)]
const X_prime = 0.66*(Xn/Xs)*(c2/S)^2
const beta = 1e-2*(r_cap*rho0*S*S/c2)

##temporal parameters
const dt = 0.2*h/c1
const t_end = 0.5
const dt_frame = max(dt, t_end/200)

const CELL_W = 25e-3
const CELL_Z1 =  120e-3
const CELL_Z0 = -10e-3
const HE_LEVEL = 23e-3
const wwall = 3.0*dr

##superleak parameters
const superleak_w = 15e-3
const superleak_y = 1e-3
const superleak_d = 2e-3
const superleak_k = 1e6 # [1/s]

##heater parameters
const heater_r = 4.0e-3
const heater_y = 8.0e-3
const heater_t = 0.02
const w_dot = 80.0 # [W/m]


##particle types
const FLUID = 0.
const ADIABATIC = 1.
const COOLER = 2.

#=
Declare variables to be stored in a Particle
=#

@with_kw mutable struct Particle <: AbstractParticle
	x::RealVector #position
    v::RealVector = VEC0
    a::RealVector = VEC0 #acceleration

	vs::RealVector = VEC0 #velocity
    as::RealVector = VEC0
	
	rho::Float64 = 0. #density
    c_rho::Float64 = rho0

	type::Float64 #particle type
    P::Float64 = 0.

    DS::Float64 = 0.
    DT::Float64 = 0.
    Xn::Float64 = Xn
    Xs::Float64 = Xs

    aS::Float64 = 0.
    vns::RealVector = VEC0
    vn::RealVector = VEC0
end

const MIRROR_MAT = RealMatrix(
    -1.0, 0.0, 0.0,
     0.0, 1.0, 0.0, 
     0.0, 0.0, 1.0
)

function Corner(cx::Float64, cy::Float64, r::Float64, quartal::Symbol)::SmoothedParticles.Shape
    s = @match quartal begin
        :bl => Rectangle(cx-r,cy-r,cx,cy)
        :br => Rectangle(cx,cy-r,cx+r,cy)
        :tl => Rectangle(cx-r,cy,cx,cy+r)
        :tr => Rectangle(cx,cy,cx+r,cy+r)
    end
    return s - Circle(cx,cy,r)       
end

function Sketch()::SmoothedParticles.Shape
    p = Polygon((0.,0.), (0.,23.7), (17.1,23.7), (17.1,16.7), (9.9,16.7), (9.9,0.))
    p -= Corner(10.7, 13.0, 10.7, :tl)
    p += Corner(11.9, 14.7,  2.0, :tl)
    p -= Corner(15.1, 18.7,  2.0, :br)
    p += Circle(4.95, 0., 4.95)
    #pridej kapilaru
    p += Rectangle(15.6, 23.7, 17.1, 33.7)
    return p
end

function make_geometry()::ParticleSystem
    grid = Grid(dr, :hexagonal)
    p = Sketch()
    left_pipe = Transform(p; A = 1e-3*MAT1, b = -(r_cap + 17.1e-3 + 0.5*dr)*VECX)
    right_pipe = Transform(p; A = 1e-3*MIRROR_MAT, b = (r_cap + 17.1e-3 + 0.5*dr)*VECX)
    pipes = left_pipe + right_pipe
    box = Rectangle(-CELL_W, CELL_Z0, CELL_W, CELL_Z1)
    inner_wall = Specification(BoundaryLayer(box - pipes, grid, wwall), x -> (-CELL_W < x[1] < CELL_W) && (CELL_Z0 < x[2] < CELL_Z1))
    outer_wall = Specification(BoundaryLayer(box, grid, wwall), x -> x[2] < CELL_Z1)
    helium = Specification(box - (pipes + inner_wall), x -> (x[2] < HE_LEVEL))
    sys = ParticleSystem(Particle, outer_wall, h)
    generate_particles!(sys, grid, inner_wall, x -> Particle(x=x, type=ADIABATIC))
    generate_particles!(sys, grid, outer_wall, x -> Particle(x=x, type=COOLER))
    generate_particles!(sys, grid, helium, x -> Particle(x=x, type=FLUID))
    create_cell_list!(sys)
    apply!(sys, find_c_rho!, self=true)
    for p in sys.particles
        #hydrostatic increment
        p.c_rho += rho0*abs(g[2])/c1^2*(HE_LEVEL - p.x[2])
        p.rho = p.c_rho
    end
    apply!(sys, find_rho!, self=true)
    apply!(sys, find_P!)
    apply!(sys, internal_force!)
    return sys
end

#=
Define particle interactions
=#

function cross_superleak(p::Particle, q::Particle)::Bool
    if (abs(p.x[1]) < superleak_w) && (abs(q.x[1]) < superleak_w)
        if (p.x[2] > superleak_y) && (q.x[2] <= superleak_y)
            return true
        end
        if (q.x[2] > superleak_y) && (p.x[2] <= superleak_y)
            return true
        end
    end
    return false
end

@inbounds function find_c_rho!(p::Particle, q::Particle, r::Float64)
    p.c_rho -= m*wendland2(h,r)
end

function reset_rho!(p::Particle)
    p.rho = p.c_rho
end

@inbounds function find_rho!(p::Particle, q::Particle, r::Float64)
    p.rho += m*wendland2(h,r)
end

@inbounds function balance_of_entropy!(p::Particle, q::Particle, r::Float64)
	ker = m*rDwendland2(h,r)
    x_pq = p.x - q.x
    j_p = p.Xs/p.rho*(S + p.DS)*p.vns
    j_q = p.Xs/q.rho*(S + q.DS)*q.vns
    #entropy production from viscosity
    p.aS += -4.0*ker*mu/((T + p.DT)*p.rho*q.rho)*dot(p.vn - q.vn, x_pq)^2/(r*r + 0.01*h*h) 
    if q.type != ADIABATIC
        p.aS += -ker*dot(j_p + j_q, x_pq)
        #entropy stabilization
        if !cross_superleak(p,q)
            p.aS += ker*beta/(p.rho*q.rho)*(p.DT - q.DT)*(1.0 + (T + q.DT)/(T + p.DT))
        end
    end
end

@inbounds function update_S!(p::Particle)
    if p.type == FLUID
        p.DS += dt*p.aS
    end
    p.aS = 0.0
end

@inbounds function find_P!(p::Particle)
    p.P = c1^2*(p.rho - rho0)
    p.Xn = Xn + X_prime*p.DS
    p.Xs = Xs - X_prime*p.DS
    p.DT = (c2/S)^2*p.Xn/p.Xs*p.DS
end

@inbounds function get_normal_velocity!(p::Particle)
    p.vn = (p.v - p.Xs*p.vs)/p.Xn
    p.vns = (p.v - p.vs)/p.Xn
end

@inbounds function internal_force!(p::Particle, q::Particle, r::Float64)
    ker = m*rDwendland2(h,r)
    x_pq = p.x - q.x
    #pressure force
    P_force = -ker*(p.P/p.rho^2 + q.P/q.rho^2)*x_pq
    p.a += P_force
    p.as += P_force
    #viscous force
    p.a += 8.0*ker*mu/(p.rho*q.rho)*dot(p.vn - q.vn, x_pq)/(r*r + 0.01*h*h)*x_pq
    if q.type != ADIABATIC
        #thermomechanical force
        p.as -= (S + p.DS)*ker/p.rho*(p.DT - q.DT)*x_pq
    end
    #convective terms
    p.a += -ker*p.Xn*p.Xs/p.rho*dot(x_pq, p.vns)*p.vns
    p.a += -ker*q.Xn*q.Xs/q.rho*dot(x_pq, q.vns)*q.vns
    p.as += -ker*p.Xn/p.rho*dot(p.vn - q.vn, p.vns)*x_pq
end

function move!(p::Particle)
    if p.type == FLUID
        p.x += 0.5*dt*p.v
    end
    p.a = VEC0
    p.as = VEC0
    p.aS = 0.
end

function superleak!(p::Particle)
    if abs(p.x[1]) < superleak_w
        d = (p.x[2] - superleak_y)/superleak_d 
        if abs(d) < 1.0
            k = 15/16*(1-d*d)*(1-d*d)*superleak_k
            if abs(d) < 0.510
            	k = Inf
            end
            B = 1.0/(1.0 + 0.5*k*dt)
            p.v = B*p.v + (1-B)*Xs*p.vs
        end
    end
end

function accelerate!(p::Particle)
    if p.type == FLUID
	    p.v += 0.5*dt*(p.a + g)
        p.vs += 0.5*dt*(p.as + g)
    end
    if p.type == COOLER
        p.vs += 0.5*dt*(p.as - p.a)
    end
end

function heat_source!(p::Particle, t::Float64)
    if p.type == FLUID
        r = sqrt(p.x[1]^2 + (p.x[2] - heater_y)^2)
        dE = 1.0/rho0*wendland2(heater_r, r)*w_dot*dt
        p.DS += min(1.0, t/heater_t)*dE/T
    end
end 

function is_leaving(p::Particle)::Bool
    return (p.type == FLUID) && (abs(p.x[1]) < r_cap) && (p.x[2] <= 33.7e-3 < p.x[2] + p.v[2]*dt_frame)
end

function get_speed(sys::ParticleSystem)::Float64
    v_jet = 0.0
    for p in sys.particles
        if is_leaving(p)
            v_jet += m/(2.0*rho0*dt_frame*r_cap)
        end
    end
    return v_jet
end 

function get_temperature(sys::ParticleSystem)::Float64
    x = RealVector(0., heater_y, 0.)
    T_cell = SmoothedParticles.sum(
        sys, 
        (p,r) -> m/rho0*p.DT*wendland2(h,r),
        x
    )
    return T_cell
end

#=
Put everything into a time loop
=#
function main()
	sys = make_geometry()
	out = new_pvd_file("results/fountain_mu1")
    csv_data = open("results/fountain_mu1/fountain.csv", "w")
    write(csv_data, "t,T_cell,v_jet\n")
    v_theory = w_dot/(2.0*T*S*rho0*r_cap)
    @show v_theory
    t = 0
	for k = 0 : Int64(round(t_end/dt))
		if (k %  Int64(round(dt_frame/dt)) == 0)
            t = k*dt
			@printf("t = %.6e\n", t)
            println("===========================")
			println("N = ", length(sys.particles))
            T_cell = get_temperature(sys)
            v_jet  = get_speed(sys)
            @show T_cell
            @show v_jet
            write(csv_data, string(t,",",T_cell,",",v_jet,"\n"))
			save_frame!(out, sys, :P, :type, :v, :vs, :vn, :vns, :DT)
            println()
		end
		apply!(sys, accelerate!)
        apply!(sys, get_normal_velocity!)
		apply!(sys, move!)
		create_cell_list!(sys)
        apply!(sys, p -> heat_source!(p,t))
        apply!(sys, superleak!)
		apply!(sys, balance_of_entropy!)
        apply!(sys, update_S!)
		apply!(sys, move!)
		create_cell_list!(sys)
        apply!(sys, reset_rho!)
        apply!(sys, find_rho!, self=true)
        apply!(sys, find_P!)
		apply!(sys, internal_force!)
		apply!(sys, accelerate!)
	end
	save_pvd_file(out)
    close(csv_data)
end ##function main

end ##module