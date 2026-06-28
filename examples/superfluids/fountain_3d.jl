module fountain
using SmoothedParticles
using Parameters
using Printf
using Match

#=
Declare constant parameters
=#

const LOAD_LAST_STATE = false
const SAVE_LAST_STATE = true

const N = 10
const d = 1.55e-3
const dr = d/N       # average particle distance
const h = 3.0*dr        # size of kernel support
const rho0 = 145.2352   # fluid density
const m = rho0*dr^3     # particle mass
const c1 = 40.0         # numerical speed of first sound
const c2 = 20.37        # [m/s]
const g = -9.8*VECZ         #gravitational acceleration
const mu = 5e-4       #dynamic viscosity of helium

const Xn = 28.09/rho0
const Xs = 1.0 - Xn

const T = 1.65          # [K]
const S = 335.0        # [J/(kg*K)]
const X_prime = 0.66*(Xn/Xs)*(c2/S)^2
const beta = 5e-2*(d*rho0*S*S/c2)

##temporal parameters
const dt = 0.2*h/c1
const t_end = 0.2
const dt_frame = max(dt, t_end/100)
const t_heater= 0.01

const CELL_R0 = 0.5*d
const CELL_R1 = 4.5*d
const CELL_Z0 = d + 0.5*(CELL_R1 - CELL_R0)
const CELL_Z1 = 6.0*d
const CAP_WIDTH = 0.45e-3/2
const CAP_HEIGHT = 5.0*d

const HE_LEVEL = CELL_Z1
const CRYO_R = 6.0*d
const CRYO_Z0 = 0.0
const CRYO_Z1 = CELL_Z1 + 8e-2
const WALL_THICKNESS = 3.0*dr

##superleak parameters
const SUPERLEAK_R = d
const SUPERLEAK_Z = CELL_Z0
const SUPERLEAK_THICKNESS = 0.5*d
const SUPERLEAK_FRICTION = 1e6 # [1/s]

##heater parameters
const HEATER_R = 0.5*d
const HEATER_Z = 0.5*CELL_Z0 + 0.5*CELL_Z1
const heater_t = 0.02
const w_dot = 0.05 # [W]


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
    p = Rectangle(CELL_R0, CELL_Z0, CELL_R1, CELL_Z1)
    p += Circle(0.5*CELL_R0 + 0.5*CELL_R1, CELL_Z0, 0.5*(CELL_R1-CELL_R0))
    ellipse = Circle(0., 0., CELL_R1-CELL_R0)
    TRANSFORM_MAT = RealMatrix(
        1.0, 0.0, 0.0,
        0.0, 0.2, 0.0, 
        0.0, 0.0, 1.0)
    ellipse = Transform(ellipse, A = TRANSFORM_MAT, b = RealVector(CELL_R0, CELL_Z1, 0.))
    p += Specification(ellipse, x -> x[1] > CELL_R0)
    p += Rectangle(CELL_R0, CELL_Z1, CELL_R0 + CAP_WIDTH, CELL_Z1 + CAP_HEIGHT)
    return p
end

function BodyOfRevolution(s::SmoothedParticles.Shape)::SmoothedParticles.Shape
    box = SmoothedParticles.boundarybox(s)
    z_min = box.x2_min
    z_max = box.x2_max
    R = max(abs(box.x1_min), abs(box.x1_max))
    dom = Box(-R, -R, z_min, R, R, z_max)
    return Specification(dom, x -> is_inside(RealVector(sqrt(x[1]*x[1] + x[2]*x[2]), x[3] ,0.), s))
end

function detect_collisions(_::Particle, _::Particle, r::Float64)
    if r < 0.1*dr
        throw(string("Particles found with rel. distance = ", r/dr))
    end
end

function make_geometry()::ParticleSystem
    grid = Grid(dr, :bodycentered)
    p = Sketch()
    @info "created sketch"
    cell = BodyOfRevolution(p)
    @info "cell generated"
    cryostat = BodyOfRevolution(Rectangle(0., CRYO_Z0, CRYO_R, CRYO_Z1))
    cell_wall = Specification(BoundaryLayer(cryostat - cell, grid, WALL_THICKNESS), x -> (x[1]^2 + x[2]^2 < CRYO_R^2) && (CRYO_Z0 < x[3] < CRYO_Z1))
    cryostat_wall = Specification(BoundaryLayer(cryostat, grid, WALL_THICKNESS), x -> x[3] < CELL_Z1 + CAP_HEIGHT)
    helium = Specification(cryostat - cell, x -> (x[3] < HE_LEVEL))
    sys = ParticleSystem(Particle, cryostat + cryostat_wall, h)
    @info "sys made"
    if LOAD_LAST_STATE
        @info "loading particles from file"
        import_particles!(sys, "results/fountain_3d_state/frame0.vtp", x -> Particle(x=x, type=FLUID))
    else
        generate_particles!(sys, grid, cell_wall, x -> Particle(x=x, type=ADIABATIC))
        generate_particles!(sys, grid, cryostat_wall, x -> Particle(x=x, type=COOLER))
        @info "wall particles generated"
        N_wall = length(sys.particles)
        @show N_wall
        generate_particles!(sys, grid, helium, x -> Particle(x=x, type=FLUID))
        @info "fluid particles generated"
        N_fluid = length(sys.particles) - N_wall
        @show N_fluid
        create_cell_list!(sys)
        @info "cell list created"
        apply!(sys, detect_collisions)
        @info "no collision found"
        apply!(sys, find_c_rho!, self=true)
        for p in sys.particles
            #hydrostatic increment
            p.c_rho += rho0*abs(g[3])/c1^2*(HE_LEVEL - p.x[3])
            p.rho = p.c_rho
        end
        apply!(sys, find_rho!, self=true)
        apply!(sys, find_P!)
        apply!(sys, internal_force!)
    end
    return sys
end

#=
Define particle interactions
=#

function cross_superleak(p::Particle, q::Particle)::Bool
    if (p.x[1]^2 + p.x[2]^2 < SUPERLEAK_R^2)
        if (p.x[3] > SUPERLEAK_Z) && (q.x[3] <= SUPERLEAK_Z)
            return true
        end
        if (q.x[3] > SUPERLEAK_Z) && (p.x[3] <= SUPERLEAK_Z)
            return true
        end
    end
    return false
end

@inbounds function find_c_rho!(p::Particle, q::Particle, r::Float64)
    p.c_rho -= m*wendland3(h,r)
end

function reset_rho!(p::Particle)
    p.rho = p.c_rho
end

@inbounds function find_rho!(p::Particle, q::Particle, r::Float64)
    p.rho += m*wendland3(h,r)
end

@inbounds function balance_of_entropy!(p::Particle, q::Particle, r::Float64)
	ker = m*rDwendland3(h,r)
    x_pq = p.x - q.x
    j_p = p.Xs/p.rho*(S + p.DS)*p.vns
    j_q = p.Xs/q.rho*(S + q.DS)*q.vns
    #entropy production from viscosity
    p.aS += -6.0*ker*mu/(T*p.rho*q.rho)*dot(p.vn - q.vn, x_pq)^2/(r*r + 0.01*h*h) 
    if q.type != ADIABATIC
        p.aS += -ker*dot(j_p + j_q, x_pq)
        #entropy stabilization
        if !cross_superleak(p,q)
            p.aS += ker*beta/(p.rho*q.rho)*(p.DT - q.DT)*2  #*(1.0  + (T + q.DT)/(T + p.DT))
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
    ker = m*rDwendland3(h,r)
    x_pq = p.x - q.x
    #pressure force
    P_force = -ker*(p.P/p.rho^2 + q.P/q.rho^2)*x_pq
    p.a += P_force
    p.as += P_force
    #viscous force
    p.a += 12.0*ker*mu/(p.rho*q.rho)*dot(p.vn - q.vn, x_pq)/(r*r + 0.01*h*h)*x_pq
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
    if p.x[1]^2 + p.x[2]^2 < SUPERLEAK_R^2
        z_ = (p.x[3] - SUPERLEAK_Z)/SUPERLEAK_THICKNESS 
        if abs(z_) < 1.0
            k = 15/16*(1-z_)*(1-z_*z_)*SUPERLEAK_FRICTION
            if abs(z_) < 0.2
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
        r = sqrt(p.x[1]^2 + p.x[2]^2 + (p.x[3] - HEATER_Z)^2)
        dE = 1.0/rho0*wendland3(HEATER_R, r)*w_dot*dt
        p.DS += min(1.0, t/t_heater)*dE/T
    end
end 

function is_leaving(p::Particle)::Bool
    return (p.type == FLUID) && (p.x[1]^2 + p.x[2]^2 <= 0.25*d^2) && (p.x[3] <= CELL_Z1 + CAP_HEIGHT < p.x[3] + p.v[3]*dt_frame)
end

function get_speed(sys::ParticleSystem)::Float64
    v_jet = 0.0
    for p in sys.particles
        if is_leaving(p)
            v_jet += m/(pi*rho0*dt_frame*0.25*d*d)
        end
    end
    return v_jet
end 

function get_temperature(sys::ParticleSystem)::Float64
    x = RealVector(0., 0., HEATER_Z)
    T_cell = SmoothedParticles.sum(
        sys, 
        (p,r) -> m/rho0*p.DT*wendland3(h,r),
        x
    )
    return T_cell
end

#=
Put everything into a time loop
=#
function main()
    @info "making geometry"
	sys = make_geometry()
	out = new_pvd_file("results/fountain_3d")
    csv_data = open("results/fountain_3d/fountain.csv", "w")
    if LOAD_LAST_STATE
        N = 0
        while ispath("results/fountain_3d/frame"*string(N)*".vtp")
            N += 1
        end
        out.frame = N
    end
    write(csv_data, "t,T_cell,v_jet\n")
    v_theory = w_dot/(pi*T*S*rho0*0.25*d*d)
    @show v_theory
	for k = 0 : Int64(round(t_end/dt))
        t = k*dt
		if (k %  Int64(round(dt_frame/dt)) == 0)
			@printf("t = %.6e\n", t)
            println("===========================")
			println("N = ", length(sys.particles))
            T_cell = get_temperature(sys)
            v_jet  = get_speed(sys)
            @show T_cell
            @show v_jet
            write(csv_data, string(t,",",T_cell,",",v_jet,"\n"))
			save_frame!(out, sys, :type, :v, :vs, :DT)
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
    if SAVE_LAST_STATE
        open("results/fountain_3d/result.pvd", "w") do file
            write(file, "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n")
            write(file, "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n")
            write(file, "  <Collection>\n")
            for i in 0:out.frame-1
                write(file, "      <DataSet timestep=\""*string(i)*"\" part=\"0\" file=\"frame"*string(i)*".vtp\"/>\n")
            end
            write(file, "  </Collection>\n")
            write(file, "</VTKFile>\n")
        end
        state_out = new_pvd_file("results/fountain_3d_state")
        save_frame!(state_out, sys, :v, :a, :vs, :as, :rho, :c_rho, :type, :P, :DS, :DT, :Xn, :Xs, :aS, :vns, :vn)
        save_pvd_file(state_out)
    end
end ##function main

end ##module