module RaBe
using SmoothedParticles
using Parameters
using CSV
using DataFrames

const save_last_state = true

# Physical parameters
# -------------------

const gamma = 1.4
const c_sound = 20.0
const contrast = 0.1
const Pr = 1.0#0.71
@show Pr
const Ra = 1e6
@show Ra
const Lbox = 1.0
const Hbox = 1.0

# Derived physical parameters
# ---------------------------

const capacity = 1.0/(gamma - 1.0)
const T_top = c_sound^2/gamma
@show T_top
const T_dif = contrast*T_top
const T_bot = T_top + T_dif
@show T_bot
const gravity = T_top/T_dif # = 1.0/contrast
const polytropic_index = gravity/T_dif - 1.0
const viscosity = sqrt(Pr/Ra)
const conductivity = gamma*capacity/sqrt(Pr*Ra)

const temperature_probes = [RealVector((Lbox*1/5,Hbox/2,0.0)), RealVector((Lbox*2/5,Hbox/2,0.0)), RealVector((Lbox*3/5, Hbox/2, 0.0)), RealVector((Lbox*4/5, Hbox/2, 0.0))]

# Numerical parameters
# --------------------

const dr = 5e-3
const h = 6.5*dr
const h_rho = 6.5*dr
const h_T = 6.5*dr
const dt = 0.1*h/c_sound
const t_end_max = 0.1
const t_end_min = 0.0
const t_static_end = 0.1 
const t_relax = 0.01
const dt_static_frame = max(dt, t_static_end/50)
const dt_frame = max(dt, t_end_max/500)
const wall_width = h + dr
const SPH_eps = 0.01*dr
const P_inf = T_top - h*T_dif
const Q_tol = 0.001 # relative tolerance when comparing heat fluxes


# Particle types
# --------------

const FLUID = 0.0
const WALL = 1.0
const HEATER = 2.0
const COOLER = 3.0

# Particle definition
# -------------------

@with_kw mutable struct Particle <: AbstractParticle
    x::RealVector             # position
    m::Float64 = 0.0          # mass
    v::RealVector = VEC0      # velocity
    a::RealVector = VEC0      # acceleration
    rho::Float64 = 0.0        # density
    rho_c::Float64 = 0.0      # density integration constant
    P::Float64 = 0.0          # pressure
    s::Float64 = 0.0          # entropy per mass
    type::Float64 = FLUID     # particle type (fluid or wall)
    T::Float64 = T_top        # temperature
    theta::Float64 = 0.0      # normalized temperature
    Q::RealVector = VEC0      # heat flux
    rho_new::Float64 = 0.0
    QC::Float64 = 0.0 # heat received from the cooler 
    QH::Float64 = 0.0 # heat received from the heater 
end

function particle(x::RealVector)::Particle
    type = FLUID
    if !(0.0 < x[1] < 1.0)
        type = WALL
    end
    if x[2] <= 0.0
        type = HEATER
    end
    if x[2] >= 1.0
        type = COOLER
    end
    p = Particle(x=x, type=type)
    p.T = T_bot - x[2]*T_dif
    p.rho = (p.T/T_top)^polytropic_index
    p.P = p.rho*p.T - P_inf
    #p.T = T_top
    #p.P = gravity*(1.0 + wall_width - x[2])
    #p.rho = (p.P + P_inf)/p.T
    p.s = 1.0/(gamma - 1.0)*log(abs(p.T/T_top)) - log(p.rho)
    p.m = p.rho*dr*dr
    p.rho_c = p.rho
    return p
end

# Dynamics
# --------

function reset_rho!(p::Particle)
    p.rho_new = p.rho_c
    p.Q = VEC0
end

function find_rho!(p::Particle, q::Particle, r::Float64)
    p.rho_new += q.m*wendland2(h_rho,r)
    qT = q.T
    if q.type == HEATER
        qT = T_bot
    end
    if q.type == COOLER
        qT = T_top
    end
    if (p.type == FLUID || p.type == WALL)
        p.Q -= conductivity*p.T*(q.m/q.rho)*(1.0 - p.T/qT)*rDwendland2(h,r)*(p.x - q.x)
    end
end

function find_PT!(p::Particle)
    p.rho = p.rho_new
    p.T = T_top*exp((gamma - 1.0)*p.s)*p.rho^(gamma-1.0)
    p.P = p.rho*p.T - P_inf
    p.theta = (p.T - T_top)/T_dif
    p.a = VEC0
end

function find_a!(p::Particle, q::Particle, r::Float64)
    
    # force computation
    # =================
    ker = rDwendland2(h,r)
    x_pq = p.x - q.x
    p.a -= q.m*rDwendland2(h_rho,r)*(p.P/p.rho^2 + q.P/q.rho^2)*x_pq # pressure
    p.a += 8.0*ker*q.m*viscosity/(p.rho*q.rho)*dot(p.v - q.v, x_pq)/(r^2 + SPH_eps^2)*x_pq # viscosity

    # diffusion step
    # ==============
    qQ = q.Q
    if (p.type == FLUID || p.type == WALL)
        if (q.type == HEATER || q.type == COOLER)
            qQ = p.Q
        end
        ds = dt*ker*q.m/(p.T*p.rho*q.rho)*dot(p.Q + qQ, x_pq) # diffusion
        p.s -= ds
        if q.type == HEATER
            p.QH += -ds * T_bot * p.m
        elseif q.type == COOLER 
            p.QC += -ds * T_top * p.m
        end
        p.s -= 4.0*dt*ker*q.m*viscosity/(p.T*p.rho*q.rho)*(dot(p.v - q.v, x_pq)^2)/(r^2 + SPH_eps^2) # dissipation
    end
end


"""
Measure and reset the heat fluxes from the cooler and heater to the fluid particles. 
Δt is the period during which heat fluxes were accumulated in the fluid particles.
"""
function measure_heat_fluxes(sys::ParticleSystem, Δt::Float64)
    q_COOLER = 0.0
    q_HEATER = 0.0
    for p in sys.particles
        if p.type == FLUID
            q_COOLER += p.QC
            q_HEATER += p.QH
            # reset the cummulative heats 
            p.QC = 0.0 
            p.QH = 0.0
        end
    end
    return q_COOLER/Δt, q_HEATER/Δt
end

function move!(p::Particle)
    if p.type == FLUID
        p.x += dt*p.v
    end
end

function accelerate!(p::Particle)
    if p.type == FLUID
        p.v += 0.5*dt*(p.a - gravity*VECY)
    end
end

# Time Integrator
# ---------------

function update!(sys::ParticleSystem)
    apply!(sys, accelerate!)
    apply!(sys, move!)
    create_cell_list!(sys)
    apply!(sys, reset_rho!)
    apply!(sys, find_rho!, self=true)
    apply!(sys, find_PT!)
    apply!(sys, find_a!)
    apply!(sys, accelerate!)
end

# Define the ParticleSystem
# -------------------------

function makesys(;load_initial_state = false)::ParticleSystem
    grid = Grid(dr, :square)
    box = Rectangle(0.0, 0.0, Lbox, Hbox)
    frame = BoundaryLayer(box, grid, wall_width)
    fluid = Rectangle(0.0, 0.0, Lbox, Hbox)
    sys = ParticleSystem(Particle, frame, h)
    
    if load_initial_state
        import_particles!(sys, "RaBe_statefile/frame0.vtp", x -> Particle(x=x))
    else
        generate_particles!(sys, grid, fluid, particle)
        generate_particles!(sys, grid, frame, particle)
    end
    
    create_cell_list!(sys)
    if !load_initial_state
        apply!(sys, find_rho_c!, self=true)
    end
    apply!(sys, reset_rho!)
    apply!(sys, find_rho!, self=true)
    apply!(sys, find_PT!)
    apply!(sys, find_a!)
    return sys
end

function find_rho_c!(p::Particle, q::Particle, r::Float64)
    p.rho_c -= q.m*wendland2(h_rho,r)
end

# Total energy and entropy 
# ------------------------

function energy(sys::ParticleSystem)::Float64
    out = 0.0
    for p in sys.particles
        pot = p.m*gravity*(p.x[2] - 0.5)
        kin = 0.5*p.m*dot(p.v,p.v)
        int = p.m*capacity*(p.T - T_top) + p.m*P_inf*(1.0/p.rho - 1.0) 
        out += pot + kin + int
    end
    return out
end

function entropy(sys::ParticleSystem)::Float64
    return sum(p -> p.m*p.s, sys.particles)
end

function angular_momentum(sys::ParticleSystem)::Float64
    return sum(p -> -p.m*(p.v[1]*(-(p.x[2]-Hbox))+p.v[2]*(p.x[1]-Lbox)), sys.particles)
end

function moment_of_inertia(sys::ParticleSystem)::Float64
    return sum(p -> p.m*((p.x[2]-Hbox)^2+(p.x[1]-Lbox)^2), sys.particles)
end

function measure_temperature(r::RealVector, sys::ParticleSystem)::Float64
    return sum(p->p.theta*wendland2(h_T, norm(r-p.x))*p.m/p.rho, sys.particles)
end

# The main loop
# -------------

function main(;load_initial_state = false, initial_frame = 0, max_real_hrs= Inf)
    initial_real_time = time() # measures the initial unix time
    @show initial_real_time
    @show max_real_hrs

    is_unstable = (polytropic_index < 1.0/(gamma - 1.0))
    @show polytropic_index
    @show is_unstable
    dist_from_boussinesq = gravity/(contrast*c_sound*c_sound)
    @show dist_from_boussinesq

    # Measure the static heat flux
    folder_name = joinpath("RaBe_Ra"*string(Ra), "static")
    csv_file = joinpath(folder_name, "stat_Q_results.csv")
    Q_static = 0.0
    if load_initial_state
        df = DataFrame(CSV.File(csv_file))
        Qs_COOLER_static = df[!, "Q_COOLER"]
        Qs_HEATER_static = df[!, "Q_HEATER"]
        Q_static = 0.5*(Qs_HEATER_static[end]-Qs_COOLER_static[end])
        @show Q_static
    else # measure the heat flux
        out = new_pvd_file(folder_name)
        sys = makesys()
        create_cell_list!(sys)
        Q_COOLER = 0.0
        Q_HEATER = 0.0
        t_averaging = 0.0
        k = 0
        Qs_COOLER = Float64[]
        Qs_HEATER = Float64[]
        ts = Float64[]

        stat_HEATER_Q_reached = false
        stat_COOLER_Q_reached = false
        HEATER_COOLER_Q_equal = false
        while (!stat_COOLER_Q_reached || !stat_HEATER_Q_reached || !HEATER_COOLER_Q_equal) && (k < t_static_end/dt)
            t = k*dt

            # Measure the heats
            t_averaging += dt
        
            if (k %  Int64(round(dt_static_frame/dt)) == 0) 
                @show t
                N = length(sys.particles)
                @show N

                new_Q_COOLER, new_Q_HEATER = measure_heat_fluxes(sys, t_averaging)
                if abs((Q_COOLER-new_Q_COOLER)/new_Q_COOLER) < Q_tol
                    @info "Stationary COOLER heat flux reached."
                    stat_COOLER_Q_reached = true 
                end
                if abs((Q_HEATER-new_Q_HEATER)/new_Q_HEATER) < Q_tol
                    @info "Stationary HEATER heat flux reached."
                    stat_HEATER_Q_reached = true 
                end
                if abs((new_Q_HEATER-(-new_Q_COOLER))/new_Q_HEATER) < Q_tol
                    @info "HEATER and COOLER heat fluxes equal."
                    HEATER_COOLER_Q_equal = true
                end
                Q_COOLER = new_Q_COOLER
                Q_HEATER = new_Q_HEATER
                t_averaging = 0.0
                @show Q_COOLER
                @show Q_HEATER
                push!(Qs_COOLER, Q_COOLER)
                push!(Qs_HEATER, Q_HEATER)
                push!(ts, t)

                println()
                save_frame!(out, sys, :type, :P, :v, :rho, :theta, :Q)
            end
            apply!(sys, reset_rho!)
            apply!(sys, find_rho!, self=true)
            apply!(sys, find_PT!)
            apply!(sys, find_a!)

            k += 1
        end
        Q_static = 0.5*(Q_HEATER - Q_COOLER) # static heat flux for the Nusselt number
        @show Q_static
        Q_Fourier = conductivity * (T_bot-T_top) * Lbox
        @show Q_Fourier
        df = DataFrame(times = ts, Q_COOLER = Qs_COOLER, Q_HEATER = Qs_HEATER)
        CSV.write(csv_file, df)
        save_pvd_file(out)
    end
    
    # The actual simulation
    @info "Dynamic simulation"
    folder_name = "RaBe_Ra"*string(Ra)
    out = new_pvd_file(folder_name)
    if load_initial_state
        out.frame = initial_frame
    end
    csv_file = joinpath(folder_name, "RaBe_results.csv")

    sys = makesys()
    E0 = energy(sys)
    S0 = entropy(sys)
    @show E0
    @show S0
    ts = Float64[]
    Es = Float64[]
    Ss = Float64[]
    Nus = Float64[]
    Qs_COOLER = Float64[]
    Qs_HEATER = Float64[]
    Q_COOLER = 0.0
    Q_HEATER = 0.0
    t_averaging = 0.0
    temperatures = Array{Float64}[]
    stat_HEATER_Q_reached = false
    stat_COOLER_Q_reached = false
    HEATER_COOLER_Q_equal = false
    k = 0
    terminate = false
    while (!terminate) && (!stat_COOLER_Q_reached || !stat_HEATER_Q_reached || !HEATER_COOLER_Q_equal || (k < t_end_min/dt)) && (k < t_end_max/dt) 
        t = k*dt

        # Measure the heats
        t_averaging += dt

        # Move and update
        if (k %  Int64(round(dt_frame/dt)) == 0)
            @show t
            push!(ts, t)
            dE = energy(sys) - E0
            push!(Es, dE)
            dS = entropy(sys) - S0
            push!(Ss, dS)
            N = length(sys.particles)
            @show dE
            @show dS
            @show N

            new_Q_COOLER, new_Q_HEATER = measure_heat_fluxes(sys, t_averaging)
            if abs((Q_COOLER-new_Q_COOLER)/new_Q_COOLER) < Q_tol
                @info "Stationary COOLER heat flux reached."
                stat_COOLER_Q_reached = true 
            end
            if abs((Q_HEATER-new_Q_HEATER)/new_Q_HEATER) < Q_tol
                @info "Stationary HEATER heat flux reached."
                stat_HEATER_Q_reached = true 
            end
            if abs((new_Q_HEATER-(-new_Q_COOLER))/new_Q_HEATER) < Q_tol
                @info "HEATER and COOLER heat fluxes equal."
                HEATER_COOLER_Q_equal = true
            end
            Q_COOLER = new_Q_COOLER
            Q_HEATER = new_Q_HEATER
            push!(Qs_COOLER, Q_COOLER)
            push!(Qs_HEATER, Q_HEATER)
            t_averaging = 0.0
            @show Q_COOLER
            @show Q_HEATER

            Nu = 0.5*(-Q_COOLER+Q_HEATER)/Q_static
            @show Nu
            push!(Nus, Nu)

            push!(temperatures, [measure_temperature(probe, sys) for probe in temperature_probes])

            println()
            save_frame!(out, sys, :type, :P, :v, :rho, :theta)
        end
        update!(sys)
        k += 1

        time_elapsed = (time() - initial_real_time)/3600 #in hrs
        @show time_elapsed
        if time_elapsed > 0.9 * max_real_hrs
            @info "Ending the simulation so that it does not run out of time, at "*string(time_elapsed)*"hrs"
            terminate = true 
        end
    end
    save_pvd_file(out)

    # Calculating angular velocity
    ω = angular_momentum(sys)/moment_of_inertia(sys)
    @show ω

    # Save the ts, Es, and dS
    Qs = zeros(length(ts))
    for i in 2:length(Qs)
        Qs[i] = (Es[i]-Es[i-1])/(ts[i]-ts[i-1])
    end
    T1_column = [temperatures[i][1] for i in 1:length(temperatures)]
    T2_column = [temperatures[i][2] for i in 1:length(temperatures)]
    T3_column = [temperatures[i][3] for i in 1:length(temperatures)]
    T4_column = [temperatures[i][4] for i in 1:length(temperatures)]
    df = DataFrame(
        times = ts, 
        energy = Es, 
        entropy = Ss, 
        Q_COOLER = Qs_COOLER, 
        Q_HEATER = Qs_HEATER, 
        Nu = Nus, 
        T1 = T1_column, 
        T2 = T2_column, 
        T3 = T3_column, 
        T4 = T4_column)
    if load_initial_state 
        CSV.write(csv_file, df, append=true)
    else
        CSV.write(csv_file, df)
    end

    if save_last_state
        state_file = new_pvd_file("RaBe_statefile")
        save_frame!(state_file, sys, :type, :v, :s, :m, :rho, :rho_c)
        save_pvd_file(state_file)
    end
    return
end


end
