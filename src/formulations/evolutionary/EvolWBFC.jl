"""
Fully compressible, well-balanced flow around a mountain with the Witch of Agnesi profile:

h(x) = (h_m a²) / (x² + a²),

The compressible Euler equations with the hydrostatic background (ρ̄, θ̄, P̄)(z) subtracted analytically:

Dρ/Dt = -ρ ∇·v,     Dθ/Dt = 0,     P = (ρ θ / θ_0)^γ,
Dv/Dt = -(1/ρ) ∇δP - g δρ/ρ ê_z,     δP = P - P̄(z),   δρ = ρ - ρ̄(z).

No approximation is involved in the last line (dP̄/dz = -ρ̄ g). Only perturbations are discretised:
- ρ is prognostic (continuity equation with velocity differences), θ is a label, P comes from the full equation of state,
- ∇·v and ∇δP are SPH sums of differences with a renormalised kernel gradient (exact for linear fields),
so every term vanishes identically at rest and for a uniform horizontal flow, on any particle arrangement.

"""

module EvolWBFC

export run_sim

using DrWatson
using Dates
@quickactivate "SPH"

using Printf
using SmoothedParticles
using DataFrames
using Plots
using JLD2
using HDF5
using Glob
using CSV
using Parameters
using LinearAlgebra

# ==============
# INCLUDE CORE SCRIPTS
# ==============

include(srcdir("core", "domains", "evolutionary_domain.jl"))
include(srcdir("core", "diagnostics.jl"))
include(srcdir("core", "time_loop.jl"))
include(srcdir("io", "data_storage.jl"))
include(srcdir("io", "checkpoint.jl"))
include(srcdir("core", "bc", "ebc_shared.jl"))

# ==============
# INCLUDE FORMULATION-SPECIFIC SCRIPTS
# ==============

include(srcdir("core", "bc", "ebc_fc.jl"))

# ==============
# NUMERICAL CONSTANTS
# ==============

const L_TOL = 1e-6              # the moment matrix is inverted only if det > L_TOL tr²
const H_MAX_FACTOR = 1.5        # the smoothing length may exceed its background value at most by this factor
const ZETA_DEFAULT = 1.0        # strength of the acoustic dissipation in the mass balance (sim parameter ζ, 1 = Riemann solver)

# ==============
# ATMOSPHERE (constants of the background state)
# ==============

struct Atmosphere
        g::Float64
        R_mass::Float64
        γ::Float64
        T_bg::Float64
        ρ0::Float64
        K::Float64      # g / (R_mass T_bg), inverse of the density scale height
        θ_0::Float64    # P_0^((γ-1)/γ) / R_mass with P_0 = ρ0 R_mass T_bg
end

function Atmosphere(global_params::Dict)
        @unpack g, R_mass, γ, T_bg, ρ0 = global_params
        P_0 = ρ0 * R_mass * T_bg
        return Atmosphere(g, R_mass, γ, T_bg, ρ0, g / (R_mass * T_bg), (P_0^((γ - 1) / γ)) / R_mass)
end

# ==============
# PARTICLE CONSTRUCTOR
# ==============

mutable struct Particle <: AbstractParticle
        h::Float64              # smoothing length
        x::RealVector           # position
        m::Float64              # mass
        v::RealVector           # velocity
        Dv::RealVector          # acceleration
        ρ::Float64              # total density (prognostic)
        ρ_bg::Float64           # background density
        δρ::Float64             # density perturbation
        Dρ::Float64             # density rate
        P::Float64              # total pressure
        P_bg::Float64           # background pressure
        δP::Float64             # pressure perturbation
        c::Float64              # local speed of sound
        θ::Float64              # total potential temperature (label)
        θ_bg::Float64           # background potential temperature
        δθ::Float64             # potential temperature perturbation
        T::Float64              # total temperature
        L11::Float64            # moment matrix Σ V_b x_ab ⊗ x_ab |w'|/r while summing, its inverse afterwards
        L12::Float64
        L22::Float64
        grad_u::RealVector      # ∇u (x-velocity gradient)
        grad_w::RealVector      # ∇w (y-velocity gradient)
        grad_δP::RealVector     # ∇δP
        w_sum::Float64          # Shepard sums of the wall and open-boundary particles
        v_sum::RealVector
        P_sum::Float64
        θ_sum::Float64
        v_f::RealVector         # averaged velocity of the fluid next to a wall particle (before the reflection)
        type::Float64           # particle type
        id::Int64               # the index of the particle in sys.particles
        row_y::Float64          # flat-grid altitude of the lattice row of the particle
        wrap::Float64           # length after which the row repeats (number of its particles × its lattice spacing)
        x_wrap::Float64         # the particle re-enters at the inflow once it is behind this position

        function Particle(x::RealVector, v::RealVector, type::Float64, global_params::Dict, sim_params::Dict)
                @unpack η, dr = sim_params
                atm = Atmosphere(global_params)

                obj = new(
                        η * dr,         # h
                        x,              # x
                        0.0,            # m
                        v,              # v
                        VEC0,           # Dv
                        0.0,            # ρ
                        0.0,            # ρ_bg
                        0.0,            # δρ
                        0.0,            # Dρ
                        0.0,            # P
                        0.0,            # P_bg
                        0.0,            # δP
                        0.0,            # c
                        0.0,            # θ
                        0.0,            # θ_bg
                        0.0,            # δθ
                        0.0,            # T
                        0.0,            # L11
                        0.0,            # L12
                        0.0,            # L22
                        VEC0,           # grad_u
                        VEC0,           # grad_w
                        VEC0,           # grad_δP
                        0.0,            # w_sum
                        VEC0,           # v_sum
                        0.0,            # P_sum
                        0.0,            # θ_sum
                        VEC0,           # v_f
                        type,           # type
                        0,              # the index of the particle in sys.particles
                        0.0,            # row_y
                        0.0,            # wrap
                        Inf,            # x_wrap
                )

                # initialization: the hydrostatic background at the height of the particle
                obj.m = atm.ρ0 * dr^2
                set_background!(obj, atm)
                obj.ρ = obj.ρ_bg
                obj.θ = obj.θ_bg
                update_pressure!(obj, atm)

                return obj
        end
end

# ==============
# Background values & equation of state
# ==============

function background_density(y::Float64, atm::Atmosphere)
        return atm.ρ0 * exp(-atm.K * y)
end

function background_pot_temperature(y::Float64, atm::Atmosphere)
        return atm.T_bg * exp(((atm.γ - 1.0) / atm.γ) * atm.K * y)
end

function equation_of_state(ρ::Float64, θ::Float64, atm::Atmosphere)
        return (ρ * θ / atm.θ_0)^atm.γ
end

# the background pressure goes through the very same function as the total one: δP = +0.0 whenever ρ = ρ̄ and θ = θ̄
function set_background!(p::Particle, atm::Atmosphere)
        p.ρ_bg = background_density(p.x[2], atm)
        p.θ_bg = background_pot_temperature(p.x[2], atm)
        p.P_bg = equation_of_state(p.ρ_bg, p.θ_bg, atm)
end

function update_pressure!(p::Particle, atm::Atmosphere)
        p.δρ = p.ρ - p.ρ_bg
        p.δθ = p.θ - p.θ_bg
        p.P = equation_of_state(p.ρ, p.θ, atm)
        p.δP = p.P - p.P_bg
        p.c = sqrt(atm.γ * p.P / p.ρ)
        p.T = p.P / (atm.R_mass * p.ρ)
end

# ==============
# Smoothing length
# ==============

# h follows the particle volume, h = η (m/ρ)^(1/2); at rest this is η times the local lattice spacing
function update_smoothing!(p::Particle, η::Float64, h_top::Float64)
        if p.type != MOUNTAIN
                h_new = η * sqrt(p.m / p.ρ)
                h_new = min(h_new, H_MAX_FACTOR * η * sqrt(p.m / p.ρ_bg))
                p.h = min(h_new, h_top)
        end
end

# ==============
# Additional forces: Rayleigh damping
# ==============

function damping_structure(z::Float64, v::RealVector, z_t::Float64, z_β::Float64, γ_r::Float64, v_initial::Float64)
        z_bottom = z_t - z_β
        if z >= z_bottom
                ζ = (z - z_bottom) / z_β

                profile = (sin(π / 2 * ζ))^2

                # friction force: relaxation towards the inflow velocity
                return -γ_r * profile * (v - v_initial * VECX)
        else
                return VEC0
        end
end

# ==============
# Kernel-gradient renormalisation, velocity gradient & mass balance
# ==============

function reset_gradients!(p::Particle)
        p.L11 = 0.0
        p.L12 = 0.0
        p.L22 = 0.0
        p.grad_u = VEC0
        p.grad_w = VEC0
        p.Dρ = 0.0
end

@inbounds function compute_gradients!(p::Particle, q::Particle, r::Float64, ζ::Float64)
        if p.type == FLUID
                x_pq = p.x - q.x
                ker = (q.m / q.ρ) * rDwendland2(0.5 * (p.h + q.h), r)

                # moment matrix A = Σ_b V_b (x_b - x_a) ⊗ ∇w_ab (symmetric, positive definite)
                p.L11 += -ker * x_pq[1] * x_pq[1]
                p.L12 += -ker * x_pq[1] * x_pq[2]
                p.L22 += -ker * x_pq[2] * x_pq[2]

                # Σ_b V_b (f_b - f_a) ∇w_ab: zero for a uniform velocity
                p.grad_u += (q.v[1] - p.v[1]) * ker * x_pq
                p.grad_w += (q.v[2] - p.v[2]) * ker * x_pq

                # acoustic dissipation of the mass balance (the density flux of a Riemann solver, U* - Ū = (P_L - P_R)/(2ρc)),
                # applied to δP only: it vanishes at rest and removes the odd-even pressure modes of the difference operators
                p.Dρ += ζ * (p.δP - q.δP) * r * ker / (0.5 * (p.c + q.c))
        end
end

# L = A⁻¹; with it the SPH gradient of a linear field is exact on any particle arrangement
function renormalise(p::Particle, f::RealVector)
        return RealVector(p.L11 * f[1] + p.L12 * f[2], p.L12 * f[1] + p.L22 * f[2], 0.0)
end

function finalize_gradients!(p::Particle)
        if p.type == FLUID
                det = p.L11 * p.L22 - p.L12 * p.L12
                tr = p.L11 + p.L22
                if det > L_TOL * tr * tr
                        p.L11, p.L12, p.L22 = p.L22 / det, -p.L12 / det, p.L11 / det
                elseif tr > 0.0
                        # degenerate support (e.g. collinear neighbours): isotropic normalisation only
                        p.L11, p.L12, p.L22 = 2.0 / tr, 0.0, 2.0 / tr
                else
                        p.L11, p.L12, p.L22 = 1.0, 0.0, 1.0
                end

                p.grad_u = renormalise(p, p.grad_u)
                p.grad_w = renormalise(p, p.grad_w)

                # balance of mass
                p.Dρ += -p.ρ * (p.grad_u[1] + p.grad_w[2])
        end
end

function compute_density!(p::Particle, dt::Float64, rho_floor::Float64, atm::Atmosphere)
        if p.type == FLUID
                p.ρ = max(p.ρ + dt * p.Dρ, rho_floor)
                set_background!(p, atm)
                update_pressure!(p, atm)
        end
end

# ==============
# Momentum balance
# ==============

function reset_acceleration!(p::Particle)
        p.Dv = VEC0
        p.grad_δP = VEC0
end

@inbounds function balance_of_momentum!(p::Particle, q::Particle, r::Float64, α::Float64, β::Float64, ϵ::Float64, rho_floor::Float64)
        if p.type == FLUID
                x_pq = p.x - q.x
                v_pq = p.v - q.v
                dot_product = SmoothedParticles.dot(x_pq, v_pq)

                h_ij = 0.5 * (p.h + q.h)
                ker = rDwendland2(h_ij, r)

                # Σ_b V_b (δP_b - δP_a) ∇w_ab: zero for a uniform δP, whatever its sign
                p.grad_δP += (q.m / q.ρ) * (q.δP - p.δP) * ker * x_pq

                # artificial viscous force
                if dot_product < 0.0
                        prho = max(p.ρ, rho_floor)
                        qrho = max(q.ρ, rho_floor)
                        c_ij = 0.5 * (p.c + q.c)
                        ρ_ij = 0.5 * (prho + qrho)
                        μ_ij = (h_ij * dot_product) / (r * r + ϵ * h_ij * h_ij)
                        π_ij = (-α * c_ij * μ_ij + β * μ_ij * μ_ij) / ρ_ij

                        # artificial viscous force
                        p.Dv += -q.m * π_ij * ker * x_pq
                end
        end
end

# Dv/Dt = -(1/ρ) ∇δP - g δρ/ρ ê_z (+ AV, which is already accumulated in p.Dv)
function finalize_momentum!(p::Particle, g::Float64)
        if p.type == FLUID
                p.grad_δP = renormalise(p, p.grad_δP)
                p.Dv += -(1.0 / p.ρ) * p.grad_δP - g * (p.δρ / p.ρ) * VECY
        end
end

# ==============
# Move & accelerate
# ==============

function move!(p::Particle, dt::Float64)
        if p.type != MOUNTAIN
                p.x += dt * p.v
        end
end

function accelerate!(p::Particle, dt::Float64, z_t::Float64, z_β::Float64, γ_r::Float64, v_initial::Float64)
        if p.type == FLUID
                p.v += 0.5 * dt * (p.Dv + damping_structure(p.x[2], p.v, z_t, z_β, γ_r, v_initial)) # this is a vector sum
        end
end

# ==============
# State evaluation (shared by the initialisation and the Verlet step)
# ==============

# the boundary particles before the mass balance: wall velocities and the state of the open boundaries (one pass)
function boundary_velocity_pass!(p::Particle, q::Particle, r::Float64, v_initial::Float64)
        compute_wall_velocity!(p, q, r)
        compute_open_boundary!(p, q, r, v_initial)
end

# the boundary particles after the mass balance: wall pressures and the state of the open boundaries (one pass)
function boundary_pressure_pass!(p::Particle, q::Particle, r::Float64, v_initial::Float64, g::Float64)
        compute_wall_pressure!(p, q, r, g)
        compute_open_boundary!(p, q, r, v_initial)
end

# boundary values → velocity gradient → density (advanced by dt) → pressure → boundary pressure → acceleration
function evaluate_state!(sys, dt::Float64, global_params, sim_params)
        # unpack all parameters
        @unpack a, dom_height = global_params
        @unpack h_m, v_initial = sim_params
        @unpack rho_floor, ϵ, α, β = sim_params
        ζ = get(sim_params, :ζ, ZETA_DEFAULT)
        @unpack η, dr = sim_params

        # compute derived parameters
        atm = Atmosphere(global_params)
        g, γ, K = atm.g, atm.γ, atm.K
        h_top = η * dr * exp(K * dom_height / 2.0) * (4 / 3)^(1 / 4)

        # smoothing length of the current densities (kept for all the sums of this evaluation)
        apply!(sys, p -> update_smoothing!(p, η, h_top))

        # velocities of the wall particles, state of the open boundaries
        apply!(sys, p -> reset_boundary_sums!(p))
        apply!(sys, (p, q, r) -> boundary_velocity_pass!(p, q, r, v_initial))
        apply!(sys, p -> finalize_wall_velocity!(p, h_m, a, dom_height))
        apply!(sys, p -> finalize_open_boundary!(p, v_initial, atm))

        # renormalisation matrix, velocity gradient and balance of mass
        apply!(sys, p -> reset_gradients!(p))
        apply!(sys, (p, q, r) -> compute_gradients!(p, q, r, ζ))
        apply!(sys, p -> finalize_gradients!(p))

        # compute density and pressure
        apply!(sys, p -> compute_density!(p, dt, rho_floor, atm))

        # pressure perturbation of the wall particles, state of the open boundaries with the new pressures
        apply!(sys, p -> reset_boundary_sums!(p))
        apply!(sys, (p, q, r) -> boundary_pressure_pass!(p, q, r, v_initial, g))
        apply!(sys, p -> finalize_wall_pressure!(p, γ))
        apply!(sys, p -> finalize_open_boundary!(p, v_initial, atm))

        # compute acceleration (balance of momentum)
        apply!(sys, p -> reset_acceleration!(p))
        apply!(sys, (p, q, r) -> balance_of_momentum!(p, q, r, α, β, ϵ, rho_floor))
        apply!(sys, p -> finalize_momentum!(p, g))
end

"""
    inert_wall_depth(global_params, sim_params)

Depth (on the flat grid) below which a wall particle cannot be reached by the fluid: h ≤ H_MAX_FACTOR η dr near
the ground. The terrain-following map lifts the wall by at most h_m and compresses it by (1 - h_m / dom_height).
"""
function inert_wall_depth(global_params, sim_params)
        @unpack dom_height = global_params
        @unpack η, dr, h_m = sim_params
        return 1.05 * (H_MAX_FACTOR * η * dr + h_m) / (1.0 - h_m / dom_height)
end

"""
    lid_height(global_params, sim_params)

Thickness of the rigid lid above the domain: the kernel support of the uppermost fluid particles (h ≤ h_top).
Without the lid the uppermost particles have a one-sided support and the scheme is unstable there.
"""
function lid_height(global_params, sim_params)
        @unpack g, R_mass, T_bg, dom_height = global_params
        @unpack η, dr = sim_params
        K = g / (R_mass * T_bg)
        return 1.05 * η * dr * exp(K * dom_height / 2.0) * (4 / 3)^(1 / 4)
end

"""
    initialize_state!(sys, global_params, sim_params)

Fresh initialisation of the physical fields at t = 0: the hydrostatic background on every particle.
"""
function initialize_state!(sys, global_params, sim_params)
        @unpack g, R_mass, T_bg, dom_height = global_params
        @unpack η, dr, v_initial = sim_params
        K = g / (R_mass * T_bg)
        h_top = η * dr * exp(K * dom_height / 2.0) * (4 / 3)^(1 / 4)

        # the rows of the lattice flow through the domain as a whole (the fixed inflow ghosts join the moving buffer)
        setup_rows!(sys, global_params, sim_params)
        create_cell_list!(sys)

        for p in sys.particles
                # make_system sets h from the horizontal lattice spacing; here it follows the particle volume, also in the walls
                # (h_top is the edge of the cells of the particle system)
                p.h = min(η * sqrt(p.m / p.ρ_bg), h_top)

                # first guess of the fluid velocity seen by the walls (exact for the uniform inflow)
                if p.type == MOUNTAIN
                        p.v_f = v_initial * VECX
                end
        end

        evaluate_state!(sys, 0.0, global_params, sim_params)
end

# ==============
# Verlet step
# ==============

function verlet_step!(sys, global_params, sim_params)
        # unpack all parameters
        @unpack g, R_mass, γ, T_bg, ρ0, N = global_params
        @unpack dom_height, dom_length, a, z_t = global_params
        @unpack h_m, z_β = sim_params
        @unpack η, dr, dt_rel, t_end, γ_r_rel, v_initial = sim_params

        # compute derived parameters
        atm = Atmosphere(global_params)
        h0 = η * dr
        c = sqrt(65e3 * (γ) / ρ0)
        dt = dt_rel * h0 / c
        γ_r = γ_r_rel * N
        K = atm.K

        # half-step acceleration & drift
        apply!(sys, p -> accelerate!(p, dt, z_t, z_β, γ_r, v_initial))
        apply!(sys, p -> move!(p, dt))

        # lifecycle management: particles behind the outflow buffer re-enter at the inflow, the buffers get the background state
        manage_flow_through!(sys, global_params, sim_params)
        apply!(sys, p -> set_buffer_values!(p, v_initial, atm))

        # create the cell list after particles move and re-enter
        create_cell_list!(sys)

        # density (balance of mass over dt), pressure and acceleration (balance of momentum)
        evaluate_state!(sys, dt, global_params, sim_params)
        apply!(sys, p -> accelerate!(p, dt, z_t, z_β, γ_r, v_initial))
end

# ==============
# Main Entry Point
# ==============

function run_sim(global_params::Dict, sim_params::Dict; restart_dir::String="")
        # ==============
        # Parameters initialization
        # ==============

        # unpack all parameters
        @unpack γ, ρ0 = global_params
        @unpack η, dr, dt_rel, t_end = sim_params

        # compute derived parameters
        h0 = η * dr
        c = sqrt(65e3 * (γ) / ρ0)
        dt = dt_rel * h0 / c
        dt_frame = t_end / 100

        # create the particle system (needed for domain geometry even on restart)
        # the deep part of the wall, which the fluid cannot reach, is not generated
        sys = make_system(Particle, global_params, sim_params; wall_depth=inert_wall_depth(global_params, sim_params), lid_height=lid_height(global_params, sim_params))

        # ==============
        # Restart or fresh initialization
        # ==============

        k_start = 0
        frame_counter = 0
        average_velocities, maximum_velocities, energies = initialize_diagnostics_arrays()

        if !isempty(restart_dir)
                run_dir = validate_restart_dir(restart_dir)
                ckpt = load_checkpoint(run_dir)
                if ckpt === nothing
                        error("Restart requested but no checkpoint found in: $run_dir")
                end
                restore_particles!(sys, ckpt.particles)
                k_start = ckpt.k
                frame_counter = ckpt.frame_counter
                average_velocities = ckpt.average_velocities
                maximum_velocities = ckpt.maximum_velocities
                energies = ckpt.energies
                println("Resumed from checkpoint: t=$(ckpt.t), step=$k_start")
        else
                # ==============
                # Fresh initialization of the physical fields
                # ==============

                initialize_state!(sys, global_params, sim_params)

                # initialize simulation output directory with metadata
                run_dir = initialize_run_directory(sim_params)
        end

        # ==============
        # Time loop
        # ==============

        # use the chosen verlet_step! to advance in time
        step_function!(sys) = verlet_step!(sys, global_params, sim_params)

        # iterate in time
        time_loop(global_params, run_dir, step_function!, t_end, dt, dt_frame, sys,
                  average_velocities, maximum_velocities, energies;
                  k_start=k_start, frame_counter=frame_counter)

        # create a xdmf file to use with the generated h5 files
        generate_sph_xdmf(run_dir)

        # write CSVs once at the end, create diagnostic plots
        finalize_diagnostics(run_dir, average_velocities, maximum_velocities, energies)

        println("\nEnd of the road 🏵.")
        println("Simulation $(basename(run_dir)) complete.")
        println("All artifacts (data, plots, metadata) securely saved to:\n  $run_dir")

        return run_dir
end
end # module
