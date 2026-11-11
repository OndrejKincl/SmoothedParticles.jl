using SmoothedParticles
using Parameters

# ==============
# Boundary treatment of the fully compressible perturbation formulation (EvolWBFC)
#
# walls (MOUNTAIN, also the lid):  free-slip mirror of the Shepard-averaged velocity of the neighbouring fluid,
#                                  δP by Shepard extrapolation with the hydrostatic correction of the perturbation
# open boundaries (inflow and      characteristic condition: the acoustic wave leaving the domain is taken from the
# outflow buffers):               neighbouring fluid, the one entering it is zero; what the flow carries (δθ, w) comes
#                                  from the fluid at the outflow and from the background at the inflow
# ==============

"""
    reset_boundary_sums!(p::AbstractParticle)

Resets the Shepard sums of the wall and open-boundary particles.
"""
function reset_boundary_sums!(p::AbstractParticle)
	if p.type != FLUID
		p.w_sum = 0.0
		p.v_sum = VEC0
		p.P_sum = 0.0
		p.θ_sum = 0.0
	end
end

# ==============
# Walls
# ==============

"""
    compute_wall_velocity!(p::AbstractParticle, q::AbstractParticle, r::Float64)

Accumulates the Shepard average of the velocity of the non-wall neighbours of the wall particle p.
The sum is taken relative to the stored averaged velocity `v_f` of the wall particle, so that a uniform flow is reproduced exactly.
"""
function compute_wall_velocity!(p::AbstractParticle, q::AbstractParticle, r::Float64)
	@inbounds if p.type == MOUNTAIN && q.type != MOUNTAIN
		ker = wendland2(0.5 * (p.h + q.h), r)
		p.w_sum += ker
		p.v_sum += ker * (q.v - p.v_f)
	end
end

"""
    finalize_wall_velocity!(p::AbstractParticle, h_m::Float64, a::Float64, dom_height::Float64)

Free-slip wall: the averaged fluid velocity with its normal component reflected, so that the tangential
velocity is continuous across the wall (no spurious divergence when the fluid slides along it) and the
normal velocity vanishes on the surface.
"""
function finalize_wall_velocity!(p::AbstractParticle, h_m::Float64, a::Float64, dom_height::Float64)
	if p.type == MOUNTAIN
		if p.w_sum > 0.0
			p.v_f += p.v_sum / p.w_sum
			# the ground follows the mountain, the lid above the domain is flat
			n = p.x[2] > dom_height ? VECY : mountain_normal(p.x[1], h_m, a)
			p.v = p.v_f - 2.0 * SmoothedParticles.dot(p.v_f, n) * n
		end
	end
end

"""
    compute_wall_pressure!(p::AbstractParticle, q::AbstractParticle, r::Float64, g::Float64)

Accumulates the Shepard extrapolation of δP to the wall particle p. The momentum balance normal to the wall
at rest relative to it reads ∂δP/∂z = -g δρ, which gives the correction term.
"""
function compute_wall_pressure!(p::AbstractParticle, q::AbstractParticle, r::Float64, g::Float64)
	@inbounds if p.type == MOUNTAIN && q.type != MOUNTAIN
		ker = wendland2(0.5 * (p.h + q.h), r)
		p.w_sum += ker
		p.P_sum += ker * (q.δP + g * q.δρ * (q.x[2] - p.x[2]))
	end
end

"""
    finalize_wall_pressure!(p::AbstractParticle, γ::Float64)

Sets δP (and the total pressure) of the wall particles.
"""
function finalize_wall_pressure!(p::AbstractParticle, γ::Float64)
	if p.type == MOUNTAIN
		p.δP = p.w_sum > 0.0 ? p.P_sum / p.w_sum : 0.0
		p.P = p.P_bg + p.δP
		p.c = sqrt(γ * p.P / p.ρ)
	end
end

# ==============
# Open boundaries
# ==============

"""
    compute_open_boundary!(p::AbstractParticle, q::AbstractParticle, r::Float64, v_initial::Float64)

Accumulates, over the FLUID neighbours of the open-boundary particle p, the characteristic variable of the acoustic
wave which leaves the domain through it, J = δP ± ρc (u - v_initial) (+ at the outflow, - at the inflow), and at the
outflow also what the flow carries out (w, δθ).
"""
function compute_open_boundary!(p::AbstractParticle, q::AbstractParticle, r::Float64, v_initial::Float64)
	@inbounds if q.type == FLUID && (p.type == OUTFLOW_GHOST || p.type == INFLOW_INCOMING)
		ker = wendland2(0.5 * (p.h + q.h), r)
		direction = p.type == OUTFLOW_GHOST ? 1.0 : -1.0
		p.w_sum += ker
		p.P_sum += ker * (q.δP + direction * q.ρ * q.c * (q.v[1] - v_initial))
		if p.type == OUTFLOW_GHOST
			p.v_sum += ker * q.v
			p.θ_sum += ker * q.δθ
		end
	end
end

"""
    finalize_open_boundary!(p::AbstractParticle, v_initial::Float64, atm)

State of the open-boundary particles: background + the outgoing acoustic wave (δP = J/2, δu = ±J/(2ρc)), no incoming one.
Without FLUID neighbours the particle keeps the background state moving with the inflow velocity.
"""
function finalize_open_boundary!(p::AbstractParticle, v_initial::Float64, atm)
	if p.type == OUTFLOW_GHOST || p.type == INFLOW_INCOMING
		δP = 0.0
		δθ = 0.0
		u = v_initial
		w = 0.0
		if p.w_sum > 0.0
			J = p.P_sum / p.w_sum
			ρc = p.ρ_bg * sqrt(atm.γ * p.P_bg / p.ρ_bg)
			δP = 0.5 * J
			if p.type == OUTFLOW_GHOST
				u = v_initial + 0.5 * J / ρc
				w = p.v_sum[2] / p.w_sum
				δθ = p.θ_sum / p.w_sum
			else
				u = v_initial - 0.5 * J / ρc
			end
		end

		p.v = RealVector(u, w, 0.0)
		p.θ = p.θ_bg + δθ
		# ρ from P = P̄ + δP and θ; written with ratios, so that ρ = ρ̄ exactly when δP = 0 and δθ = 0
		p.ρ = p.ρ_bg * ((p.P_bg + δP) / p.P_bg)^(1.0 / atm.γ) * (p.θ_bg / p.θ)
		update_pressure!(p, atm)
	end
end

"""
    set_buffer_values!(p::AbstractParticle, v_initial::Float64, atm)

Assigns the background state moving with the inflow velocity to the particles of the inflow and outflow buffers
(the outgoing acoustic wave is added by `finalize_open_boundary!`).
"""
function set_buffer_values!(p::AbstractParticle, v_initial::Float64, atm)
	if p.type == INFLOW_INCOMING || p.type == OUTFLOW_GHOST
		p.v = v_initial * VECX
		p.Dv = VEC0
		p.Dρ = 0.0
		set_background!(p, atm)
		p.ρ = p.ρ_bg
		p.θ = p.θ_bg
		update_pressure!(p, atm)
	end
end

# ==============
# Flow-through of the lattice rows
#
# Every row of the lattice (inflow buffer, fluid, outflow buffer) moves through the domain as a whole. A particle
# which leaves the rear of the outflow buffer re-enters at the rear of the inflow buffer, shifted by a whole number
# of the lattice spacings of its row. In a uniform flow the particles therefore always sit on the moving lattice:
# no seam and no pair of nearly coincident particles appears at the inflow, and no fluid particle runs into a
# fixed ghost at the outflow.
# ==============

"""
    setup_rows!(sys::ParticleSystem, global_params::Dict, sim_params::Dict)

Finds the rows of the lattice among the non-wall particles and stores, in every particle, the flat-grid altitude of
its row (`row_y`), the length after which the row repeats (`wrap`) and the position behind which the particle
re-enters (`x_wrap`). The fixed inflow ghosts become a part of the moving inflow buffer.
"""
function setup_rows!(sys::ParticleSystem, global_params::Dict, sim_params::Dict)
	@unpack g, R_mass, T_bg, dom_height, a = global_params
	@unpack dr, h_m = sim_params
	K = g / (R_mass * T_bg)
	a_factor = (4 / 3)^(1 / 4)

	# altitude on the flat grid: the inverse of the terrain-following map of WitchExpGrid
	function flat_altitude(x::RealVector)
		hx = (h_m * a^2) / (x[1]^2 + a^2)
		return (x[2] - hx) / (1.0 - hx / dom_height)
	end

	idx = [i for (i, p) in enumerate(sys.particles) if p.type != MOUNTAIN]
	y_flat = Dict(i => flat_altitude(sys.particles[i].x) for i in idx)
	sort!(idx; by=i -> y_flat[i])

	n_rows = 0
	irregular = 0
	first = 1
	while first <= length(idx)
		# the rows are at least b_factor·dr apart, round-off of the inverse map is far below 1e-6 dr
		last = first
		while last < length(idx) && y_flat[idx[last+1]] - y_flat[idx[first]] < 1e-6 * dr
			last += 1
		end
		row = idx[first:last]
		row_y = sum(y_flat[i] for i in row) / length(row)
		a_row = dr * exp(K * max(0.0, row_y) / 2.0) * a_factor
		x_min, x_max = extrema(sys.particles[i].x[1] for i in row)

		# the row has to be a gap-free piece of the lattice, otherwise the shift does not close it
		if abs((x_max - x_min) - (length(row) - 1) * a_row) > 1e-6 * a_row
			irregular += 1
		end
		for i in row
			p = sys.particles[i]
			p.row_y = row_y
			p.wrap = length(row) * a_row
			p.x_wrap = x_max + 0.5 * a_row
			if p.type == INFLOW_GHOST
				p.type = INFLOW_INCOMING
			end
		end
		n_rows += 1
		first = last + 1
	end
	if irregular > 0
		@warn "$irregular of $n_rows lattice rows are not gap-free: the re-entering particles will not sit on the lattice"
	end
	return n_rows
end

"""
    manage_flow_through!(sys::ParticleSystem, global_params::Dict, sim_params::Dict)

Moves the particles which left the outflow buffer to the rear of the inflow buffer (onto the terrain-following level
of their row) and sets the type of every non-wall particle by its position.
"""
function manage_flow_through!(sys::ParticleSystem, global_params::Dict, sim_params::Dict)
	@unpack dom_height, dom_length, a = global_params
	@unpack h_m = sim_params

	for p in sys.particles
		if p.type != MOUNTAIN
			if p.x[1] > p.x_wrap
				x = p.x[1] - p.wrap
				hx = (h_m * a^2) / (x^2 + a^2)
				p.x = RealVector(x, p.row_y + hx * ((dom_height - p.row_y) / dom_height), 0.0)
			end

			if p.x[1] > dom_length / 2.0
				p.type = OUTFLOW_GHOST
			elseif p.x[1] > -dom_length / 2.0
				p.type = FLUID
			else
				p.type = INFLOW_INCOMING
			end
		end
	end
end
