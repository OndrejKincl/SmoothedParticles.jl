using SmoothedParticles
using Parameters

import SmoothedParticles: Grid2, covering, is_inside, boundarybox, RealVector, Shape

# declare particle types
const FLUID = 0.0
const MOUNTAIN = 1.0
const INFLOW_GHOST = 2.0
const INFLOW_INCOMING = 3.0
const OUTFLOW_GHOST = 4.0

"""
    ExpGrid(dr, K)

Creates an isotropically spaced grid whose spacing decays exponentially.
"""
mutable struct ExpGrid <: Grid2
	dr::Float64
	K::Float64
end

function covering(grid::ExpGrid, s::Shape)::Vector{RealVector}
	rect = boundarybox(s)
	dr_base = grid.dr
	K = grid.K

	# hexagonal factors 
	a_factor = (4/3)^(1/4)
	b_factor = (3/4)^(1/4)

	xs = RealVector[]

	# upper half
	y = 0.0
	layer = 0
	while y <= rect.x2_max
		s_local = dr_base * exp(K * y / 2.0)
		a_local = s_local * a_factor
		b_local = s_local * b_factor

		if y >= rect.x2_min
			i_min = Int(floor(rect.x1_min / a_local))
			i_max = Int(ceil(rect.x1_max / a_local))
			offset = iseven(layer) ? 0.0 : 0.5 * a_local

			for i in i_min:i_max
				pos = RealVector((i * a_local) + offset, y, 0.0)
				if is_inside(pos, s)
					push!(xs, pos)
				end
			end
		end

		y += b_local
		layer += 1
	end
	# lower half
	y = 0.0
	layer = 0
	while true
		# step down into the boundary
		layer -= 1

		s_local = dr_base * exp(K * y / 2.0)
		b_local = s_local * b_factor
		a_local = s_local * a_factor

		y -= b_local 

		if y < rect.x2_min
			break
		end

		if y <= rect.x2_max
			i_min = Int(floor(rect.x1_min / a_local))
			i_max = Int(ceil(rect.x1_max / a_local))
			offset = iseven(layer) ? 0.0 : 0.5 * a_local

			for i in i_min:i_max
				pos = RealVector((i * a_local) + offset, y, 0.0)
				if is_inside(pos, s)
					push!(xs, pos)
				end
			end
		end
	end   
	return xs
end

"""
    WitchExpGrid(dr, K, h_m, a, H)

Creates an isotropically spaced grid with the Witch of Agnesi mountain profile. The vertical spacing decays exponentially 
"""
mutable struct WitchExpGrid <: Grid2
	dr::Float64
	K::Float64
	h_m::Float64
	a::Float64
	H::Float64
end

function covering(grid::WitchExpGrid, s::Shape)::Vector{RealVector}
	rect = boundarybox(s) 

	dr_base = grid.dr
	K = grid.K
	h_m = grid.h_m
	a = grid.a
	H = grid.H

	# hexagonal factors 
	a_factor = (4/3)^(1/4)
	b_factor = (3/4)^(1/4)

	xs = RealVector[]

	# upper half
	y = 0.0
	layer = 0
	while y <= rect.x2_max
		s_local = dr_base * exp(K * y / 2.0)
		a_local = s_local * a_factor
		b_local = s_local * b_factor

		if y >= rect.x2_min
			i_min = Int(floor(rect.x1_min / a_local))
			i_max = Int(ceil(rect.x1_max / a_local))
			offset = iseven(layer) ? 0.0 : 0.5 * a_local

			for i in i_min:i_max
				pos_flat = RealVector((i * a_local) + offset, y, 0.0)

				# check if it belongs in the geometric shape
				if is_inside(pos_flat, s)

					# apply the terrain-following transformation
					hx = (h_m * a^2) / (pos_flat[1]^2 + a^2)
					y_warped = pos_flat[2] + hx * ((H - pos_flat[2]) / H)

					pos_warped = RealVector(pos_flat[1], y_warped, 0.0)
					push!(xs, pos_warped)
				end
			end
		end

		y += b_local
		layer += 1
	end

	# lower half
	y = 0.0
	layer = 0
	while true
		# step down into the boundary
		layer -= 1

		s_local = dr_base * exp(K * y / 2.0)
		b_local = s_local * b_factor
		a_local = s_local * a_factor

		y -= b_local 

		if y < rect.x2_min
			break
		end

		if y <= rect.x2_max
			i_min = Int(floor(rect.x1_min / a_local))
			i_max = Int(ceil(rect.x1_max / a_local))
			offset = iseven(layer) ? 0.0 : 0.5 * a_local

			for i in i_min:i_max
				pos_flat = RealVector((i * a_local) + offset, y, 0.0)

				if is_inside(pos_flat, s)
					# Keep mapping the ground so the wall bends with the mountain
					hx = (h_m * a^2) / (pos_flat[1]^2 + a^2)
					y_warped = pos_flat[2] + hx * ((H - pos_flat[2]) / H)

					pos_warped = RealVector(pos_flat[1], y_warped, 0.0)
					push!(xs, pos_warped)
				end
			end
		end
	end   
	return xs
end

"""
make_system(Particle::Type, global_params::Dict, sim_params::Dict; wall_depth=Inf, lid_height=0.0)

Creates the system's geometry and places particles of type `Particle in the positions.
The wall below the ground is generated only down to the (flat-grid) depth `wall_depth`; the default keeps the whole fence.
A rigid lid of wall particles of thickness `lid_height` is placed above the domain if `lid_height > 0`; the default leaves the top open.
"""
function make_system(Particle::Type, global_params::Dict, sim_params::Dict; wall_depth::Float64=Inf, lid_height::Float64=0.0)
	# unpack needed parameters

	@unpack dom_height, dom_length, a = global_params
	@unpack η, dr, dt_rel, t_end, γ_r_rel, v_initial, h_m = sim_params
	@unpack g, R_mass, T_bg = global_params

	K = g / (R_mass * T_bg)

	# derive bc_width from the maximal spacing
	a_factor = (4/3)^(1/4)
	s_max = dr * exp(K * dom_height / 2.0)
	bc_width = 6 * s_max * a_factor 

	# place the particles into WitchExpGrid (terrain following grid)
	grid = WitchExpGrid(dr, K, h_m, a, dom_height)
	domain = Rectangle(- dom_length / 2.0, 0.0, dom_length / 2.0, dom_height)
	fence = BoundaryLayer(domain, grid, bc_width)

	# ==============
	# Extrapolating Boundary Conditions setup
	# ==============

	# horizontal spacing local_a(y) of our grid changes, we need columns
	local_a(y) = dr * exp(K * max(0.0, y) / 2.0) * a_factor	

	## INFLOW ##
	# place 2 layers of incoming particles at the inflow (part of the fence)
	inflow_incoming_spec = Specification(fence, x -> 
				      (-dom_length/2 - 2.1*local_a(x[2]) < x[1] <= -dom_length/2) && 
				      (0 <= x[2] <= dom_height)
				      )

	# place 3 layers of ghost particles at the inflow (part of the fence)
	inflow_ghost_spec = Specification(fence, x -> 
				   (-dom_length/2 - 5.1*local_a(x[2]) < x[1] <= -dom_length/2 - 2.1*local_a(x[2])) &&
				   (0 <= x[2] <= dom_height)
				   )

	## OUTFLOW ##
	# place 3 layers of ghost particles at the outflow (part of the fence)
	outflow_ghost_spec = Specification(fence, x -> 
				    (dom_length/2 <= x[1] < dom_length/2 + 3.1*local_a(x[2])) && 
				    (0 <= x[2] <= dom_height)
				    )
	# specify the top lid
	top_lid_spec = Specification(fence, x -> (x[2] > dom_height))

	## THE REST ##
	# what is not specified is the fence's wall
	fence_wall = fence - inflow_incoming_spec - inflow_ghost_spec - outflow_ghost_spec - top_lid_spec

	# the fence is bc_width deep also below the ground, where the spacing keeps shrinking:
	# most of these particles are out of reach of the fluid and can be left out
	if isfinite(wall_depth)
		fence_wall = Specification(fence_wall, x -> x[2] >= -wall_depth)
	end

	# compute the largest smoothing length
	h_top = η * s_max * a_factor
	
	# initialize the system
	sys = ParticleSystem(Particle, domain + fence, h_top)

	# initial velocity
	initial_velocity = v_initial * VECX

	# domain + walls
	generate_particles!(sys, grid, domain, x -> Particle(x, initial_velocity, FLUID, global_params, sim_params))
	generate_particles!(sys, grid, fence_wall, x -> Particle(x, VEC0, MOUNTAIN, global_params, sim_params))

	# outflow and inflow
	generate_particles!(sys, grid, inflow_incoming_spec, x -> Particle(x, initial_velocity, INFLOW_INCOMING, global_params, sim_params))
	generate_particles!(sys, grid, inflow_ghost_spec, x -> Particle(x, initial_velocity, INFLOW_GHOST, global_params, sim_params))
	generate_particles!(sys, grid, outflow_ghost_spec, x -> Particle(x, initial_velocity, OUTFLOW_GHOST, global_params, sim_params))

	# rigid lid (generated last, so that the order of the other particles does not depend on it)
	if lid_height > 0.0
		lid_spec = Specification(fence, x -> (dom_height < x[2] <= dom_height + lid_height))
		generate_particles!(sys, grid, lid_spec, x -> Particle(x, VEC0, MOUNTAIN, global_params, sim_params))
	end

	# set the smoothing length using the particle's altitude
	for p in sys.particles
		p.h = η * dr * exp(K * max(0.0, p.x[2]) / 2.0) * a_factor
	end

	create_cell_list!(sys)

	# initialize the particle.id
	for i=1:length(sys.particles)
		sys.particles[i].id = i
	end

	return sys
end


