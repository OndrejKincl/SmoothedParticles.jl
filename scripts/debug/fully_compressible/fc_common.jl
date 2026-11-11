# fc_common.jl — shared helpers of the EvolWBFC test scripts (T1–T5)

include(joinpath(@__DIR__, "..", "frozen_background", "common.jl"))
include(srcdir("formulations", "evolutionary", "EvolWBFC.jl"))
const M = EvolWBFC

"the particle system exactly as EvolWBFC.run_sim builds it"
function make_fc_system(global_params::Dict, sim_params::Dict)
	sys = M.make_system(M.Particle, global_params, sim_params; wall_depth=M.inert_wall_depth(global_params, sim_params),
		lid_height=M.lid_height(global_params, sim_params))
	@printf("particle system: %d particles (%d FLUID, %d wall)\n", length(sys.particles),
		count(p -> p.type == M.FLUID, sys.particles), count(p -> p.type == M.MOUNTAIN, sys.particles))
	return sys
end

"the FLUID particle closest to (x, z)"
function nearest_fluid(sys, x::Float64, z::Float64)
	best, best_d = 0, Inf
	for (i, p) in enumerate(sys.particles)
		p.type == M.FLUID || continue
		dd = hypot(p.x[1] - x, p.x[2] - z)
		if dd < best_d
			best, best_d = i, dd
		end
	end
	return best
end

"output directories of the EvolWBFC tests"
fc_datadir() = mkpath(datadir("fully_compressible"))
function save_fc_plot(name::String, fig)
	out_dir = mkpath(plotsdir("fully_compressible"))
	path = joinpath(out_dir, name)
	save(path, fig; px_per_unit = 2)
	println("plot saved: $path")
	return path
end
