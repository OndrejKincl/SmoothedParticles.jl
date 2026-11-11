# T1_balance.jl — exact balance of EvolWBFC
#
#   (a) rest state over the mountain (v_initial = 0):
#       every FLUID velocity bitwise +0.0 and ρ = ρ̄(z) bitwise after every step
#   (b) uniform flow U over flat terrain (h_m = 0), including the inflow buffer, the ghosts and the
#       teleport from the outflow to the inflow:
#       v = (U, +0.0) and ρ = ρ̄(z) bitwise after every step
#
# Usage:
#   julia -t auto scripts/debug/fully_compressible/T1_balance.jl [nsteps]

include(joinpath(@__DIR__, "fc_common.jl"))

nsteps = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 300

function run_case(label::String, nsteps::Int; overrides...)
	global_params, sim_params = load_params("EvolWBFC"; overrides...)
	d = derived(global_params, sim_params)
	U = sim_params[:v_initial]
	println("\n", "="^70, "\n$label: dr = $(sim_params[:dr]), h_m = $(sim_params[:h_m]), U = $U m/s, $nsteps steps of $(round(d.dt; digits=4)) s\n", "="^70)

	sys = make_fc_system(global_params, sim_params)
	M.initialize_state!(sys, global_params, sim_params)
	x0 = Dict(p.id => p.x for p in sys.particles if p.type == M.FLUID)

	bad_v = 0       # FLUID particle-steps with v ≠ (U, +0.0) bitwise
	bad_ρ = 0       # FLUID particle-steps with ρ ≠ ρ̄(z) bitwise
	dv_max = 0.0
	δP_max = 0.0
	t_wall = @elapsed for k in 1:nsteps
		M.verlet_step!(sys, global_params, sim_params)
		for p in sys.particles
			p.type == M.FLUID || continue
			exact = reinterpret(UInt64, p.v[1]) == reinterpret(UInt64, U) && reinterpret(UInt64, p.v[2]) == 0x0000000000000000
			bad_v += !exact
			bad_ρ += p.ρ != p.ρ_bg
			dv_max = max(dv_max, hypot(p.v[1] - U, p.v[2]))
			δP_max = max(δP_max, abs(p.δP))
		end
	end

	n_fluid = count(p -> p.type == M.FLUID, sys.particles)
	n_new = count(p -> p.type == M.FLUID && !haskey(x0, p.id), sys.particles)
	n_teleported = count(p -> p.type != M.FLUID && haskey(x0, p.id), sys.particles) + n_new
	@printf("FLUID particles at the end: %d (initially %d), particles that left through the outflow or entered through the inflow: %d\n",
		n_fluid, length(x0), n_teleported)
	@printf("FLUID particle-steps with v ≠ (U, +0.0) bitwise: %d,  with ρ ≠ ρ̄(z) bitwise: %d\n", bad_v, bad_ρ)
	@printf("max |v − U ê_x| = %.3e m/s,  max |δP| = %.3e Pa,  %.3f s per step\n", dv_max, δP_max, t_wall / nsteps)
	println(bad_v == 0 && bad_ρ == 0 ? "RESULT ($label): exactly balanced." : "RESULT ($label): NOT exactly balanced.")
	return bad_v == 0 && bad_ρ == 0
end

ok_a = run_case("(a) rest state over the mountain", nsteps; v_initial=0.0)
ok_b = run_case("(b) uniform flow over flat terrain", nsteps; h_m=0.0)
println()
println(ok_a && ok_b ? "T1 RESULT: both states are preserved bitwise." : "T1 RESULT: failed.")
