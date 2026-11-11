# T6_restart.jl — is a restarted EvolWBFC run bitwise identical to an uninterrupted one?
#
# Mountain configuration (v_initial = 20), 2n steps:
#   A: 0 → 2n in one go
#   B: 0 → n, then restart from checkpoint.jld2 → 2n
# The final checkpoints are compared bitwise, field by field.
#
# Usage:
#   julia -t auto scripts/debug/fully_compressible/T6_restart.jl [n]

include(joinpath(@__DIR__, "fc_common.jl"))
using JLD2

n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 150
global_params, sim_params = load_params("fully_compressible/T6_restart")
d = derived(global_params, sim_params)

# with t_end = n·dt and 2n·dt the checkpoints fall on the steps n and 2n
sim_A = copy(sim_params); sim_A[:t_end] = 2n * d.dt
sim_B = copy(sim_params); sim_B[:t_end] = n * d.dt

println("="^70, "\nT6 — A: 0 → $(2n) steps\n", "="^70)
dir_A = M.run_sim(global_params, sim_A)
sleep(1.1)
println("="^70, "\nT6 — B: 0 → $n steps\n", "="^70)
dir_B = M.run_sim(global_params, sim_B)
println("="^70, "\nT6 — B: restart, $n → $(2n) steps\n", "="^70)
M.run_sim(global_params, sim_A; restart_dir=dir_B)

ckpt_A = load(joinpath(dir_A, "checkpoint.jld2"))
ckpt_B = load(joinpath(dir_B, "checkpoint.jld2"))
pa, pb = ckpt_A["particles"], ckpt_B["particles"]
println()
println("final checkpoints: A at step $(ckpt_A["k"]), B at step $(ckpt_B["k"]); particles: $(length(pa)) / $(length(pb))")

bits(v::Float64) = (reinterpret(UInt64, v),)
bits(v::RealVector) = (reinterpret(UInt64, v[1]), reinterpret(UInt64, v[2]), reinterpret(UInt64, v[3]))
bits(v::Int64) = (v,)
all_equal = length(pa) == length(pb) && ckpt_A["k"] == ckpt_B["k"]
for name in fieldnames(M.Particle)
	n_diff = count(i -> bits(getfield(pa[i], name)) != bits(getfield(pb[i], name)), 1:min(length(pa), length(pb)))
	global all_equal &= n_diff == 0
	n_diff == 0 || @printf("  field %-14s differs for %d particles\n", name, n_diff)
end
println(all_equal ? "T6 RESULT: the restarted run is bitwise identical to the uninterrupted one (all $(length(fieldnames(M.Particle))) particle fields)." :
        "T6 RESULT: the restarted run DIFFERS from the uninterrupted one.")
