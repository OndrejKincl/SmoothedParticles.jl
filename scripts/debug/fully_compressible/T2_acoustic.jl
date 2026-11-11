# T2_acoustic.jl — compressional force and speed of sound in EvolWBFC
#
# Atmosphere at rest over flat terrain. A Gaussian density excess (ρ ← ρ̄ (1 + ε exp(-r²/σ²)), θ unchanged,
# hence δP = c² δρ > 0) is released. It must expand and radiate a ring travelling with c = sqrt(γ R T_bg).
# The ring is timed by the passage of the maximum of δP at the FLUID particles of a horizontal strip.
#
# Usage:
#   julia -t auto scripts/debug/fully_compressible/T2_acoustic.jl [t_end]

include(joinpath(@__DIR__, "fc_common.jl"))
include(joinpath(@__DIR__, "..", "frozen_background", "plot_style.jl"))

t_end = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 60.0
global_params, sim_params = load_params("EvolWBFC"; v_initial=0.0, h_m=0.0, t_end=t_end)
d = derived(global_params, sim_params)
nsteps = round(Int, t_end / d.dt)

sys = make_fc_system(global_params, sim_params)
M.initialize_state!(sys, global_params, sim_params)

# ── the pulse ─────────────────────────────────────────────────────────────────
x_c, z_c = 100e3, 8e3
ε, σ = 1e-4, 2.5e3
for p in sys.particles
	if p.type == M.FLUID
		p.ρ = p.ρ_bg * (1.0 + ε * exp(-((p.x[1] - x_c)^2 + (p.x[2] - z_c)^2) / σ^2))
	end
end
M.evaluate_state!(sys, 0.0, global_params, sim_params)

centre = sys.particles[nearest_fluid(sys, x_c, z_c)]
println("T2 — acoustic pulse, EvolWBFC: dr = $(sim_params[:dr]), η = $(sim_params[:η]), $nsteps steps of $(round(d.dt; digits=4)) s")
@printf("pulse at (%.0f, %.0f) m: δρ/ρ̄ = %.1e, σ = %.0f m;  at the centre δP = %.4f Pa (γ P̄ ε = %.4f Pa), outward acceleration next to it: max|Dv| = %.3e m/s²\n",
	x_c, z_c, ε, σ, centre.δP, global_params[:γ] * centre.P_bg * ε, maximum(SmoothedParticles.norm(p.Dv) for p in sys.particles if p.type == M.FLUID))

# ── time integration: time of the maximum of δP at every particle of the strip ─
strip = [i for (i, p) in enumerate(sys.particles) if p.type == M.FLUID && abs(p.x[2] - z_c) < 1.0e3 && 2σ < abs(p.x[1] - x_c) < 40e3]
dist = [abs(sys.particles[i].x[1] - x_c) for i in strip]
peak = zeros(length(strip))
t_peak = zeros(length(strip))
δP_centre = Float64[centre.δP]
times = Float64[0.0]
for k in 1:nsteps
	M.verlet_step!(sys, global_params, sim_params)
	for (j, i) in enumerate(strip)
		if sys.particles[i].δP > peak[j]
			peak[j] = sys.particles[i].δP
			t_peak[j] = k * d.dt
		end
	end
	push!(δP_centre, centre.δP)
	push!(times, k * d.dt)
end

# particles which the maximum has passed well before the end of the run
sel = [j for j in eachindex(strip) if 0.0 < t_peak[j] < t_end - 10.0]
t, r = t_peak[sel], dist[sel]
slope = sum((t .- mean(t)) .* (r .- mean(r))) / sum((t .- mean(t)) .^ 2)
intercept = mean(r) - slope * mean(t)
println()
@printf("ring speed (least squares over %d particles, %.1f–%.1f km from the centre): %.1f m/s\n", length(sel), minimum(r) / 1e3, maximum(r) / 1e3, slope)
@printf("c = sqrt(γ R T_bg) = %.1f m/s   →   ratio %.4f\n", d.c_s, slope / d.c_s)
@printf("δP at the centre: %.4f Pa at t = 0, minimum %.4f Pa at t = %.1f s, %.4f Pa at the end\n", δP_centre[1], minimum(δP_centre), times[argmin(δP_centre)], δP_centre[end])
@printf("max |v| at the end: %.3e m/s\n", maximum(SmoothedParticles.norm(p.v) for p in sys.particles if p.type == M.FLUID))

# ── figure ────────────────────────────────────────────────────────────────────
set_style!()
fig = Figure(size=(1150, 470))
ax1 = Axis(fig[1, 1]; xlabel="t [s]", ylabel="δP [Pa]", title="pressure perturbation at the centre of the pulse")
hlines!(ax1, [0.0]; color=C_INK2, linewidth=1)
lines!(ax1, times, δP_centre; color=C_NEW)
ax2 = Axis(fig[1, 2]; xlabel="time of the maximum of δP [s]", ylabel="distance from the centre [km]", title="passage of the ring, strip |z − z_c| < 1 km")
scatter!(ax2, t, r ./ 1e3; color=C_NEW, markersize=8, strokecolor=C_SURFACE, strokewidth=1, label="FLUID particles")
ts = [minimum(t), maximum(t)]
lines!(ax2, ts, (intercept .+ slope .* ts) ./ 1e3; color=C_INK, label="fit: $(round(slope; digits=1)) m/s")
lines!(ax2, ts, (intercept .+ slope * ts[1] .+ d.c_s .* (ts .- ts[1])) ./ 1e3; color=C_REF, linestyle=:dash, label="c = sqrt(γRT) = $(round(d.c_s; digits=1)) m/s")
axislegend(ax2; position=:rb)
Label(fig[0, :], "T2 — acoustic pulse, EvolWBFC (dr = $(Int(sim_params[:dr])) m, η = $(sim_params[:η]))"; fontsize=18, halign=:left, tellwidth=false)
save_fc_plot("T2_acoustic.png", fig)
