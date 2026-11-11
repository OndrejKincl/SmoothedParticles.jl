# T3_gravity_wave.jl — buoyancy frequency in EvolWBFC
#
# Atmosphere at rest over flat terrain, between the ground and the rigid lid at z = D. The first vertical
# mode of the internal gravity waves of the isothermal atmosphere is set as a vertical displacement
#
#   ξ = ξ_0 exp(z/(2H)) sin(π z/D) cos(k_x (x − x_c)) exp(−(x − x_c)²/(2σ_x²)),
#
# i.e. δθ = −ξ dθ̄/dz at constant pressure (ρ = ρ̄ θ̄/θ, so δP = 0 and nothing acoustic is excited at
# leading order), with the fluid at rest. It must oscillate with
#
#   ω² = N² k_x² / (k_x² + (π/D)² + 1/(4H²)),   N² = g²/(c_p T_bg),   H = R T_bg / g.
#
# The period is read from the zero crossings of δθ and w at mid-height in the centre of the packet.
#
# Usage:
#   julia -t auto scripts/debug/fully_compressible/T3_gravity_wave.jl [t_end] [eta] [alpha]

include(joinpath(@__DIR__, "fc_common.jl"))
include(joinpath(@__DIR__, "..", "frozen_background", "plot_style.jl"))
using CSV, DataFrames

t_end = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 900.0
global_params, sim_params = load_params("EvolWBFC"; v_initial=0.0, h_m=0.0, t_end=t_end)
if length(ARGS) >= 2
	sim_params[:η] = parse(Float64, ARGS[2])
end
if length(ARGS) >= 3
	sim_params[:α] = parse(Float64, ARGS[3])
	sim_params[:β] = 2.0 * sim_params[:α]
end
tag = (length(ARGS) >= 2 ? "_eta$(ARGS[2])" : "") * (length(ARGS) >= 3 ? "_alpha$(ARGS[3])" : "")
d = derived(global_params, sim_params)
atm = M.Atmosphere(global_params)
nsteps = round(Int, t_end / d.dt)

sys = make_fc_system(global_params, sim_params)
M.initialize_state!(sys, global_params, sim_params)

# ── the wave packet ───────────────────────────────────────────────────────────
g, R_mass, γ, T_bg = atm.g, atm.R_mass, atm.γ, atm.T_bg
N = g / sqrt(global_params[:cp] * T_bg)
H = R_mass * T_bg / g
D = global_params[:dom_height]

x_c = 100e3
λ_x, σ_x = 20e3, 30e3
ξ_0 = 1.0
k_x = 2π / λ_x
k_z = π / D
for p in sys.particles
	if p.type == M.FLUID
		x, z = p.x[1], p.x[2]
		ξ = ξ_0 * exp(z / (2H)) * sin(k_z * z) * cos(k_x * (x - x_c)) * exp(-(x - x_c)^2 / (2σ_x^2))
		p.θ = p.θ_bg - ξ * p.θ_bg * N^2 / g
		p.ρ = p.ρ_bg * p.θ_bg / p.θ
	end
end
M.evaluate_state!(sys, 0.0, global_params, sim_params)

# theory
ω = N * k_x / sqrt(k_x^2 + k_z^2 + 1.0 / (4.0 * H^2))

z_p = 0.5 * D
probe = sys.particles[nearest_fluid(sys, x_c, z_p)]
spacing = sim_params[:dr] * exp(atm.K * z_p / 2.0)
println("T3 — internal gravity wave, EvolWBFC: dr = $(sim_params[:dr]), η = $(sim_params[:η]), α = $(sim_params[:α]), $nsteps steps of $(round(d.dt; digits=4)) s")
@printf("first vertical mode between z = 0 and the lid at %.0f km: λ_x = %.0f km (%.1f lattice spacings at the probe), σ_x = %.0f km, ξ_0 = %.1f m\n",
	D / 1e3, λ_x / 1e3, λ_x / spacing, σ_x / 1e3, ξ_0)
@printf("probe at (%.1f, %.2f) km: δθ = %.5f K, δP = %.2e Pa, max |δP| over the fluid = %.2e Pa\n", probe.x[1] / 1e3, probe.x[2] / 1e3, probe.δθ, probe.δP,
	maximum(abs(p.δP) for p in sys.particles if p.type == M.FLUID))
@printf("theory: N = %.5f rad/s, ω = %.5f rad/s = %.3f N, period %.1f s\n", N, ω, ω / N, 2π / ω)

# ── time integration ──────────────────────────────────────────────────────────
record = DataFrame(t=Float64[0.0], δθ=Float64[probe.δθ], w=Float64[0.0], z=Float64[probe.x[2]], δP=Float64[probe.δP], w_max=Float64[0.0])
t_wall = @elapsed for k in 1:nsteps
	M.verlet_step!(sys, global_params, sim_params)
	push!(record, (k * d.dt, probe.δθ, probe.v[2], probe.x[2], probe.δP, maximum(abs(p.v[2]) for p in sys.particles if p.type == M.FLUID)))
end
CSV.write(joinpath(fc_datadir(), "T3_probe$(tag).csv"), record)

# zero crossings by linear interpolation
function crossings(t, f)
	out = Float64[]
	for i in 1:length(f)-1
		if f[i] == 0.0 || f[i] * f[i+1] < 0.0
			push!(out, t[i] - f[i] * (t[i+1] - t[i]) / (f[i+1] - f[i]))
		end
	end
	return out
end
c_θ = crossings(record.t, record.δθ)
c_w = crossings(record.t[2:end], record.w[2:end])
half_periods = vcat(diff(c_θ), diff(c_w))

println()
@printf("%.3f s per step\n", t_wall / nsteps)
println("zero crossings of δθ at the probe [s]: ", join((@sprintf("%.1f", c) for c in c_θ), ", "))
println("zero crossings of w  at the probe [s]: ", join((@sprintf("%.1f", c) for c in c_w), ", "))
if !isempty(half_periods)
	T_meas = 2.0 * mean(half_periods)
	@printf("measured period: %.1f s (from %d half periods, %.1f–%.1f s)   →   ω = %.5f rad/s = %.3f N\n",
		T_meas, length(half_periods), 2minimum(half_periods), 2maximum(half_periods), 2π / T_meas, 2π / T_meas / N)
	@printf("theory:          %.1f s   →   measured/theory = %.3f\n", 2π / ω, T_meas / (2π / ω))
end
# amplitude after one period, where δθ returns to an extremum of the same sign
i_T = argmin(abs.(record.t .- 2π / ω))
lo, hi = max(1, i_T - round(Int, 0.15 * 2π / ω / d.dt)), min(nrow(record), i_T + round(Int, 0.15 * 2π / ω / d.dt))
if hi > lo
	back = sign(record.δθ[1]) * maximum(sign(record.δθ[1]) .* record.δθ[lo:hi])
	@printf("δθ at the probe: %.5f K at t = 0, extremum after one period: %.5f K (ratio %.3f)\n", record.δθ[1], back, back / record.δθ[1])
end
@printf("max |w| over the fluid: %.4e m/s,  vertical excursion of the probe: %.2f m\n", maximum(record.w_max), maximum(record.z) - minimum(record.z))

# ── figure ────────────────────────────────────────────────────────────────────
set_style!()
fig = Figure(size=(1150, 470))
ax1 = Axis(fig[1, 1]; xlabel="t [s]", ylabel="δθ [K]", title="potential-temperature perturbation at the probe")
hlines!(ax1, [0.0]; color=C_INK2, linewidth=1)
lines!(ax1, record.t, record.δθ; color=C_NEW, label="EvolWBFC")
lines!(ax1, record.t, record.δθ[1] .* cos.(ω .* record.t); color=C_REF, linestyle=:dash, label="δθ(0) cos(ωt), linear theory without damping")
axislegend(ax1; position=:rt)
ax2 = Axis(fig[1, 2]; xlabel="t [s]", ylabel="w [m/s]", title="vertical velocity at the probe")
hlines!(ax2, [0.0]; color=C_INK2, linewidth=1)
lines!(ax2, record.t, record.w; color=C_NEW)
Label(fig[0, :], "T3 — internal gravity wave, EvolWBFC (dr = $(Int(sim_params[:dr])) m, η = $(sim_params[:η]), α = $(sim_params[:α]))";
	fontsize=18, halign=:left, tellwidth=false)
save_fc_plot("T3_gravity_wave$(tag).png", fig)
