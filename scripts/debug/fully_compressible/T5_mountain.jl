# T5_mountain.jl — the mountain configuration with EvolWBFC
#
# Runs the thesis mountain configuration (dr = 520, v_initial = 20, h_m = 100, t_end = 1200) through
# EvolWBFC.run_sim and summarises the frames: time series of the vertical velocity, of the deviation from
# the inflow velocity and of the particle spacing, and the w-field at a few frames.
#
# Usage:
#   julia -t auto scripts/debug/fully_compressible/T5_mountain.jl [t_end] [key=value …]
#
# e.g. `γ_r_rel=0.2` switches the sponge layer on, `α=0.05 β=0.1` changes the artificial viscosity.

include(joinpath(@__DIR__, "fc_common.jl"))
include(joinpath(@__DIR__, "..", "frozen_background", "plot_style.jl"))
using HDF5
using CSV, DataFrames

positional = [a for a in ARGS if !occursin("=", a)]
overrides = Dict{Symbol,Any}()
for a in ARGS
	occursin("=", a) || continue
	key, value = split(a, "="; limit=2)
	overrides[Symbol(key)] = parse(Float64, value)
end
tag = isempty(overrides) ? "" : "_" * join(("$(k)$(v)" for (k, v) in sort(collect(overrides); by=first)), "_")

global_params, sim_params = load_params("fully_compressible/T5_mountain$(tag)"; overrides...)
if length(positional) >= 1
	sim_params[:t_end] = parse(Float64, positional[1])
end
d = derived(global_params, sim_params)
U, dr = sim_params[:v_initial], sim_params[:dr]
println("T5 — mountain run, EvolWBFC: dr = $dr, U = $U, h_m = $(sim_params[:h_m]), α = $(sim_params[:α]), γ_r_rel = $(sim_params[:γ_r_rel]), t_end = $(sim_params[:t_end]) s")

t_wall = @elapsed run_dir = M.run_sim(global_params, sim_params)

# ── frames ────────────────────────────────────────────────────────────────────
function read_fluid(path)
	h5open(path, "r") do fid
		mask = read(fid["types"]) .== M.FLUID
		pos = read(fid["positions"])
		vel = read(fid["velocities"])
		(t = read(HDF5.attributes(fid)["time"]), x = pos[1, mask], z = pos[2, mask], u = vel[1, mask], w = vel[2, mask],
			δθ = read(fid["pot_temperatures_pert"])[mask])
	end
end

# smallest distance to another FLUID particle, in units of the local lattice spacing s(z) = dr·exp(Kz/2)
function nearest_neighbour(x, z)
	order = sortperm(x)
	xs, zs = x[order], z[order]
	dmin = fill(Inf, length(xs))
	Threads.@threads for i in eachindex(xs)
		s = dr * exp(d.K * max(zs[i], 0.0) / 2.0)
		best = Inf
		j = i - 1
		while j >= 1 && xs[i] - xs[j] < min(3s, best)
			best = min(best, hypot(xs[i] - xs[j], zs[i] - zs[j]))
			j -= 1
		end
		j = i + 1
		while j <= length(xs) && xs[j] - xs[i] < min(3s, best)
			best = min(best, hypot(xs[i] - xs[j], zs[i] - zs[j]))
			j += 1
		end
		dmin[i] = best / s
	end
	return dmin
end

frames = sort(filter(f -> occursin(r"^frame_\d{4}\.h5$", basename(f)), readdir(run_dir; join=true)))
series = DataFrame(frame=Int[], t=Float64[], n_fluid=Int[], w_max=Float64[], w_rms=Float64[], dv_max=Float64[], δθ_max=Float64[], nn_min=Float64[])
for (k, path) in enumerate(frames)
	f = read_fluid(path)
	push!(series, (k - 1, f.t, length(f.x), maximum(abs, f.w), sqrt(mean(abs2, f.w)), maximum(hypot.(f.u .- U, f.w)), maximum(abs, f.δθ), minimum(nearest_neighbour(f.x, f.z))))
end
CSV.write(joinpath(fc_datadir(), "T5_series$(tag).csv"), series)

steps = (length(frames) - 1) * max(1, round(Int, sim_params[:t_end] / 100 / d.dt))
println()
@printf("%d steps in %.1f s of wall time (%.3f s per step, including frames and checkpoints)\n", steps, t_wall, (mtime(frames[end]) - mtime(frames[1])) / steps)
@printf("  %5s %8s %8s %12s %12s %14s %12s %10s\n", "frame", "t [s]", "N_fluid", "max|w|", "rms w", "max|v − U ê_x|", "max|δθ| [K]", "nn min/s")
for row in eachrow(series)
	(row.frame % 10 == 0 || row.frame == series.frame[end]) || continue
	@printf("  %5d %8.1f %8d %12.4e %12.4e %14.4e %12.4e %10.3f\n", row.frame, row.t, row.n_fluid, row.w_max, row.w_rms, row.dv_max, row.δθ_max, row.nn_min)
end

# ── figures ───────────────────────────────────────────────────────────────────
set_style!()
fig = Figure(size=(1200, 420))
ax1 = Axis(fig[1, 1]; xlabel="t [s]", ylabel="[m/s]", yscale=log10, title="vertical velocity of the FLUID particles")
lines!(ax1, series.t[2:end], series.w_max[2:end]; color=C_NEW, label="max |w|")
lines!(ax1, series.t[2:end], series.w_rms[2:end]; color=C_OLD, label="r.m.s. w")
axislegend(ax1; position=:rb)
ax2 = Axis(fig[1, 2]; xlabel="t [s]", ylabel="distance / local lattice spacing", title="smallest nearest-neighbour distance")
lines!(ax2, series.t, series.nn_min; color=C_NEW)
ylims!(ax2, 0.0, 1.0)
Label(fig[0, :], "T5 — mountain run, EvolWBFC (dr = $(Int(dr)) m, U = $(Int(U)) m/s, h_m = $(Int(sim_params[:h_m])) m, α = $(sim_params[:α]), γ_r_rel = $(sim_params[:γ_r_rel]))";
	fontsize=18, halign=:left, tellwidth=false)
save_fc_plot("T5_series$(tag).png", fig)

shown = [k for k in (20, 40, 60, 80, 100) if k < length(frames)]
if !isempty(shown)
	fields = [read_fluid(frames[k+1]) for k in shown]
	w_lim = maximum(quantile(abs.(f.w), 0.995) for f in fields)
	fig = Figure(size=(1150, 250 * length(shown) + 110))
	for (row, f) in enumerate(fields)
		ax = Axis(fig[row, 1]; xlabel=row == length(shown) ? "x [km]" : "", ylabel="z [km]", title="t = $(round(Int, f.t)) s", xgridvisible=false, ygridvisible=false)
		scatter!(ax, f.x ./ 1e3, f.z ./ 1e3; color=f.w, colormap=DIVERGING, colorrange=(-w_lim, w_lim), markersize=4)
		xlims!(ax, -100, 100)
		ylims!(ax, 0, 26)
	end
	Colorbar(fig[1:length(shown), 2]; colormap=DIVERGING, limits=(-w_lim, w_lim), label="w [m/s]")
	Label(fig[0, :], "T5 — vertical velocity, EvolWBFC (colour range ±$(round(w_lim; sigdigits=2)) m/s)"; fontsize=18, halign=:left, tellwidth=false)
	save_fc_plot("T5_w_fields$(tag).png", fig)
end
println("run directory: $run_dir")
