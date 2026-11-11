using DrWatson
@quickactivate "SPH"
using HDF5
using JLD2
using CairoMakie
using Colors
using Printf
using Statistics
using TOML

# Flow field and potential-temperature perturbation of an evolutionary run:
# three panels with w, u − U and δθ, interpolated from the FLUID particles onto a regular grid
# (the companion of emerging_wave.jl, which shows w with the isentropes).
#
# Usage: julia scripts/plotting/wave_perturbations.jl <run_dir> <frame_file.h5> <formulation> [global_params.toml]

const FLUID = 0.0
const MOUNTAIN = 1.0

# ── 2-D Wendland C² kernel with the support radius h ──────────────────────────
@inline function wendland(h::Float64, r::Float64)::Float64
    q = r / h
    q >= 1.0 && return 0.0
    return (1.0 - q)^4 * (1.0 + 4.0 * q)
end

# ── Shepard interpolation onto a regular grid, with the smoothing length of every particle ──
function shepard_grid(xp, zp, hp, fields, x_grid, z_grid)
    Nx, Nz = length(x_grid), length(z_grid)
    dx, dz = Float64(step(x_grid)), Float64(step(z_grid))
    sums = [zeros(Nx, Nz) for _ in fields]
    weight = zeros(Nx, Nz)

    for j in eachindex(xp)
        i_cen = round(Int, (xp[j] - first(x_grid)) / dx) + 1
        k_cen = round(Int, (zp[j] - first(z_grid)) / dz) + 1
        ix_half, iz_half = ceil(Int, hp[j] / dx) + 1, ceil(Int, hp[j] / dz) + 1
        for i in max(1, i_cen - ix_half):min(Nx, i_cen + ix_half)
            for k in max(1, k_cen - iz_half):min(Nz, k_cen + iz_half)
                w = wendland(hp[j], hypot(x_grid[i] - xp[j], z_grid[k] - zp[j]))
                weight[i, k] += w
                for (f, s) in zip(fields, sums)
                    s[i, k] += w * f[j]
                end
            end
        end
    end
    return [[weight[i, k] > 1e-6 ? Float32(s[i, k] / weight[i, k]) : NaN32 for i in 1:Nx, k in 1:Nz] for s in sums]
end

# ── Arguments ─────────────────────────────────────────────────────────────────
length(ARGS) >= 3 || error("Usage: julia wave_perturbations.jl <run_dir> <frame_file.h5> <formulation> [global_params.toml]")
run_dir     = ARGS[1]
frame_file  = ARGS[2]
formulation = ARGS[3]
global_file = length(ARGS) >= 4 ? ARGS[4] : projectdir("config", "global_params.toml")

frame_num = replace(basename(frame_file), r"frame_(\d+)\.h5" => s"\1")

meta = load(joinpath(run_dir, "metadata.jld2"))
dr, eta, U = Float64(meta["dr"]), Float64(meta["η"]), Float64(meta["v_initial"])
glob = TOML.parsefile(global_file)
K = glob["g"] / (glob["R_mass"] * glob["T_bg"])

# ── Frame ─────────────────────────────────────────────────────────────────────
t_snap, x_f, z_f, u_f, w_f, δθ_f, x_b, z_b = h5open(joinpath(run_dir, frame_file), "r") do fid
    pos   = read(fid["positions"])
    vel   = read(fid["velocities"])
    types = read(fid["types"])
    δθ    = read(fid["pot_temperatures_pert"])
    mask_fluid = types .== FLUID
    mask_bnd   = (types .== MOUNTAIN) .& (pos[2, :] .< 2.0e3)

    (read(HDF5.attributes(fid)["time"]),
     pos[1, mask_fluid], pos[2, mask_fluid],
     vel[1, mask_fluid] .- U, vel[2, mask_fluid], δθ[mask_fluid],
     pos[1, mask_bnd], pos[2, mask_bnd])
end

# interpolation radius: twice the lattice spacing at the height of the particle
h_f = 2.0 .* dr .* exp.(K .* max.(z_f, 0.0) ./ 2.0)

x_min, x_max = -60_000.0, 90_000.0
z_min, z_max = 0.0, 20_000.0
x_grid = range(x_min, x_max; length=450)
z_grid = range(z_min, z_max; length=200)

println("Interpolating w, u − U and δθ ...")
w_grid, u_grid, δθ_grid = shepard_grid(x_f, z_f, h_f, (w_f, u_f, δθ_f), x_grid, z_grid)

# ── Figure ────────────────────────────────────────────────────────────────────
tol_diverging = cgrad([
    colorant"#882255", colorant"#AA4466", colorant"#CC6677", colorant"#E8A5AE",
    colorant"#FFFFFF", colorant"#B8DCF0", colorant"#88CCEE", colorant"#5577BB", colorant"#332288"
])

# symmetric colour range from the data, rounded to one significant digit
function colour_limit(grid)
    values = abs.(filter(isfinite, vec(grid)))
    isempty(values) && return 1.0
    lim = quantile(values, 0.998)
    lim > 0.0 || return 1.0
    return round(lim; sigdigits=1)
end

panels = [(w_grid, L"w\;[\mathrm{m/s}]"), (u_grid, L"u - U\;[\mathrm{m/s}]"), (δθ_grid, L"\delta\theta\;[\mathrm{K}]")]

fig = Figure(size=(1000, 1150), fontsize=20)
for (row, (grid, label)) in enumerate(panels)
    ax = Axis(fig[row, 1], xlabel=row == length(panels) ? "distance x [km]" : "", ylabel="height z [km]")
    lim = colour_limit(grid)
    hm = heatmap!(ax, collect(x_grid) ./ 1e3, collect(z_grid) ./ 1e3, grid;
                  colormap=tol_diverging, colorrange=(-lim, lim), interpolate=true)
    scatter!(ax, x_b ./ 1e3, z_b ./ 1e3, color=(:darkgray, 1.0), markersize=4.5)
    xlims!(ax, x_min / 1e3, x_max / 1e3)
    ylims!(ax, -0.4, z_max / 1e3)
    Colorbar(fig[row, 2], hm, label=label)
end
Label(fig[0, :], @sprintf("%s,  t = %.0f s", formulation, t_snap), fontsize=22, halign=:left, tellwidth=false)

out_dir = plotsdir("evolutionary")
mkpath(out_dir)
out_filename = "fields_$(formulation)_$(frame_num).png"
save(joinpath(out_dir, out_filename), fig; px_per_unit = 3)
println("Saved -> $(joinpath(out_dir, out_filename))")
