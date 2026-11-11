using DrWatson
@quickactivate "SPH"

# load all packages used by formulations
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
using TOML

# ============================================================
# dummy parameters 
# ============================================================

dummy_global = Dict{Symbol, Any}(
	:g          => 9.81,
	:R_mass     => 287.05,
	:cp         => 1004.675,
	:cv         => 717.625,
	:γ          => 1.4,
	:R_gas      => 8.314,
	:T_bg       => 250.0,
	:ρ0         => 1.393,
	:N          => 0.14,
	:dom_height => 5e3,  # use a small domain
	:dom_length => 10e3, # use a small domain
	:a          => 5e3,
	:z_t        => 5e3,
)

dummy_sim = Dict{Symbol, Any}(
	:dr         => 1000.0,
	:dt_rel     => 0.1,
	:t_end      => 1.0,       # short simulation time
	:rho_floor  => 1e-6,
	:P_floor    => 1e-10,
	:η          => 1.8,
	:ϵ          => 0.01,
	:α          => 0.15,
	:β          => 0.3,
	:γ_r_rel    => 0.0,
	:h_m        => 100.0,
	:z_β        => 2.5e3,
	:v_initial  => 20.0,
	:model      => "precompile",
)

# ============================================================
# reach formulation needs to bu run 
# ============================================================

formulations = [
	"EvolWBSPH",
	"EvolWBPA",
	"EvolWBPTH",
	"EvolWBFC",
]

for name in formulations
	model_path = srcdir("formulations", "evolutionary", "$(name).jl")
	if isfile(model_path)
		include(model_path)
		M = getfield(Main, Symbol(name))
		println("Running dummy sim: $name")
		try
			M.run_sim(dummy_global, dummy_sim)
		catch e
			@warn "Precompile run failed for $name" exception=(e, catch_backtrace())
		end
		println("Precompiled: $name")
	else
		println("Warning: $model_path not found")
	end
end

println("-------sysimage precompilation complete-------")
