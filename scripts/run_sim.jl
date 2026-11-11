using DrWatson
using Dates
@quickactivate "SPH"
using TOML

# ============================================================
# Runs a simulation with the given .toml config files 
# Includes some useful debugging info (learned the hard way)
# ============================================================


# ============================================================
# Parsing the input and loading
# works from the REPL (useful for local tuning and debugging): suffices to include this file without any further arguments
# the handling is done using DrWatson
# when submitting from command line (e.g. on a cluster) there is no REPL, so the .tomls are passed manually
# ============================================================

if isempty(ARGS)
	# REPL mode
	println("Using default project configs.")
	global_file = projectdir("config", "global_params.toml")
	sweep_file  = projectdir("config", "sim_params.toml")
	restart_dir = ""
else
	# command line mode
	if length(ARGS) < 2
		println("Usage: julia scripts/run_sim.jl global_params.toml sim_params.toml [restart_dir]")
		exit(1)
	end
	global_file = ARGS[1]
	sweep_file  = ARGS[2]
	restart_dir = length(ARGS) >= 3 ? ARGS[3] : ""
end

println("\n" * "="^60)
println("=== CONFIGURATION LOADED")
println("="^60)
println("Global Config:      $global_file")
println("Simulation Config:  $sweep_file")
if !isempty(restart_dir)
	println("Restart from:       $restart_dir")
end
println("-"^60)

# ============================================================
# Configuration loading & Mapping
# ============================================================

raw_global = TOML.parsefile(global_file)
raw_sweep  = TOML.parsefile(sweep_file)

# mapping the keys from TOML configs to symbols from the model source code 
# this would not be needed if one chose to name the variables using ascii only...
key_map = Dict(
	# physics
	"g" => :g,
	"R_mass" => :R_mass,
	"cp" => :cp,
	"cv" => :cv,
	"gamma" => :γ,
	"R_gas" => :R_gas,
	"T_bg" => :T_bg,
	"rho0" => :ρ0,
	"N" => :N,

	# geometry
	"dom_height" => :dom_height,
	"dom_length" => :dom_length,
	"h_m" => :h_m,
	"a" => :a,
	"z_t" => :z_t,
	"z_beta" => :z_β,

	# spatial and temporal resolution
	"dr" => :dr,
	"dt_rel" => :dt_rel,
	"t_end" => :t_end,

	# numerics
	"rho_floor" => :rho_floor,
	"P_floor" => :P_floor,
	"eta" => :η, 
	"epsilon" => :ϵ, 
	"alpha" => :α,
	"beta" => :β,
	"gamma_r_rel" => :γ_r_rel,
	"zeta" => :ζ,

	# initial conditions
	"v_initial" => :v_initial,

	# formulation
	"formulation" => :formulation
)

# produce a valid config dictionary for the model to use
function process_config(raw_dict)
	processed = Dict{Symbol, Any}()
	for (key, value) in raw_dict
		if value isa Dict
			for (sub_k, sub_v) in value
				final_k = get(key_map, sub_k, Symbol(sub_k))
				processed[final_k] = sub_v
			end
		else
			final_k = get(key_map, key, Symbol(key))
			processed[final_k] = value
		end
	end
	return processed
end

# build the global_params config
global_params = process_config(raw_global)

# build the sim_params config
temp_sweep = process_config(raw_sweep)
sweep_params_lists = Dict{Symbol, Vector}()
model_list = String[]

# we differentiate between the formulation entries and other paramaters 
# from the .toml file

for (k, v) in temp_sweep
	if k == :formulation
		# formulation
		global model_list = (v isa Vector) ? v : [v]
	else
		# other paramateres
		sweep_params_lists[k] = (v isa Vector) ? v : [v]
	end
end

if isempty(model_list)
	println("Error: No 'formulation' specified in sweep config!")
	exit(1)
end

# ============================================================
# Run the simulation with the chosen models and parameters
# ============================================================

println("\n Starting the batch job. Models to run: $model_list")

for model_name in model_list
	println("\n\n" * "="^60)
	println("=== LOADING MODEL: $model_name")
	println("="^60)
	# routing to the model source code (try evolutionary first, then stationary)
	model_path = srcdir("formulations", "evolutionary", "$(model_name).jl")
	if !isfile(model_path)
		model_path = srcdir("formulations", "stationary", "$(model_name).jl")
	end
	if !isfile(model_path)
		println("     Warning: Model file not found in evolutionary/ or stationary/.")
		exit(1)
	end

	# source the module
	include(model_path)

	# extract the module name from the source file
	try
		global ModelModule = getfield(Main, Symbol(model_name))
	catch
		println("    Error: Module '$model_name' not found inside file.")
		exit(1)
	end

	# sweep through the parameters for the chosen formulation
	# a nice addition: if we add the model name to the params, DrWatson saves it in the folder name

	current_sweep_lists = copy(sweep_params_lists)
	current_sweep_lists[:model] = [model_name] 

	# gather all the simulations to be run with the single formulation
	job_list = dict_list(current_sweep_lists)
	println("    Simulations to run: $(length(job_list))")

	for (i, sim_params) in enumerate(job_list)
		# Clear separator for each run
		println("\n" * "-"^60)
		println("### [Run $i/$(length(job_list))]")
		println("-"^60)

		println("Parameters:")
		for (k, v) in sort(collect(sim_params), by=x->x[1])
			# Aligned printing
			println(rpad("  __ $k", 20), " = $v")
		end
		println("") # Empty line for spacing before logs
		flush(stdout) # force the output

		# run the model (pass restart_dir if resuming from checkpoint)
		output_path = ModelModule.run_sim(global_params, sim_params; restart_dir=restart_dir)
		println("\n>>> Saved: $output_path")

		# for job chaining: write the run_dir to CHAIN_FILE so subsequent jobs know where to restart
		chain_file = get(ENV, "CHAIN_FILE", "")
		if !isempty(chain_file)
			open(chain_file, "w") do io
				write(io, output_path)
			end
			println("Chain file updated: $chain_file -> $output_path")
		end
	end
end

println("\n" * "="^60)
println("=== END OF THE ROAD. ALL RUNS FINISHED.")
println("="^60)
