"""
	AbstractParticle

Abstract supertype for smoothed particles. Any structure with `AbstractParticle`
supertype is expected to:

* be mutable
* have fields `x::Float64` and `y::Float64`
"""
abstract type AbstractParticle end


"""
	@define_particle(name::Symbol, vars::Symbol...)

Macro for quick definition of particles. Defines a `struct` with given `name`
that is a subtype of [`AbstractParticle`](@ref) and has fields `x`, `y` and
other fields specified in vars. Comes equipped with a default constructor that
takes `x`, `y` as positional and `vars` as optionl keyword arguments (default
to `0.0`)

# Example:
```
@define_particle Particle Dx Dy rho
```
is equivalent to
```
mutable struct Particle <: AbstractParticle
	x::Float64
	y::Float64
	Dx::Float64
	Dy::Float64
	rho::Float64
	Particle(x, y; Dx = 0.0, Dy = 0.0, rho = 0.0) = new(x, y, Dx, Dy, rho)
end
```
"""
macro define_particle(name::Symbol, vars::Symbol...)
	code = "mutable struct "*string(name)*" <: AbstractParticle\n"
	code *= "\tx::Float64\n"
	code *= "\ty::Float64\n"
	for var in vars
		code *= "\t"*string(var)*"::Float64\n"
	end
	code *= "\t"*string(name)*"(x,y;"
	for var in vars
		code *= string(var)*"= 0.,"
	end
	code *= ") = new(x,y,"
	for var in vars
		code *= string(var)*","
	end
	code *= ")\nend"
	return Meta.parse(code)
end


"""
	ParticleSystem(T::Type;
                   xlims::Tuple{Float64, Float64},
                   ylims::Tuple{Float64, Float64},
                   h::Float64)

An immutable struct that contains all vital information about the simulation.
The constructor specifies that:
* the simulation will use particles of type `T <: AbstractParticle`,
* it will take place inside a cartesian product `xlims` times `ylims`,
* Particles are considered neighbours if their distance is less than `h`.
Individual particles can be accessed through attribute `ParticleSystem.particles`
"""
struct ParticleSystem{T <: AbstractParticle}
	dr::Float64
	h::Float64

	xmin::Float64
	xmax::Float64
	ymin::Float64
	ymax::Float64

	i_phase::Int64
	j_phase::Int64
	i_max::Int64
	j_max::Int64
	key_max::Int64

	particle_type::DataType
	particles::Vector{T}
	cell_list::Vector{Vector{Union{Missing, T}}}
	cell_lock::Vector{Threads.ReentrantLock}	#locks used in parallel cell_list generation

	ParticleSystem(T::Type, xrange::Tuple{Float64, Float64}, yrange::Tuple{Float64, Float64}, dr::Float64, h::Float64 = 2.0*dr) =
	begin
			@assert(dr > 0.0, "invalid ParticleSystem declaration! (dr must be a positive float)")
			@assert(h  > 0.0, "invalid ParticleSystem declaration! (h must be a positive float)")
			@assert(T <: AbstractParticle, "invalid ParticleSystem declaration! ("*T*" is not an AbstractParticle subtype)")
			@assert(hasfield(T, :x)  , "invalid ParticleSystem declaration! (particles must have a field `x::Float64`)")
			@assert(hasfield(T, :y)  , "invalid ParticleSystem declaration! (particles must have a field `y::Float64`)")

			xmin = minimum(xrange)
			xmax = maximum(xrange)
			ymin = minimum(yrange)
			ymax = maximum(yrange)
			i_phase = Int64(floor(xmin/h))
			j_phase = Int64(floor(ymin/h))
			i_max 	= Int64(floor(xmax/h)) - i_phase + 1
			j_max 	= Int64(floor(ymax/h)) - j_phase + 1
			key_max   = i_max*j_max

			cell_list = Vector{Vector{Union{Missing, T}}}(undef, key_max)
			cell_lock = Vector{Threads.ReentrantLock}(undef, key_max)
			for key in eachindex(cell_list)
				cell_list[key] = Vector{Union{Missing, T}}(missing, Int64(ceil((h/dr + 1.0)^2)))
				cell_lock[key] = Threads.ReentrantLock()
			end

			return new{T}(dr, h, xmin, xmax, ymin, ymax, i_phase, j_phase, i_max, j_max, key_max, T, Vector{T}(), cell_list, cell_lock)
	end
end

ParticleSystem(T::Type;
			   xlims::Tuple{Float64, Float64},
			   ylims::Tuple{Float64, Float64},
			   h::Float64
) = ParticleSystem(T, xlims, ylims, h/4, h)



"""
	ScalarField(sys::ParticleSystem, varS::Symbol, name::String = string(varS))

Creates an abstract array whose ``n``-th element is the value of scalar `varS`
of ``n``-th particle in `sys`. The `name` will be used when saved into a file.
"""
struct ScalarField <: AbstractArray{Float64, 1}
	sys::ParticleSystem
	varS::Symbol
	name::String
	ScalarField(sys::ParticleSystem, varS::Symbol, name::String = string(varS)) = new(sys, varS, name)
end

Base.size(f::ScalarField) = (length(f.sys.particles),)
Base.IndexStyle(::Type{<:ScalarField}) = IndexLinear()
Base.getindex(f::ScalarField, i::Int64) = getproperty(f.sys.particles[i], f.varS)
Base.setindex!(f::ScalarField, val::Float64, k::Int64) = setproperty!(f.sys.particles[k], f.varS, val)

"""
	VectorField(sys::ParticleSystem, varS::Tuple{Symbol, Symbol}, name::String = string(varS))

Creates an abstract array whose ``n``-th element is the value of vector `varS`
of ``n``-th particle in `sys`.  The `name` will be used when saved into a file.
"""
struct VectorField <: AbstractArray{Float64, 2}
	sys::ParticleSystem
	varS::Tuple{Symbol, Symbol}
	name::String
	VectorField(sys::ParticleSystem, varS::Tuple{Symbol, Symbol}, name::String = string(varS)) = new(sys, varS, name)
end

Base.size(f::VectorField) = (3, length(f.sys.particles))

Base.getindex(f::VectorField, I::Vararg{Int64, 2}) = begin
	if I[1] != 3
		return getproperty(f.sys.particles[I[2]], f.varS[I[1]])
	end
	return 0.
end

Base.setindex!(f::VectorField, val::Float64, I::Vararg{Int64, 2}) = begin
	if I[1] != 3
		setproperty!(f.sys.particles[I[2]], val, f.varS[I[1]])
	end
end

const DataField = Union{ScalarField, VectorField}

function unary_void(p::AbstractParticle)
end

function binary_void(p::AbstractParticle, q::AbstractParticle)
end
