using  StaticArrays

const Vec2 = SArray{Tuple{2},Float64,1,2}
const Mat2 = SArray{Tuple{2,2},Float64,2,4}
const INIT_CELL_SIZE = 0

"""
	AbstractParticle

Abstract supertype for smoothed particles. Any structure with `AbstractParticle`
supertype is expected to:

* be mutable
* have field `x::Vec2` (particle position)
"""
abstract type AbstractParticle end

"""
    Shape

Supertype for geometrical shapes.
"""
abstract type Shape end

"""
	ParticleSystem(T::Type; h::Float64)

An immutable struct that contains all vital information about the simulation.
The constructor specifies that:
* the simulation will use particles of type `T <: AbstractParticle`,
* Particles are considered neighbours if their distance is less than `h`.
* 
"""
struct ParticleSystem{T <: AbstractParticle}
	h::Float64

	x1_min::Float64
	x2_min::Float64
	x1_max::Float64
	x2_max::Float64

	i_phase::Int64
	j_phase::Int64
	i_max::Int64
	j_max::Int64
	key_max::Int64

	particle_type::DataType
	particles::Vector{T}
	cell_list::Vector{Vector{Union{Missing, T}}}
	cell_lock::Vector{Threads.ReentrantLock}	#locks used in parallel cell_list generation

	ParticleSystem(T::DataType, dom::Shape, h::Float64) =
	begin
			@assert(h  > 0.0, "invalid ParticleSystem declaration! (h must be a positive float)")
			@assert(T <: AbstractParticle, "invalid ParticleSystem declaration! ("*string(T)*" is not an AbstractParticle subtype)")
			@assert(hasfield(T, :x)  , "invalid ParticleSystem declaration! (particles must have a field `x::Vec2`)")
			rect = boundarybox(dom)
			x1_min = rect.x1_min
			x2_min = rect.x2_min
			x1_max = rect.x1_max
			x2_max = rect.x2_max
			i_phase = Int64(floor(x1_min/h))
			j_phase = Int64(floor(x2_min/h))
			i_max 	= Int64(floor(x1_max/h)) - i_phase + 1
			j_max 	= Int64(floor(x2_max/h)) - j_phase + 1
			key_max = i_max*j_max

			particles = Vector{T}()
			cell_list = Vector{Vector{Union{Missing, T}}}(undef, key_max)
			cell_lock = Vector{Threads.ReentrantLock}(undef, key_max)
			for key in 1:key_max
				cell_list[key] = Vector{Union{Missing, T}}(missing, INIT_CELL_SIZE)
				cell_lock[key] = Threads.ReentrantLock()
			end
			return new{T}(h, 
				x1_min, x2_min, 
				x1_max, x2_max, 
				i_phase, j_phase, 
				i_max, j_max, key_max, 
				T, particles,
				cell_list, cell_lock
			)
	end
end

"""
	ParticleField(sys::ParticleSystem, varS::Symbol)

Creates an abstract array whose ``n``-th element is the value of scalar `varS`
of ``n``-th particle in `sys`.
"""
struct ParticleField <: AbstractArray{Float64, 1}
	sys::ParticleSystem
	varS::Symbol
end

Base.size(f::ParticleField) = (length(f.sys.particles),)
Base.IndexStyle(::Type{<:ParticleField}) = IndexLinear()
Base.getindex(f::ParticleField, i::Int64) = getproperty(f.sys.particles[i], f.varS)
Base.setindex!(f::ParticleField, val::Any, k::Int64) = setproperty!(f.sys.particles[k], f.varS, val)