using  StaticArrays

const RealVector = SArray{Tuple{3},Float64,1,3}
const RealMatrix = SArray{Tuple{3,3},Float64,2,9}
const VECX = RealVector(1., 0., 0.)
const VECY = RealVector(0., 1., 0.)
const VECZ = RealVector(0., 0., 1.)
const VEC0 = zero(RealVector)
const MAT0 = zero(RealMatrix)
const MAT1 = RealMatrix(1., 0., 0., 
                        0., 1., 0.,
						0., 0., 1.)

"""
	AbstractParticle

Abstract supertype for smoothed particles. Any structure with `AbstractParticle`
supertype is expected to:

* be mutable
* have field `x::RealVector` (particle position)
"""
abstract type AbstractParticle end

"""
    Shape

Supertype for geometrical shapes.
"""
abstract type Shape end

#cell for storing particle indeces
mutable struct Cell
	entries::Vector{Int64}
	lock::Threads.ReentrantLock
	Cell() = new(Int64[], Threads.ReentrantLock())
end

Base.size(c::Cell) = (length(c.entries),)
Base.IndexStyle(::Type{<:Cell}) = IndexLinear()
Base.getindex(c::Cell, i::Int64) = getindex(c.entries, i)
Base.setindex!(c::Cell, val::Any, i::Int64) = setindex!(c.entries, val, i)

"""
	ParticleSystem(T::Type, domain::Shape, h::Float64)

Struct that contains all vital information about the simulation.
The constructor specifies that:
* the simulation will use particles of type `T <: AbstractParticle`,
* Particles outside of the 'domain' can be disregarded (and will be automatically removed).
* Particles are considered neighbours if their distance is less than `h`.
Please, do not make 'domain' unnecessarily large (has negative impact on performance).
"""
struct ParticleSystem{T <: AbstractParticle}
	h::Float64
	domain::Shape

	#values for computing location keys
	key_phase::NTuple{3, Int64}
	key_lim::NTuple{3, Int64}
	key_max::Int64
	key_diff::Vector{Int64}

	particles::Vector{T}
	cell_list::Vector{Cell}
	removal_cell::Cell

	ParticleSystem(T::DataType, domain::Shape, h::Float64) =
	begin
			@assert(h  > 0.0, "invalid ParticleSystem declaration! (h must be a positive float)")
			@assert(T <: AbstractParticle, "invalid ParticleSystem declaration! ("*string(T)*" is not an AbstractParticle subtype)")
			@assert(hasfield(T, :x) && (attribute_type(T, :x) == RealVector), "invalid ParticleSystem declaration! (particles must have a field `x::RealVector`)")

			box = boundarybox(domain)
			x_min = (box.x1_min, box.x2_min, box.x3_min)
			x_max = (box.x1_max, box.x2_max, box.x3_max)
			key_phase = Int64.(floor.(x_min./h))
			key_lim = Int64.(floor.(x_max./h)) .- key_phase .+ 1
			key_max = prod(key_lim)
			key_diff = Int64[]
			if key_lim[3] == 1
				@info("2D sim")
				#simulation in 2d
				for di in -1:1, dj in -1:1
					push!(key_diff, di + key_lim[1]*dj)
				end
			else
				@info("3D sim")
				#simulation in 3d
				for di in -1:1, dj in -1:1, dk in -1:1
					push!(key_diff, di + key_lim[1]*(dj + key_lim[2]*dk))
				end
			end
			particles = Vector{T}()
			cell_list = [Cell() for _ in 1:key_max]
			removal_cell = Cell()
			return new{T}(
				h, box,
				key_phase, key_lim, key_max, key_diff,
				particles, cell_list, removal_cell
			)
	end
end

get_particle_type(::ParticleSystem{T}) where T = T

#Find the cell index for given coordinate x::RealVector
function find_key(sys::ParticleSystem, x::RealVector)::Int64
	try
		i = 1 + Int64(floor(x[1]/sys.h)) - sys.key_phase[1]
		j = 1 + Int64(floor(x[2]/sys.h)) - sys.key_phase[2]
		k = 1 + Int64(floor(x[3]/sys.h)) - sys.key_phase[3]
		return i + sys.key_lim[1]*(j-1) + sys.key_lim[1]*sys.key_lim[2]*(k-1)
	catch #Int64(x) fails when x is nan or infinity
		return -1
	end
end

"""
	ParticleField(sys::ParticleSystem, varS::Symbol)

Creates an abstract array whose ``n``-th element is the value of scalar `varS`
of ``n``-th particle in `sys`.

!!! warning "Warning"
    The indentity of ``n``-th particle in `ParticleSystem` may change
    when particles are added or removed.
"""
struct ParticleField <: AbstractArray{Float64, 1}
	sys::ParticleSystem
	varS::Symbol
end

Base.size(f::ParticleField) = (length(f.sys.particles),)
Base.getindex(f::ParticleField, i::Int64) = getproperty(f.sys.particles[i], f.varS)
Base.setindex!(f::ParticleField, val::Any, k::Int64) = setproperty!(f.sys.particles[k], f.varS, val)


function attribute_type(type::DataType, var::Symbol)
    ind = findfirst(s -> s == var, fieldnames(type))
    if typeof(ind) == Nothing
        throw("Variable "*string(var)* 
            " does not exist!")
    end
    return fieldtypes(type)[ind]
end