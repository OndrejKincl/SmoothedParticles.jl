using SparseArrays

#fast scalar product and norm
function dot(x::RealVector, y::RealVector)::Float64
	return x[1]*y[1] + x[2]*y[2] + x[3]*y[3]
end

function norm(x::RealVector)::Float64
	return sqrt(x[1]*x[1] + x[2]*x[2] + x[3]*x[3])
end

function dist(p::AbstractParticle, q::AbstractParticle)::Float64
	return norm(p.x - q.x)
end

#Find vacations, where particles are missing.
function find_vacation!(a::AbstractArray)::Int64
	#find first non-missing entry
	for i in eachindex(a)
		if a[i] == 0
			return i
		end
	end
	#if a is full, expand
	i = length(a) + 1
	resize!(a, i)
	return i
end

function add_index!(cell::Cell, i::Int64)
	lock(cell.lock)
	try
		ind = find_vacation!(cell.entries)
		cell[ind] = i
		#reorder
		while ind > 1 && cell[ind-1] < cell[ind]
			temp = cell[ind]
			cell[ind] = cell[ind-1]
			cell[ind-1] = temp
			ind -= 1
		end
	finally
		unlock(cell.lock)
	end
end


"""
	create_cell_list!(sys::ParticleSystem)

Create the cell list for given particle system `sys`. This function should be
always called after updating positions. Without updated cell list, applying
binary particle operators or assembling matrices will lead to incorrect results.
"""
@inline function create_cell_list!(sys::ParticleSystem)
	
	#declare all entries null
	Threads.@threads for cell in sys.cell_list
		for k in eachindex(cell.entries)
			cell.entries[k] = 0
		end
	end
	for k in eachindex(sys.removal_cell.entries)
		sys.removal_cell.entries[k] = 0
	end

	#identify particles outside domain
	Threads.@threads for i in 1:length(sys.particles)
		p = sys.particles[i]
		if !is_inside(p, sys.domain)
			add_index!(sys.removal_cell, i)
		end
	end

	#remove particles in removal_cell
	begin
		i = 1
		while i <= length(sys.removal_cell.entries) && sys.removal_cell[i] != 0
			sys.particles[sys.removal_cell[i]] = sys.particles[end+1-i]
			i += 1
		end
		if i > 1
			resize!(sys.particles, length(sys.particles)+1-i)
		end
	end
	
	#fill the cell list with new entries
	Threads.@threads for i in 1:length(sys.particles)
		p = sys.particles[i]
		key = find_key(sys, p.x)
		cell = sys.cell_list[key]
		add_index!(cell, i)
	end
end


#Apply action between a cell and all neighbouring cells.
@inline function _apply_binary!(sys::ParticleSystem, action!::Function, p::AbstractParticle)
	key = find_key(sys, p.x)
	for Δkey in sys.key_diff
		neighbour_key = key + Δkey
		if 1 <= neighbour_key <= sys.key_max
			for j in sys.cell_list[neighbour_key].entries
				if j == 0
					break
				end
				q = sys.particles[j]
				r = dist(p,q)
				if (r > sys.h) || (p == q)
					continue
				end
				action!(p,q,r)
			end
		end
	end
end

"""
	apply_binary!(sys::ParticleSystem, action!::Function)

Apply a binary operator `action!(p::T, q::T, r::Float64)` between any two
neighbouring particles `p`, `q` in `sys::ParticleSystem{T}`. Value `r` is their
mutual distance. This excludes particle pairs with distance greater than `sys.h`.
This has linear complexity in number of particles and runs in parallel.

!!! warning "Warning"
    Modifying second particle `q` within `action!` can lead to race condition, so do not do this. Also, make sure that result will not depend on the order of particle evaluation, which is implementation-specific.
"""
function apply_binary!(sys::ParticleSystem, action!::Function)
	Threads.@threads for p in sys.particles
		_apply_binary!(sys, action!, p)
	end
end

"""
	apply_unary!(sys::ParticleSystem, action!::Function)

Apply a unary operator `action!(p::T)` on every particle `p` in
`sys::ParticleSystem{T}`. 
This has linear complexity in number of particles and runs in parallel.
"""
function apply_unary!(sys::ParticleSystem, action!::Function)
	Threads.@threads for p in sys.particles
		action!(p)
	end
end

"""
	apply!(sys::ParticleSystem, action!::Function; self::Bool = false)

Calls either [`apply_unary!`](@ref) or [`apply_binary!`](@ref) according to the
signature of `action!`. If `self == true`, then particle self-interaction
for binary operator is allowed.
"""
function apply!(sys::ParticleSystem, action!::Function; self::Bool = false)
	Type = get_particle_type(sys)
	if hasmethod(action!, (Type, Type, Float64))
		apply_binary!(sys, action!)
		if self
			apply_unary!(sys, (p -> action!(p, p, 0.0)))
		end
	else
		apply_unary!(sys, action!)
	end
end


"""
	assemble_vector(sys::ParticleSystem, func::Function)::Vector{Float64}

For given function `func(q::T)::Float64 where T <: AbstractParticle`, assemble a vector ``\\mathbf{v}``, such that 

```math
	v_i =  \\text{func}(p_i),
```

where ``p_i`` is the i-th particle in `sys::ParticleSystem{T}`.
"""
function assemble_vector(sys::ParticleSystem, func::Function)::Vector{Float64}
	N = length(sys.particles)
	v = zeros(N)
	Threads.@threads for index in 1:N
		v[index] = func(sys.particles[index])
	end
	return v
end

"""
	assemble_matrix(sys::ParticleSystem, func::Function)::SparseMatrixCSC{Float64}

For given function `func(p::T, q::T)::Float64` where `T <: AbstractParticle`, assemble a sparse matrix ``\\mathbb{A}``, such that

```math
	A_{ij} = \\text{func}(p_i, p_j, r_{ij}),
```

where ``p_i``, ``p_j`` are respectively the i-th and j-th particle in `sys::ParticleSystem{T}` and ``r_{ij}`` is their mutual distance.
This assumes that ``A_{ij} = 0`` for ``r_{ij} > h``.
"""
@inline function assemble_matrix(sys::ParticleSystem, func::Function)::SparseMatrixCSC{Float64}
	N = length(sys.particles)
	id = Dict(sys.particles[i] => i for i in 1:N)	#identify every particle by its order
	I = Int64[] 	#coord i of every nonzero entry
	J = Int64[] 	#coord j of every nonzero entry
	V = Float64[]   #value of entry i,j
	for p in sys.particles
		key = find_key(sys, p.x)
		for Δkey in sys.key_diff
			neighbour_key = key + Δkey
			if 1 <= neighbour_key <= sys.key_max
				for l in sys.cell_list[neighbour_key].entries
					if l == 0
						break
					end
					q = sys.particles[l]
					#filter remote particles
					r = dist(p,q)
					if (r > sys.h)
						continue
					end
					push!(I, id[p])
					push!(J, id[q])
					push!(V, func(p,q,r))
				end
			end
		end
	end
	return sparse(I, J, V, N, N)
end


"""
	sum(sys::ParticleSystem, func::Function, x::RealVector)::Float64

For given function `func(p::T, r::Float64)::Float64 where T <: AbstractParticle` it returns the sum

```math
	\\sum_{p \\in \\text{sys.particles}} \\text{func}(p, \\sqrt{(p.x - x)^2 + (p.y - y)^2}).
```

This can be useful if one needs to compute SPH interpolation at a point which is
not occupied by a particle.
"""
@inline function sum(sys::ParticleSystem, func::Function, x::RealVector)::Float64
	out = 0.0
	key = find_key(sys, x)
	for Δkey in sys.key_diff
		neighbour_key = key + Δkey
		if 1 <= neighbour_key <= sys.key_max
			for l in sys.cell_list[neighbour_key].entries
				if l == 0
					break
				end
				q = sys.particles[l]
				r = norm(x - q.x)
				if (r > sys.h)
					continue
				end
				out += func(q,r)
			end
		end
	end
	return out
end

"""
	sum(sys::ParticleSystem, func::Function, x::RealVector)::Float64

For given function `func(p::T, q::T, r::Float64)::Float64 where T <: AbstractParticle` and particle `p` it returns the sum

```math
	\\sum_{q \\in \\text{sys.particles}} \\text{func}(p, q, r).
```
"""
@inline function sum(sys::ParticleSystem, func::Function, p::AbstractParticle)::Float64
	out = 0.0
	key = find_key(sys, p.x)
	for Δkey in sys.key_diff
		neighbour_key = key + Δkey
		if 1 <= neighbour_key <= sys.key_max
			for l in sys.cell_list[neighbour_key].entries
				if l == 0
					break
				end
				q = sys.particles[l]
				r = norm(p.x - q.x)
				if (r > sys.h)
					continue
				end
				out += func(p,q,r)
			end
		end
	end
	return out
end