using SparseArrays

#fast scalar product
function dot(x::Vec2, y::Vec2)::Float64
	return x[1]*y[1] + x[2]*y[2]
end

function norm(x::Vec2)::Float64
	return sqrt(x[1]*x[1] + x[2]*x[2])
end

#Find the cell index for given coordinates x1,x2.
function find_key(sys::ParticleSystem, x1::Float64, x2::Float64, Δi::Int64 = 0, Δj::Int64 = 0)::Int64
	try
		i = 1 + Int64(floor(x1/sys.h)) + Δi - sys.i_phase
		j = 1 + Int64(floor(x2/sys.h)) + Δj - sys.j_phase
		return i + sys.i_max*(j-1)
	catch
		return -1
	end
end

#Find vacations, where particles are missing.
function find_vacation!(a::AbstractArray)::Int64
	#find first non-missing entry
	for i in eachindex(a)
		if ismissing(a[i])
			return i
		end
	end
	#if a is full, expand
	i = length(a) + 1
	resize!(a, i)
	#@warn("extending particles per cell limit to "*string(i))
	return i
end

"""
	create_cell_list!(sys::ParticleSystem)

Create the cell list for given particle system `sys`. This function should be
always called after updating positions. Without updated cell list, applying
binary particle operators or assembling matrices will lead to incorrect results.
"""
@inline function create_cell_list!(sys::ParticleSystem)
	#delete all particles outside of the domain
	remove_particles!(sys, p -> !is_inside(p, sys.domain))

	#declare all entries of every cell as missing
	Threads.@threads for cell in sys.cell_list
		for k in 1:length(cell)
			cell[k] = missing
		end
	end
	#fill the cell list with new entries
	Threads.@threads for p in sys.particles
		key = find_key(sys, p.x[1], p.x[2])
		if 1 <= key <= sys.key_max
			cell = sys.cell_list[key]
			lock(sys.cell_lock[key])
			try
				ind = find_vacation!(cell)
				cell[ind] = p
			finally
				unlock(sys.cell_lock[key])
			end
		end
	end
end


#Apply action between a cell and all neighbouring cells.
@inline function _apply_binary!(sys::ParticleSystem, action!::Function, p::AbstractParticle)
	key = find_key(sys, p.x[1], p.x[2])
	for Δkey in 0:8
		neighbour_key = key + (Δkey%3-1) + sys.i_max*(div(Δkey,3)%3-1)
		if 1 <= neighbour_key <= sys.key_max
			for q in sys.cell_list[neighbour_key]
				if ismissing(q)
					break
				end
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
	if hasmethod(action!, (sys.particle_type, sys.particle_type, Float64))
		apply_binary!(sys, action!)
		if self
			apply_unary!(sys, (p -> action!(p, p, 0.0)))
		end
	else
		apply_unary!(sys, action!)
	end
end

"""
	dist(p::AbstractParticle, q::AbstractParticle)::Float64

Calculate the distance between two particles `p` and `q`.
"""
function dist(p::AbstractParticle, q::AbstractParticle)::Float64
	return sqrt((p.x[1] - q.x[1])^2 + (p.x[2] - q.x[2])^2)
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

For given function `func(p::T, q::T)::Float64 where T <: AbstractParticle`, assemble a sparse matrix ``\\mathbb{A}``, such that

```math
	A_{ij} = \\text{func}(p_i, p_j),
```

where ``p_i``, ``p_j`` are respectively the i-th and j-th particle in `sys::ParticleSystem{T}`.
"""
@inline function assemble_matrix(sys::ParticleSystem, func::Function)::SparseMatrixCSC{Float64}
	N = length(sys.particles)
	id = Dict(sys.particles[i] => i for i in 1:N)	#identify every particle by its order
	I = Int64[] 	#coord i of every nonzero entry
	J = Int64[] 	#coord j of every nonzero entry
	V = Float64[]   #value of entry i,j
	for p in sys.particles
		for k in 0:8
			key = find_key(sys, p.x[1], p.x[2], k%3-1, div(k,3)-1)
			if 1 <= key <= sys.key_max
				for q in sys.cell_list[key]
					if ismissing(q)
						break
					end
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
	sum(sys::ParticleSystem, func::Function, x::Vec2)::Float64

For given function `func(p::T, r::Float64)::Float64 where T <: AbstractParticle` it returns the sum

```math
	\\sum_{p \\in \\text{sys.particles}} \\text{func}(p, \\sqrt{(p.x - x)^2 + (p.y - y)^2}).
```

This can be useful if one needs to compute SPH interpolation at a point which is
not occupied by a particle.
"""
@inline function sum(sys::ParticleSystem, func::Function, x::Vec2)::Float64
	out = 0.0
	for k in 0:8
		key = find_key(sys, x[1], x[2], k%3-1, div(k,3)-1)
		if 1 <= key <= sys.key_max
			for q in sys.cell_list[key]
				if ismissing(q)
					break
				end
				r = sqrt((x[1] - q.x[1])^2 + (x[2] - q.x[2])^2)
				if (r > sys.h)
					continue
				end
				out += func(q,r)
			end
		end
	end
	return out
end

@inline function sum(sys::ParticleSystem, func::Function, p::AbstractParticle)::Float64
	out = 0.0
	for k in 0:8
		key = find_key(sys, x[1], x[2], k%3-1, div(k,3)-1)
		if 1 <= key <= sys.key_max
			for q in sys.cell_list[key]
				if ismissing(q)
					break
				end
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