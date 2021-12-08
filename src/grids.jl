using Match
using Random
abstract type Grid end

"""
    Grid(dr::Float64, symm::Symbol)::Grid

Define a grid with a given characterstic length ``\\text{d}r`` and symmetry type
`symm`. Characterstic length means that one particle occupies a volume
``\\text{d}r^n``, where ``n`` is the dimension. Supported symmetry types are

2D:
1) `:square` (square grid)
2) `:hexagonal` (hexagrid, dual of isometric grid)
3) `:noisy` (square grid with random, divergence-free noise)

3D: (to be done in a future version)
"""
function Grid(dr::Float64, symm::Symbol)::Grid
    return @match symm begin
		:square 	=> Squaregrid(dr)
		:hexagonal  => Hexagrid(dr)
        :noisy      => NoisyGrid(dr)
		_ 			=> @error("Unsupported grid type: "*string(symm))
    end
end

mutable struct Squaregrid <: Grid
    dr::Float64
end

function covering(grid::Squaregrid, s::Shape)::Vector{Vec2}
    xs = Vec2[]
    rect = boundarybox(s)
    i_min = Int64(floor(rect.x1_min/grid.dr))
    j_min = Int64(floor(rect.x2_min/grid.dr))
    i_max = Int64(ceil(rect.x1_max/grid.dr))
    j_max = Int64(ceil(rect.x2_max/grid.dr))
	for i in i_min:i_max, j in j_min:j_max
        x1 = i*grid.dr
        x2 = j*grid.dr
        if is_inside(x1, x2, s)
            push!(xs, Vec2(x1, x2))
        end
	end
    return xs
end

mutable struct Hexagrid <: Grid
    dr::Float64
    a::Float64
    b::Float64
    Hexagrid(dr::Float64) = new(dr, (4/3)^(1/4)*dr, (3/4)^(1/4)*dr)
end

function covering(grid::Hexagrid, s::Shape)::Vector{Vec2}
    xs = Vec2[]
    rect = boundarybox(s)
    i_min = Int64(floor(rect.x1_min/grid.a)) - 1
    j_min = Int64(floor(rect.x2_min/grid.b))
    i_max = Int64(ceil(rect.x1_max/grid.a))
    j_max = Int64(ceil(rect.x2_max/grid.b))
	for i in i_min:i_max, j in j_min:j_max
        x1 = (i + (j%2)/2)*grid.a
        x2 = j*grid.b
        if is_inside(x1,x2,s)
            push!(xs, Vec2(x1,x2))
        end
	end
	return xs
end

mutable struct NoisyPoint
    x::Vec2
    w::Float64
end

mutable struct NoisyGrid <: Grid
    points::Dict{NTuple{2,Int64},NoisyPoint}
    max::Int64
    dr::Float64
    noise_fun::Function
    noise_rad::Int64
    rng::MersenneTwister
    NoisyGrid(dr::Float64) = begin
        max = -1
        noise_fun = (::Vec2 -> 1.0)
        noise_rad = 3
        rng = MersenneTwister(42)
        points = Dict{NTuple{2,Int64},NoisyPoint}()
        return new(points,max,dr,noise_fun,noise_rad,rng)
    end
end

function add_point!(ng::NoisyGrid, i::Int64, j::Int64)
    if !haskey(ng.points, (i,j))
        x = Vec2(ng.dr*i, ng.dr*j)
        np = NoisyPoint(x, (ng.dr^2)*ng.noise_fun(x)*(2*rand(ng.rng)-1.0))
        push!(ng.points, (i,j) => np)
        for di in -ng.noise_rad:ng.noise_rad, dj in -ng.noise_rad:ng.noise_rad
            i0 = i - di
            j0 = j - dj
            if haskey(ng.points, (i0,j0))
                np0 = ng.points[(i0,j0)]
                r = ng.dr*sqrt(di^2 + dj^2)
                h = ng.noise_rad*ng.dr
                F = ng.dr^3*rDwendland2(h, r)*Vec2(-dj,di)
                np.x += np0.w*F
                np0.x -= np.w*F
            end
        end
    end
end

function add_layer!(ng::NoisyGrid)
    ng.max += 1
    for k in -ng.max:ng.max
        add_point!(ng, ng.max,k)
        add_point!(ng,-ng.max,k)
        add_point!(ng,k, ng.max)
        add_point!(ng,k,-ng.max)
    end
end

function covering(ng::NoisyGrid, s::Shape)::Vector{Vec2}
    bb = boundarybox(s)
    _max = Int64(ceil((maximum([abs(bb.x1_min), abs(bb.x1_max),
        abs(bb.x2_min), abs(bb.x2_max)])/ng.dr))) + ng.noise_rad
    while ng.max < _max
        add_layer!(ng)
    end
    xs = Vec2[]
    for (_,p) in ng.points
        if is_inside(p.x, s)
            push!(xs, p.x)
        end
    end
    return xs
end 