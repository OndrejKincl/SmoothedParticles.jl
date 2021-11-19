using Match

"""
    Rectangle() <: Shape

Define a rectangle as a cartesian product of two intervals.
"""
struct Rectangle <: Shape
    x1_min::Float64
    x2_min::Float64
    x1_max::Float64
    x2_max::Float64
end

function is_inside(x1::Float64, x2::Float64, r::Rectangle)::Bool
    return (r.x1_min <= x1 <= r.x1_max) && (r.x2_min <= x2 <= r.x2_max)
end

function boundarybox(r::Rectangle)::Rectangle
    return r
end

"""
    Circle(x1::Float64, x2::Float64, r::Float64) <: Shape

Define a circle by specifying the center `(x1, x2)` and the radius `r`.
"""
struct Circle <: Shape
    x1::Float64
    x2::Float64
    r::Float64
    Circle(x1::Float64, x2::Float64, r::Float64) = begin
        @assert(r > 0.0, "radius of a circle must be positive")
        return new(x1,x2,r)
    end
end

function is_inside(x1::Float64, x2::Float64, c::Circle)::Bool
    return (x1 - c.x1)^2 + (x2 - c.x2)^2 <= c.r^2
end

function boundarybox(c::Circle)::Rectangle
    return Rectangle(c.x1-r, c.x2-r, c.x1+r, c.x2+r) 
end

"""
    Ellipse(x1::Float64, x2::Float64, r1::Float64, r2::Float64) <: Shape

Define an ellipse by specifying the center `(x1, x2)` and semi-major/minor
axes `r1`, `r2`.
"""
struct Ellipse <: Shape
    x1::Float64
    x2::Float64
    r1::Float64
    r2::Float64
    Ellipse(x1::Float64, x2::Float64, r1::Float64, r2::Float64) = begin
        @assert(r1 > 0.0 && r2 > 0.0, "semi-major & semi-minor axes must be positive")
        return new(x1, x2, r1, r2)
    end
end

function is_inside(x1::Float64, x2::Float64, e::Ellipse)::Bool
    return ((x1 - e.x1)/e.r1)^2 + ((x2 - e.x2)/e.r2)^2 <= 1
end

function boundarybox(e::Ellipse)::Rectangle
    return Rectangle(e.x1-r1, e.x2-r2, e.x1+r1, e.x2+r2) 
end

"""
    Polygon(v::Tuple{Float64, Float64}...) <: Shape

Define a polygon by specifying all vortices.
"""
struct Polygon <: Shape
    xs::Vector{Vec2}
    Polygon(xs::Vec2...) = begin
        degree = length(v)
        @assert(degree > 2, "number of polygon vortices must be > 2")
        return new(xs)
    end
end

function boundarybox(p::Polygon)::Rectangle
    x1_min = minimum(x -> x[1], p.xs)
    x1_max = maximum(x -> x[1], p.xs)
    x2_min = minimum(x -> x[2], p.xs)
    x2_max = maximum(x -> x[2], p.xs)
    return Rectangle(x1_min, x2_min, x1_max, x2_max) 
end


function is_inside(x1::Float64, x2::Float64, p::Polygon)::Bool
    #Warning! Risky math ahead!
    #Computes twice the winding number of any simple polygon.
    #Interior: +2 (if counter-clockwise) or -2 (if clockwise)
    #Edges:    +1 (if counter-clockwise) or -1 (if clockwise)
    #Vortices:  NaN
    #Exterior:  0
    #Relies on round(x) rounding half-integers to nearest even number
    winding2 = 0
    for k in 2:(length(xs) - 1)
        winding2 += Int64(round((
          sign((p.xs[k+1][1] - p.xs[k  ][1])*(x2 - p.xs[k  ][2]) - (p.xs[k+1][2] - p.x2[k  ][2])*(x1 - p.xs[k  ][1]))
        + sign((p.xs[1  ][1] - p.xs[k+1][1])*(x2 - p.xs[k+1][2]) - (p.xs[1  ][2] - p.x2[k+1][2])*(x1 - p.xs[k+1][1]))
        + sign((p.xs[k  ][1] - p.xs[1  ][1])*(x2 - p.xs[1  ][2]) - (p.xs[k  ][2] - p.x2[1  ][2])*(x1 - p.xs[1  ][1]))
        )/2))
    end
    for k in 1:length(xs)
        if p.x1[k] == x1 && p.x2[k] == x2
            winding2 = NaN
        end
    end
    return !(winding2 == 0)
end



"""
    BooleanUnion(s1::Shape, s2::Shape) <: Shape

Define a shape as union of two shapes. Equivalent to `s1 + s2`.
"""

struct BooleanUnion <: Shape
    s1::Shape
    s2::Shape
end

function is_inside(x1::Float64, x2::Float64, u::BooleanUnion)::Bool
    return is_inside(x1, x2, u.s1)||is_inside(x1, x2, u.s2)
end

function boundarybox(u::BooleanUnion)::Rectangle
    r1 = boundarybox(u.s1)
    r2 = boundarybox(u.s2)
    x1_min = min(r1.x1_min, r2.x1_min)
    x2_min = min(r1.x2_min, r2.x2_min)
    x1_max = max(r1.x1_max, r2.x1_max)
    x2_max = max(r1.x2_max, r2.x2_max)
    return Rectangle(x1_min, x2_min, x1_max, x2_max)
end

"""
    BooleanIntersection(s1::Shape, s2::Shape) <: Shape

Define a shape as intersection of two shapes. Equivalent to `s1 * s2`.
"""
struct BooleanIntersection <: Shape
    s1::Shape
    s2::Shape
end

function is_inside(x1::Float64, x2::Float64, i::BooleanIntersection)::Bool
    return is_inside(x1, x2, i.s1) && is_inside(x1, x2, i.s2)
end

function boundarybox(i::BooleanIntersection)::Rectangle
    r1 = boundarybox(i.s1)
    r2 = boundarybox(i.s2)
    x1_min = max(r1.x1_min, r2.x1_min)
    x2_min = max(r1.x2_min, r2.x2_min)
    x1_max = min(r1.x1_max, r2.x1_max)
    x2_max = min(r1.x2_max, r2.x2_max)
    return Rectangle(x1_min, x2_min, x1_max, x2_max)
end

"""
    BooleanDifference(s1::Shape, s2::Shape) <: Shape

Define a shape as difference of two shapes. Equivalent to `s1 - s2`.
"""
struct BooleanDifference <: Shape
    s1::Shape
    s2::Shape
end

function is_inside(x1::Float64, x2::Float64, d::BooleanDifference)::Bool
    return is_inside(x1, x2, d.s1) && !is_inside(x1, x2, d.s2)
end

function boundarybox(d::BooleanDifference)::Rectangle
    return boundarybox(d.s1)
end

"""
    Specification(s::Shape, f::Function) <: Shape

Define a shape of all `(x1, x2)` in `s`, such that `f(x1,x2) == true`.
"""
struct Specification <: Shape
    s::Shape
    f::Function
end

function is_inside(x1::Float64, x2::Float64, sp::Specification)::Bool
    return sp.f(x1,x2) && is_inside(x1, x2, sp.s)
end

function boundarybox(sp::Specification)::Rectangle
    return boundarybox(sp.s1)
end

"""
    BoundaryLayer(s::Shape, dr::Float64, width::Float64; symmetry = default_symmetry) <: Shape

Creates a layer of certain `width` around shape `s`. This requires some details
about the discretization, namely `dr` (characterstic length) and `symmetry`.

Supported values of `symmetry` are `"hexagonal"` or `"square"`.
"""
struct BoundaryLayer <: Shape
    s::Shape
    dr::Float64
    width::Float64
    symmetry::String
    BoundaryLayer(s::Shape, dr::Float64, width::Float64, symmetry = "square") = begin
        @assert(dr > 0., "second argument must be positive")
        @assert(width > 0., "third argument must be positive")
        return new(s, dr, width, symmetry)
    end
end

function is_inside(x1::Float64, x2::Float64, bl::BoundaryLayer)
    if is_inside(x1, x2, bl.s)
        return false
    end
    if bl.symmetry == "square"
        N = Int64(ceil(bl.width/bl.dr))
        for i in -N:N, j in -N:N
            dx1 = i*bl.dr
            dx2 = j*bl.dr
            if is_inside(x1 + dx1, x2 + dx2, bl.s) && (dx1^2 + dx2^2 <= bl.width^2)
                return true
            end
        end
    elseif bl.symmetry == "hexagonal"
        N = Int64(ceil(sqrt(2)*bl.width/bl.dr))
        a = (4/3)^(1/4)*bl.dr
        b = (3/4)^(1/4)*bl.dr
        for i in -N:N, j in -N:N
            dx1 = (i + 0.5*j)*a
            dx2 = j*b
            if is_inside(x1 + dx1, x2 + dx2, bl.s) && (dx1^2 + dx2^2 <= bl.width^2)
                return true
            end
        end
    end
    return false
end

function boundarybox(bl::BoundaryLayer)::Rectangle
    r = boundarybox(bl.s)
    x1_min = r.x1_min - bl.width
    x2_min = r.x2_min - bl.width
    x1_max = r.x1_max + bl.width
    x2_max = r.x2_max + bl.width
    return Rectangle(x1_min, x2_min, x1_max, x2_max)
end

#define +,-,* on shapes as equivalent to BooleanUnion, BooleanDifference and BooleanIntersection respectively
Base.:+(s1::Shape, s2::Shape) = BooleanUnion(s1, s2)
Base.:-(s1::Shape, s2::Shape) = BooleanDifference(s1, s2)
Base.:*(s1::Shape, s2::Shape) = BooleanIntersection(s1, s2)

#Function that creates regular grid inside a given shape

function grid(sys::ParticleSystem, char_function::Function, dr::Float64, symmetry = "square")::Vector{Vec2}
    return @match symmetry begin
		"square" 	=> squaregrid(sys, char_function, dr)
		"hexagonal" => hexagrid(sys, char_function, dr)
		_ 			=> error("invalid symmetry type: "*symmetry)
    end
end

function squaregrid(sys::ParticleSystem, char_function::Function, dr::Float64)::Vector{Vec2}
	N = Int64(round((sys.x1_max - sys.x1_min)/dr))
	M = Int64(round((sys.x2_max - sys.x2_min)/dr))
    xs = Vec2[]
	for j in 0:M, i in 0:N
        x1 = sys.x1_min + i*dr
        x2 = sys.x2_min + j*dr
        x = (x1, x2)
        if char_function(x)
            push!(xs, x)
        end
	end
	return xs
end

function hexagrid(sys::ParticleSystem, char_function::Function, dr::Float64)::Vector{Vec2}
	a = (4/3)^(1/4)*dr
	b = (3/4)^(1/4)*dr
	N = Int64(round((sys.x1_max - sys.x1_min)/a))
	M = Int64(round((sys.x2_max - sys.x2_min)/b))
    xs = Vec2[]
	for j in 0:M, i in 0:N
        x1 = sys.x1_min + (i + (j%2)/2)*a
        x2 = sys.x2_min + j*b
        x = Vec2(x1, x2)
        if char_function(x)
            push!(xs, x)
        end
	end
	return xs
end


"""
    generate_particles!(sys::ParticleSystem,
                        geometry::Shape,
                        constructor::Function;
                        dr::Float64 = sys.dr,
                        symmetry::String = default_symmetry)

Create particles using `constructor(x::Vec2)::AbstractParticle` at every point
inside a given shape. The density of particles will be ``\\frac{1}{\\text{d}r^2}``.

Supported values of `symmetry` are `"hexagonal"` or `"square"`.
"""
function generate_particles!(sys::ParticleSystem, geometry::Shape, constructor::Function, dr::Float64, symmetry::String = "square")
    xs = grid(sys, x -> geometry.is_inside(x[1], x[2]), dr, symmetry)
    for x in xs
        push!(sys.particles, constructor(x))
    end
end


#tools to remove particles from system

#remove the last particle from system
@inline function pop!(sys::ParticleSystem)
	p = sys.particles[end]
	# remove p from cell list
	key = find_key(sys, p.x[1], p.x[2])
	if 1 <= key <= sys.key_max
		for i in 1:length(sys.cell_list[key])
			if ismissing(sys.cell_list[key][i])
				continue
			end
			if p == sys.cell_list[key][i]
				sys.cell_list[key][i] = missing
			end
		end
	end
	#resize particlel list
	Base.pop!(sys.particles)
end

#sort array by boolean value (0 or 1) in linear time
function sort_by_bool!(a::AbstractArray, f::Function)
       i = 1
       j = length(a)
       while i < j
           while i < j && !f(a[i])
               i+=1
           end
           while i < j && f(a[j])
               j-=1
           end
           temp = a[i]
           a[i] = a[j]
           a[j] = temp
           i += 1
           j -= 1
       end
end

"""
    remove_particles!(sys::ParticleSystem, criterion::Function)

Remove all particles `p` in `sys` satisfying `criterion(p) == true`.
"""
@inline function remove_particles!(sys::ParticleSystem, criterion::Function)
	#move particles satisfying criterion to the end of the list
	sort_by_bool!(sys.particles, criterion)
	#pop last particle until it no longer satisfies the criterion
	while length(sys.particles) > 0 && criterion(sys.particles[end])
		pop!(sys)
	end
end