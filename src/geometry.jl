using Match

"""
    default_symmetry = "hexagonal"
"""
const default_symmetry = "hexagonal"

"""
    Shape

Supertype for 2d geometrical shapes. Every user-defined subtype `T <: Shape` is
expected to support method

    is_inside(x::Float64, y::Float64, s::Shape) -> Bool

which determines whether a given point `(x, y)` lies inside `s`.
"""
abstract type Shape end


"""
    Rectangle(xrange::Tuple{Float64, Float64}, yrange::Tuple{Float64, Float64}) <: Shape

Define a rectangle as a cartesian product of two intervals.
"""
struct Rectangle <: Shape
    xmin::Float64
    ymin::Float64
    xmax::Float64
    ymax::Float64
    Rectangle(xrange::Tuple{Float64, Float64}, yrange::Tuple{Float64, Float64}) = begin
        (x1, x2) = xrange
        (y1, y2) = yrange
        return new(min(x1, x2), min(y1, y2), max(x1, x2), max(y1, y2))
    end
end

function is_inside(x::Float64, y::Float64, r::Rectangle)::Bool
    return (r.xmin <= x <= r.xmax) && (r.ymin <= y <= r.ymax)
end

"""
    Circle(x::Float64, y::Float64, r::Float64) <: Shape

Define a circle by specifying the center `(x, y)` and the radius `r`.
"""
struct Circle <: Shape
    x::Float64
    y::Float64
    r::Float64
    Circle(x::Float64, y::Float64, r::Float64) = begin
        @assert(r > 0.0, "radius of a circle must be positive")
        return new(x,y,r)
    end
end

function is_inside(x::Float64, y::Float64, c::Circle)::Bool
    return (x - c.x)^2 + (y - c.y)^2 <= c.r^2
end

"""
    Ellipse(x::Float64, y::Float64, r1::Float64, r2::Float64) <: Shape

Define an ellipse by specifying the center `(x, y)` and two semi-major/minor
axes `r1`, `r2`.
"""
struct Ellipse <: Shape
    x::Float64
    y::Float64
    r1::Float64
    r2::Float64
    Ellipse(x::Float64, y::Float64, r1::Float64, r2::Float64) = begin
        @assert(r1 > 0.0 && r2 > 0.0, "semi-major/minor axis must be positive")
        return new(x, y, r1, r2)
    end
end

function is_inside(x::Float64, y::Float64, e::Ellipse)::Bool
    return ((x - e.x)/e.r1)^2 + ((y - e.y)/e.r2)^2 <= 1
end

"""
    Polygon(v::Tuple{Float64, Float64}...) <: Shape

Define a polygon by specifying all vortices.
"""
struct Polygon <: Shape
    x::Vector{Float64}
    y::Vector{Float64}
    degree::Int64
    Polygon(v::Tuple{Float64, Float64}...) = begin
        degree = length(v)
        @assert(degree > 2, "number of polygon vortices must be > 2")
        x = Vector{Float64}(undef, degree)
        y = Vector{Float64}(undef, degree)
        for i in 1:degree
            x[i] = v[i][1]
            y[i] = v[i][2]
        end
        return new(x, y, degree)
    end
end


function is_inside(x::Float64, y::Float64, p::Polygon)::Bool
    #Warning! Risky math ahead!
    #Computes twice the winding number of any simple polygon.
    #Interior: +2 (if counter-clockwise) or -2 (if clockwise)
    #Edges:    +1 (if counter-clockwise) or -1 (if clockwise)
    #Vortices:  NaN
    #Exterior:  0
    #Relies on round(x) rounding half-integers to nearest even number
    winding2 = 0
    for k in 2:(p.degree - 1)
        winding2 += Int64(round((
          sign((p.x[k+1] - p.x[k  ])*(y - p.y[k  ]) - (p.y[k+1] - p.y[k  ])*(x - p.x[k  ]))
        + sign((p.x[1  ] - p.x[k+1])*(y - p.y[k+1]) - (p.y[1  ] - p.y[k+1])*(x - p.x[k+1]))
        + sign((p.x[k  ] - p.x[1  ])*(y - p.y[1  ]) - (p.y[k  ] - p.y[1  ])*(x - p.x[1  ]))
        )/2))
    end
    for k in 1:p.degree
        if p.x[k] == x && p.y[k] == y
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

function is_inside(x::Float64, y::Float64, u::BooleanUnion)::Bool
    return is_inside(x, y, u.s1)||is_inside(x, y, u.s2)
end

"""
    BooleanIntersection(s1::Shape, s2::Shape) <: Shape

Define a shape as intersection of two shapes. Equivalent to `s1 * s2`.
"""
struct BooleanIntersection <: Shape
    s1::Shape
    s2::Shape
end

function is_inside(x::Float64, y::Float64, i::BooleanIntersection)::Bool
    return is_inside(x, y, i.s1) && is_inside(x, y, i.s2)
end

"""
    BooleanDifference(s1::Shape, s2::Shape) <: Shape

Define a shape as difference of two shapes. Equivalent to `s1 - s2`.
"""
struct BooleanDifference <: Shape
    s1::Shape
    s2::Shape
end

function is_inside(x::Float64, y::Float64, d::BooleanDifference)::Bool
    return is_inside(x, y, d.s1) && !is_inside(x, y, d.s2)
end


"""
    Specification(s::Shape, f::Function) <: Shape

Define a shape of all `(x, y)` in `s`, such that `f(x,y) == true`.
"""
struct Specification <: Shape
    s::Shape
    f::Function
end

function is_inside(x::Float64, y::Float64, sp::Specification)::Bool
    return sp.f(x,y) && is_inside(x, y, sp.s)
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
    BoundaryLayer(s::Shape, dr::Float64, width::Float64; symmetry = default_symmetry) = begin
        @assert(dr > 0., "second argument must be positive")
        @assert(width > 0., "third argument must be positive")
        return new(s, dr, width, symmetry)
    end
end

function is_inside(x::Float64, y::Float64, bl::BoundaryLayer)
    if is_inside(x, y, bl.s)
        return false
    end
    if bl.symmetry == "square"
        N = Int64(ceil(bl.width/bl.dr))
        for i in -N:N, j in -N:N
            dx = i*bl.dr
            dy = j*bl.dr
            if is_inside(x + dx, y + dy, bl.s) && (dx^2 + dy^2 <= bl.width^2)
                return true
            end
        end
    elseif bl.symmetry == "hexagonal"
        N = Int64(ceil(sqrt(2)*bl.width/bl.dr))
        a = (4/3)^(1/4)*bl.dr
        b = (3/4)^(1/4)*bl.dr
        for i in -N:N, j in -N:N
            dx = (i + 0.5*j)*a
            dy = j*b
            if is_inside(x + dx, y + dy, bl.s) && (dx^2 + dy^2 <= bl.width^2)
                return true
            end
        end
    end
    return false
end

#define +,-,* on shapes as equivalent to BooleanUnion, BooleanDifference and BooleanIntersection respectively
Base.:+(s1::Shape, s2::Shape) = BooleanUnion(s1, s2)
Base.:-(s1::Shape, s2::Shape) = BooleanDifference(s1, s2)
Base.:*(s1::Shape, s2::Shape) = BooleanIntersection(s1, s2)

#Function that creates regular grid inside a given shape

function grid(sys::ParticleSystem, char_function::Function; symmetry = default_symmetry)::Tuple{Vector{Float64}, Vector{Float64}}
    return @match symmetry begin
		"square" 	=> squaregrid(sys, char_function)
		"hexagonal" => hexagrid(sys, char_function)
		_ 			=> error("invalid symmetry type: "*symmetry)
    end
end

function squaregrid(sys::ParticleSystem, char_function::Function)::Tuple{Vector{Float64}, Vector{Float64}}
	N = Int64(round((sys.xmax - sys.xmin)/sys.dr))
	M = Int64(round((sys.ymax - sys.ymin)/sys.dr))
    x = Float64[]
    y = Float64[]
	for j in 0:M, i in 0:N
        _x = sys.xmin + i*sys.dr
        _y = sys.ymin + j*sys.dr
        if char_function(_x, _y)
            push!(x, _x)
            push!(y, _y)
        end
	end
	return (x, y)
end

function hexagrid(sys::ParticleSystem, char_function::Function)::Tuple{Vector{Float64}, Vector{Float64}}
	a = (4/3)^(1/4)*sys.dr
	b = (3/4)^(1/4)*sys.dr
	N = Int64(round((sys.xmax - sys.xmin)/a))
	M = Int64(round((sys.ymax - sys.ymin)/b))
    x = Float64[]
    y = Float64[]
	for j in 0:M, i in 0:N
        _x = sys.xmin + (i + (j%2)/2)*a
        _y = sys.ymin + j*b
        if char_function(_x, _y)
            push!(x, _x)
            push!(y, _y)
        end
	end
	return (x,y)
end

"""
    generate_particles!(sys::ParticleSystem,
                        char_function::Function,
                        constructor::Function;
                        dr::Float64 = sys.dr,
                        symmetry::String = default_symmetry)

Create particles using `constructor(x::Float64, y::Float64)::AbstractParticle` at every point,
where `char_function(x::Float64, y::Float64)::Bool` is `true`. The density of particles
will be ``\\frac{1}{\\text{d}r^2}``.

Supported values of `symmetry` are `"hexagonal"` or `"square"`.
"""
function generate_particles!(sys::ParticleSystem, char_function::Function, constructor::Function; dr::Float64 = sys.dr, symmetry::String = default_symmetry)
    (x, y) = grid(sys, char_function; symmetry = symmetry)
    for k in 1:min(length(x), length(y))
        push!(sys.particles, constructor(x[k], y[k]))
    end
end

"""
    generate_particles!(sys::ParticleSystem,
                        geometry::Shape,
                        constructor::Function;
                        dr::Float64 = sys.dr,
                        symmetry::String = default_symmetry)

Create particles using `constructor(x::Float64, y::Float64)::AbstractParticle` at every point
inside a given shape. The density of particles will be ``\\frac{1}{\\text{d}r^2}``.

Supported values of `symmetry` are `"hexagonal"` or `"square"`.
"""
function generate_particles!(sys::ParticleSystem, geometry::Shape, constructor::Function; dr::Float64 = sys.dr, symmetry::String = default_symmetry)
    generate_particles!(sys, (x,y) -> is_inside(x, y, geometry), constructor; dr = dr, symmetry = symmetry)
end

"""
    snap_to_grid(sys::ParticleSystem,
                 x::Float64, y::Float64;
                 symmetry::String = default_symmetry
                 )::Tuple{Float64, Float64}

Find the nearest point to `(x, y)` which will be inside a grid.
"""
function snap_to_grid(sys::ParticleSystem,
                      x::Float64, y::Float64;
                      symmetry::String = default_symmetry
                      )::Tuple{Float64, Float64}
    if symmetry == "square"
        x = sys.xmin + round((x - sys.xmin)/sys.dr)*sys.dr
        y = sys.ymin + round((y - sys.ymin)/sys.dr)*sys.dr
    elseif symmetry == "hexagonal"
        a = (4/3)^(1/4)*sys.dr
    	b = (3/4)^(1/4)*sys.dr
    	y = sys.ymin + round((y - sys.ymin)/b)*b
        x = sys.xmin + (j%2)/2*a + round((x - sys.xmin - (j%2)/2*a)/b)*b
    end
    return (x,y)
end


#tools to remove particles from system

#remove the last particle from system
@inline function pop!(sys::ParticleSystem)
	p = sys.particles[end]
	# remove p from cell list
	key = find_key(sys, p)
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
