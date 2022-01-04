using Match

function is_inside(p::AbstractParticle, s::Shape)::Bool
    return is_inside(p.x, s)
end

struct Box <: Shape
    x1_min::Float64
    x2_min::Float64
    x3_min::Float64
    x1_max::Float64
    x2_max::Float64
    x3_max::Float64
end

function is_inside(x::RealVector, b::Box)::Bool
    return (
        b.x1_min <= x[1] <= b.x1_max &&
        b.x2_min <= x[2] <= b.x2_max &&
        b.x3_min <= x[3] <= b.x3_max
    )
end

function boundarybox(b::Box)::Box
    return b    
end

"""
    Rectangle(x1_min, x2_min, x1_max, x2_max) <: Shape

Define a rectangle by specifying bottom left and top right corner.
"""
function Rectangle(x1_min::Float64, x2_min::Float64, x1_max::Float64, x2_max::Float64)::Box
    return Box(x1_min, x2_min, 0., x1_max, x2_max, 0.)
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
        if r <= 0.0
            @warn("Degenerate circle definition (r <= 0)!")
        end
        return new(x1,x2,r)
    end
end

function is_inside(x::RealVector, c::Circle)::Bool
    return (x[1] - c.x1)^2 + (x[2] - c.x2)^2 <= c.r^2
end

function boundarybox(c::Circle)::Box
    return Rectangle(c.x1-c.r, c.x2-c.r, c.x1+c.r, c.x2+c.r) 
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
        if r1 <= 0.
            @warn("Degenerate ellipse definition (r1 <= 0)!")
        end
        if r2 <= 0.
            @warn("Degenerate ellipse definition (r2 <= 0)!")
        end
        return new(x1, x2, r1, r2)
    end
end

function is_inside(x::RealVector, e::Ellipse)::Bool
    return ((x[1] - e.x1)/e.r1)^2 + ((x[2] - e.x2)/e.r2)^2 <= 1
end

function boundarybox(e::Ellipse)::Box
    return Rectangle(e.x1-r1, e.x2-r2, e.x1+r1, e.x2+r2) 
end



"""
    BooleanUnion(s1::Shape, s2::Shape) <: Shape

Define a shape as union of two shapes. Equivalent to `s1 + s2`.
"""

struct BooleanUnion <: Shape
    s1::Shape
    s2::Shape
end

function is_inside(x::RealVector, u::BooleanUnion)::Bool
    return is_inside(x, u.s1)||is_inside(x, u.s2)
end

function boundarybox(u::BooleanUnion)::Box
    r1 = boundarybox(u.s1)
    r2 = boundarybox(u.s2)
    x1_min = min(r1.x1_min, r2.x1_min)
    x2_min = min(r1.x2_min, r2.x2_min)
    x3_min = min(r1.x3_min, r2.x3_min)
    x1_max = max(r1.x1_max, r2.x1_max)
    x2_max = max(r1.x2_max, r2.x2_max)
    x3_max = max(r1.x3_max, r2.x3_max)
    return Box(x1_min, x2_min, x3_min, x1_max, x2_max, x3_max)
end

"""
    BooleanIntersection(s1::Shape, s2::Shape) <: Shape

Define a shape as intersection of two shapes. Equivalent to `s1 * s2`.
"""
struct BooleanIntersection <: Shape
    s1::Shape
    s2::Shape
end

function is_inside(x::RealVector, i::BooleanIntersection)::Bool
    return is_inside(x, i.s1) && is_inside(x, i.s2)
end

function boundarybox(i::BooleanIntersection)::Box
    r1 = boundarybox(i.s1)
    r2 = boundarybox(i.s2)
    x1_min = max(r1.x1_min, r2.x1_min)
    x2_min = max(r1.x2_min, r2.x2_min)
    x3_min = max(r1.x3_min, r2.x3_min)
    x1_max = min(r1.x1_max, r2.x1_max)
    x2_max = min(r1.x2_max, r2.x2_max)
    x3_max = min(r1.x3_max, r2.x3_max)
    return Box(x1_min, x2_min, x3_min, x1_max, x2_max, x3_max)
end

"""
    BooleanDifference(s1::Shape, s2::Shape) <: Shape

Define a shape as difference of two shapes. Equivalent to `s1 - s2`.
"""
struct BooleanDifference <: Shape
    s1::Shape
    s2::Shape
end

function is_inside(x::RealVector, d::BooleanDifference)::Bool
    return is_inside(x, d.s1) && !is_inside(x, d.s2)
end

function boundarybox(d::BooleanDifference)::Box
    return boundarybox(d.s1)
end

"""
    Specification(s::Shape, f::Function) <: Shape

Define a shape of all `x` in `s`, such that `f(x) == true`.
"""
struct Specification <: Shape
    s::Shape
    f::Function
end

function is_inside(x::RealVector, sp::Specification)::Bool
    return sp.f(x) && is_inside(x, sp.s)
end

function boundarybox(sp::Specification)::Box
    return boundarybox(sp.s)
end

"""
    BoundaryLayer(s::Shape, grid::Grid, width::Float64) <: Shape

Creates a layer of certain `width` around shape `s`. More specifically, a point is inside
boundary layer if it is not in `s` and at the same time has distance less than `width` to
at least one point on `grid` in `s`.
"""
struct BoundaryLayer <: Shape
    s::Shape
    dim::Int64
    xs::Vector{RealVector}
    width::Float64
    BoundaryLayer(s::Shape, grid::Grid, width::Float64) = begin
        if width <= 0.
            @warn("Degenerate boundary layer definition (width <= 0)!")
        end
        xs = covering(grid, s)
        return new(s, dimension(grid), xs, width)
    end
end

function is_inside(x::RealVector, bl::BoundaryLayer)::Bool
    if is_inside(x, bl.s)
        return false
    end
    for y in bl.xs
        if norm(x - y) < bl.width
            return true
        end
    end
    return false
end

function boundarybox(bl::BoundaryLayer)::Box
    r = boundarybox(bl.s)
    x1_min = r.x1_min - bl.width
    x2_min = r.x2_min - bl.width
    x3_min = r.x3_min - bl.width
    x1_max = r.x1_max + bl.width
    x2_max = r.x2_max + bl.width
    x3_max = r.x3_max + bl.width
    if bl.dim == 2
        return Rectangle(x1_min, x2_min, x1_max, x2_max)
    else
        return Box(x1_min, x2_min, x3_min, x1_max, x2_max, x3_max)
    end
end

#define +,-,* on shapes as equivalent to BooleanUnion, BooleanDifference and BooleanIntersection respectively
Base.:+(s1::Shape, s2::Shape) = BooleanUnion(s1, s2)
Base.:-(s1::Shape, s2::Shape) = BooleanDifference(s1, s2)
Base.:*(s1::Shape, s2::Shape) = BooleanIntersection(s1, s2)


"""
    Ball(x1::Float64, x2::Float64, x3::Float64, r::Float64) <: Shape

Define a ball by specifying the center `(x1, x2, x3)` and the radius `r`.
"""
struct Ball <: Shape
    x1::Float64
    x2::Float64
    x3::Float64
    r::Float64
    Ball(x1::Float64, x2::Float64, x3::Float64, r::Float64) = begin
        if r <= 0.0
            @warn("Degenerate ball definition (r <= 0)!")
        end
        return new(x1,x2,x3,r)
    end
end

function is_inside(x::RealVector, b::Ball)::Bool
    return (x[1] - b.x1)^2 + (x[2] - b.x2)^2  + (x[3] - b.x3)^2 <= b.r^2
end

function boundarybox(b::Ball)::Box
    return Box(b.x1-b.r, b.x2-b.r, b.x3-b.r, b.x1+b.r, b.x2+b.r, b.x3+b.r) 
end