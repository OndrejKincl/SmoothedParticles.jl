using Match
using Interpolations

function is_inside(p::AbstractParticle, s::Shape)::Bool
    return is_inside(p.x, s)
end


"""
    Box(x1_min::Float64, x2_min::Float64, x3_min::Float64,
        x1_max::Float64, x2_max::Float64, x3_max::Float64)::Shape
    
Defines a box specified by two corners.
"""
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
            @error("Degenerate circle definition (r <= 0)!")
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
            @error("Degenerate ellipse definition (r1 <= 0)!")
        end
        if r2 <= 0.
            @error("Degenerate ellipse definition (r2 <= 0)!")
        end
        return new(x1, x2, r1, r2)
    end
end

function is_inside(x::RealVector, e::Ellipse)::Bool
    return ((x[1] - e.x1)/e.r1)^2 + ((x[2] - e.x2)/e.r2)^2 <= 1
end

function boundarybox(e::Ellipse)::Box
    return Rectangle(e.x1-e.r1, e.x2-e.r2, e.x1+e.r1, e.x2+e.r2) 
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
    dxs::Vector{RealVector}
    width::Float64
    BoundaryLayer(s::Shape, grid::Grid, width::Float64) = begin
        dxs = covering(grid, Ball(0.,0.,0.,width))
        return new(s, dimension(grid), dxs, width)
    end
end

function is_inside(x::RealVector, bl::BoundaryLayer)::Bool
    if is_inside(x, bl.s)
        return false
    end
    for dx in bl.dxs
        if is_inside(x + dx, bl.s)
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
end

function is_inside(x::RealVector, b::Ball)::Bool
    return (x[1] - b.x1)^2 + (x[2] - b.x2)^2  + (x[3] - b.x3)^2 <= b.r^2
end

function boundarybox(b::Ball)::Box
    return Box(b.x1-b.r, b.x2-b.r, b.x3-b.r, b.x1+b.r, b.x2+b.r, b.x3+b.r) 
end

"""
    Ellipsoid(x1::Float64, x2::Float64, x3::Float64, r1::Float64, r2::Float64, r3::Float64) <: Shape

Define an ellipsoid by specifying the center `(x1, x2, x3)` and three radii `r1,r2,r3`.
"""
struct Ellipsoid <: Shape
    x1::Float64
    x2::Float64
    x3::Float64
    r1::Float64
    r2::Float64
    r3::Float64
end

function is_inside(x::RealVector, ell::Ellipsoid)::Bool
    return ((x[1] - ell.x1)/ell.r1)^2 + ((x[2] - ell.x2)/ell.r2)^2  + ((x[3] - ell.x3)/ell.r3)^2 <= 1.0
end

function boundarybox(ell::Ellipsoid)::Box
    return Box(ell.x1-ell.r1, ell.x2-ell.r2, ell.x3-ell.r3, ell.x1+ell.r1, ell.x2+ell.r2, ell.x3+ell.r3) 
end

"""
    Transform(s::Shape; A::RealMatrix = MAT1, b::RealVector = VEC0)

Define shape as a linear transform ``x \\to Ax + b`` applied to shape `s`.
"""
struct Transform <: Shape
    s::Shape
    A::RealMatrix
    A_inv::RealMatrix
    b::RealVector
    Transform(s::Shape; A::RealMatrix = MAT1, b::RealVector = VEC0) = new(s, A, inv(A), b)
end

function is_inside(x::RealVector, tr::Transform)::Bool
    return is_inside(tr.A_inv*(x - tr.b), tr.s)
end

function boundarybox(tr::Transform)::Box
    box = boundarybox(tr.s)
    x1 = (box.x1_min, box.x1_max)
    x2 = (box.x2_min, box.x2_max)
    x3 = (box.x3_min, box.x3_max)
    x_min = [+Inf, +Inf, +Inf]
    x_max = [-Inf, -Inf, -Inf]
    for i1 in 1:2, i2 in 1:2, i3 in 1:2
        x = tr.A*RealVector(x1[i1], x2[i2], x3[i3]) + tr.b
        x_min .= min.(x_min, x)
        x_max .= max.(x_max, x)
    end
    return Box(x_min[1], x_min[2], x_min[3], x_max[1], x_max[2], x_max[3])
end

"""
    Polygon(x::Tuple{Float64, Float64}...)
"""
struct Polygon <: Shape
    xs::Vector{Float64}
    ys::Vector{Float64}
    deg::Int64
    Polygon(x::Tuple{Float64, Float64}...) = begin
        deg = length(x)
        xs = [x[k][1] for k in 1:deg]
        ys = [x[k][2] for k in 1:deg]
        return new(xs, ys, deg)
    end
end

function boundarybox(p::Polygon)::Box
    x_min = minimum(p.xs)
    x_max = maximum(p.xs)
    y_min = minimum(p.ys)
    y_max = maximum(p.ys)
    return Rectangle(x_min, y_min, x_max, y_max)
end

function is_inside(x::RealVector, p::Polygon)::Bool
    wn = 0
    x_ = x[1]
    y_ = x[2]
    for i in 1:p.deg
        next = i%p.deg + 1
        isleft = (
              (p.xs[next] - p.xs[i])*(y_ - p.ys[i])
            - (x_ - p.xs[i])*(p.ys[next] - p.ys[i])
        )
        if (p.ys[i] <= y_ < p.ys[next]) && (isleft > 0.)
            wn += 1
        end
        if (p.ys[i] > y_ >= p.ys[next]) && (isleft < 0.)
            wn -= 1
        end
    end
    return wn != 0
end

"""
    ClosedSpline(x::Tuple{Float64, Float64}...; n::Int64 = 32)
"""
function ClosedSpline(x::Tuple{Float64, Float64}...; n::Int64 = 32)::Shape
    xs = [x[k][1] for k in 1:length(x)]
    ys = [x[k][2] for k in 1:length(x)]
    push!(xs, xs[1])
    push!(ys, ys[1])
    ts = 0. : 1/length(x) : 1.0 
    itp = Interpolations.scale(interpolate(hcat(xs,ys), (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), ts, 1:2)
    ts_fine = [i/(n-1) for i in 0:(n-1)]
    y = Tuple((itp(t, 1), itp(t,2)) for t in ts_fine)
    return Polygon(y...)
end

"""
    Cone(a1::Float64, a2::Float64, a3::Float64, b1::Float64, b2::Float64, b3::Float64, ar::Float64, br::Float64)

Makes a (truncated) cone with a basis of radius `ar` centered at point `a` and tip centered at point `b` of radius `br`.
"""

struct Cone <: Shape
    a::RealVector
    b::RealVector
    ar::Float64
    br::Float64
    len::Float64
    Cone(a1::Float64, a2::Float64, a3::Float64, b1::Float64, b2::Float64, b3::Float64, ar::Float64, br::Float64) = begin
        a = RealVector(a1,a2,a3)
        b = RealVector(b1,b2,b3)
        return new(a, b, ar, br, norm(a-b))
    end        
end

function is_inside(x::RealVector, cone::Cone)::Bool
    s = dot(x - cone.a, cone.b - cone.a)
    if !(0.0 <= s <= cone.len)
        return false
    end
    t = norm(x - s*cone.b - (1-s)*cone.a)
    return (s/cone.len*cone.br + (1.0 - s/cone.len)*cone.ar >= t) 
end

function boundarybox(cone::Cone)::Box
    R = max(cone.ar, cone.br)
    x_min = min(cone.a[1], cone.b[1]) - R
    x_max = max(cone.a[1], cone.b[1]) + R
    y_min = min(cone.a[2], cone.b[2]) - R
    y_max = max(cone.a[2], cone.b[2]) + R
    z_min = min(cone.a[3], cone.b[3]) - R
    z_max = max(cone.a[3], cone.b[3]) + R
    return Box(x_min, y_min, z_min, x_max, y_max, z_max)
end


"""
    RevolutionBody(s::Shape)

Makes a 3d body by rotating a 2d shape 's' around the z-axis.
"""

struct RevolutionBody <: Shape
    s::Shape     
end

function is_inside(x::RealVector, rev::RevolutionBody)::Bool
    r = sqrt(x[1]^2 + x[2]^2)
    return is_inside(RealVector(r,x[3],0.), rev.s)
end

function boundarybox(rev::RevolutionBody)::Box
    rect = boundarybox(rev.s)
    R = rect.x1_max
    z_min = rect.x2_min
    z_max = rect.x2_max
    return Box(-R, -R, z_min, R, R, z_max)
end