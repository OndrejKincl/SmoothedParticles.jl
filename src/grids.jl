using Match
abstract type Grid end

function Grid(dr::Float64, symm::Symbol)::Grid
    return @match symm begin
		:square 	=> Squaregrid(dr)
		:hexagonal  => Hexagrid(dr)
		_ 			=> error("Invalid Grid definition.")
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