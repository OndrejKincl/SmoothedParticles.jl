module test_geometry
include("../src/SmoothedParticles.jl")
using .SmoothedParticles
using Test
using Images


function rotmat(x::Float64)::RealMatrix
    return RealMatrix(
         cos(x),  sin(x), 0.,
        -sin(x),  cos(x), 0.,
             0.,      0., 1.
    )
end

# create a complicated geometry using the library
function make_figure()::Shape
    head = Circle(0., 0.8, 0.2) * Rectangle(-1.0, -1.0, 1.0, 0.9)
    left_eye = Ellipse(-0.1, 0.8, 0.03, 0.01)
    right_eye = Transform(left_eye; b = 0.2VECX)
    eyes = left_eye + right_eye
    mouth = Specification(
        Rectangle(-0.15, 0.6, 0.15, 1.0), 
        x -> (abs(x[2] - 0.72 - 0.02*sin(40*x[1])) < 0.01) 
    )
    head -= eyes + mouth
    body = Rectangle(-0.15, 0.1, 0.15, 0.6)
    left_hand = Transform(Rectangle(-0.15, 0.38, 0.1, 0.42); A=rotmat(pi/6))
    mirrormat = RealMatrix(-1.,0.,0.,  0.,1.,0.,  0.,0.,1.)
    right_hand = Transform(left_hand; A=mirrormat)
    hair = Polygon((-0.1,0.9), (+0.1,0.9), (0.1,0.98), (0.05, 0.94), (0.0, 0.98), (-0.05, 0.9))
    heart = ClosedSpline((0.02, 0.39), (0.07, 0.47), (0.04, 0.49), (0.02, 0.47), (0.0, 0.49), (-0.03, 0.47))
    figure = (body - heart) + head + left_hand + right_hand + hair
    return figure
end

function main()
    @testset "geo_2d" begin
        figure = make_figure()
        # plot figure using Images.jl
        res = 1000
        img = zeros(RGB, res, res)
        @test_nowarn for i in 1:res, j in 1:res
            x = (j-1)/(res-1) - 0.5
            y = 1.0 - (i-1)/(res-1)
            p = RealVector(x, y, 0.)
            if is_inside(p, figure)
                img[i,j] = RGB(1.0, 1.0, 1.0)
            end
        end
        # uncomment next line to create reference solution or to perform visual check
        # save("figure.png", img)
        #_img = RGB.(load(joinpath(@__DIR__, "figure.png")))
        #@test img == _img
    end
end


end
