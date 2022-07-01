using Test

@testset "kernels" begin
    include("test_kernels.jl")
    test_kernels.main()
end

@testset "geometry" begin
    include("test_geometry.jl")
    test_geometry.main()
end

@testset "IO" begin
    include("test_IO.jl")
    test_IO.main()
end

@testset "collision 2d" begin
    include("test_collision_2d.jl")
    test_collision_2d.main()
end