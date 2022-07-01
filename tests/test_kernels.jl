module test_kernels
include("../src/SmoothedParticles.jl")
using .SmoothedParticles
using Test

const TOL = 0.01
const N = 1000

function simpson_rule(f::Function, a::Float64, b::Float64, n::Int64=N)
    I = 0.
    h = (b-a)/n
    for i in 1:n-1
        _a = a + i*h
        _b = a + (i+1)*h
        I += h/6.0*(f(_a) + 4.0*f(0.5*(_a + _b)) + f(_b))
    end
    return I
end

function test_local_ker(dim::Int64, f, Df, rDf)
    h = 0.42
    @test f(h, 4.0) == 0.
    @test isfinite(f(h, 0.0))
    int = 0.
    if dim == 2
        int = simpson_rule(r -> 2.0*pi*r*f(h,r), 0., h)
    elseif dim == 3
        int = simpson_rule(r -> 4.0*pi*r*r*f(h,r), 0., h)
    end
    @test int ≈ 1.0 rtol=TOL
    @test Df(h, 4.0) == 0.
    @test isfinite(Df(h, 0.0))
    int = simpson_rule(r -> Df(h,r), 0.2, 0.3)
    diff = f(h, 0.3) - f(h, 0.2)
    @test int ≈ diff rtol=0.01
    @test rDf(h, 4.0) == 0.
    @test isfinite(rDf(h, 0.0))
    @test rDf(h, 0.1) ≈ Df(h, 0.1)/0.1 rtol=TOL
end

function main()
    @testset "wendland2" begin
        test_local_ker(2, wendland2, Dwendland2, rDwendland2)
    end
    @testset "spline23" begin
        test_local_ker(2, spline23, Dspline23, rDspline23)
    end
    @testset "spline24" begin
        test_local_ker(2, spline24, Dspline24, rDspline24)
    end
    @testset "wendland3" begin
        test_local_ker(3, wendland3, Dwendland3, rDwendland3)
    end
end

end