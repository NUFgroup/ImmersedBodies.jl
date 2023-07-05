module TestCurves

using ImmersedBoundaryProjection.Curves
using Test

@testset "curves" begin
    # Define a simple parameterized curve function
    f(t) = (cos(2π * t), sin(2π * t))

    # Create a ParameterizedCurve instance
    curve = ParameterizedCurve(f)

    # Test partition function with 4 points for an open curve
    segments = partition(curve, 4)
    @test length(segments.points) == 4
    @test length(segments.lengths) == 4

    expected_points = [
        SVector(cos(0), sin(0)),
        SVector(cos(2π / 3), sin(2π / 3)),
        SVector(cos(4π / 3), sin(4π / 3)),
        SVector(cos(2π), sin(2π)),
    ]
    @test all(isapprox.(segments.points, expected_points, atol=1e-7))

    # Test partition function with 4 points for a closed curve
    curve_closed = ParameterizedCurve(curve.point, curve.arclength, true)
    segments_closed = partition(curve_closed, 4)
    @test length(segments_closed.points) == 4
    @test length(segments_closed.lengths) == 4

    expected_points_closed = [
        SVector(cos(0), sin(0)),
        SVector(cos(π / 2), sin(π / 2)),
        SVector(cos(π), sin(π)),
        SVector(cos(3π / 2), sin(3π / 2)),
    ]
    @test all(isapprox.(segments_closed.points, expected_points_closed, atol=1e-7))
end

end # module
