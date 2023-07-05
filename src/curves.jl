module Curves

export Curve, Segments, arclength, partition

using StaticArrays

"""
    Curve

A mathematical curve in 2D space. Calling `curve(t)` with `0 ≤ t ≤ 1` gives an `[x y]` point
along the curve.
"""
abstract type Curve end
abstract type ClosedCurve <: Curve end
abstract type OpenCurve <: Curve end

"""
    isclosed(curve::Curve)

Return true if the `curve`'s start and end points are equal, false otherwise.
"""
isclosed(::ClosedCurve) = true
isclosed(::OpenCurve) = false

struct Segments
    points::Matrix{Float64}
    lengths::Vector{Float64}
end

"""
    arclength(curve::Curve)

The total arclength of a curve.
"""
function arclength end

"""
    partition(curve::Curve, ds::AbstractFloat) :: Segments
    partition(curve::Curve, n::Integer) :: Segments

Partition a curve into segments of approximately equal length. Either specify the target
segment length `ds` or the target segment count `n`. The curve's endpoints are preserved.
"""
function partition end

# Curves either need to define a method for ds or for n
# The other is satisfied by these defaults
function partition(curve::Curve, ds::AbstractFloat)
    n = round(Int, arclength(curve) / ds) + !isclosed(curve)
    return partition(curve, n)
end
function partition(curve::Curve, n::Integer)
    ds = arclength(curve) / (n - !isclosed(curve))
    return partition(curve, ds)
end

"""
    LineSegment((x1, y1), (x2, y2)) :: Curve

A line segment between two points.
"""
struct LineSegment <: OpenCurve
    p1::SVector{2,Float64}
    p2::SVector{2,Float64}
end

(line::LineSegment)(t) = line.p1 + t * (line.p2 - line.p1)

arclength(line::LineSegment) = hypot((line.p2 .- line.p1)...)

function Base.show(io::IO, ::MIME"text/plain", line::LineSegment)
    print(io, "LineSegment: from ", line.p1, " to ", line.p2)
    return nothing
end

function partition(line::LineSegment, n::Integer)
    @assert n > 1

    xs, ys = (range(x1, x2, n) for (x1, x2) in zip(line.p1, line.p2))
    ds = hypot(xs[2] - xs[1], ys[2] - ys[1])

    points = [xs ys]
    lengths = fill(ds, n)

    return Segments(points, lengths)
end

"""
    Circle(r=1, center=(0, 0)) :: Curve

A circle with radius `r` centered at `center`.
"""
struct Circle <: ClosedCurve
    r::Float64
    center::SVector{2,Float64}
    Circle(r=1, center=(0, 0)) = new(r, center)
end

function (circle::Circle)(t)
    s = 2 * pi * t
    return circle.center + circle.r * SVector(cos(s), sin(s))
end

arclength(circle::Circle) = 2 * pi * circle.r

function Base.show(io::IO, ::MIME"text/plain", circle::Circle)
    print(io, "Circle: radius=", circle.r, ", center=", circle.center)
    return nothing
end

function partition(circle::Circle, n::Integer)
    x0, y0 = circle.center
    r = circle.r

    t = 2pi / n * (0:(n - 1))

    xs = @. x0 + r * cos(t)
    ys = @. y0 + r * sin(t)
    ds = hypot(xs[2] - xs[1], ys[2] - ys[1])

    points = [xs ys]
    lengths = fill(ds, n)

    return Segments(points, lengths)
end

struct ParameterizedCurve{F} <: Curve
    point::F
    arclength::Float64
    closed::Bool
end

function ParameterizedCurve(f)
    # f(t) -> (x(t), y(t))
    arclength = 0 # get arclength
    itp = Interpolation() # define interpolation

    function point(t)
        s = itp(t)
        return SVector(f(t))
    end

    closed = true # is first point approx equal last point?

    return ParameterizedCurve(point, arclength, closed)
end

isclosed(curve::ParameterizedCurve) = curve.closed
arclength(curve::ParameterizedCurve) = curve.arclength

(curve::ParameterizedCurve)(t) = curve.point(t)

function partition(curve::ParameterizedCurve, n::Integer)
    # n is the number of points you should return

    # Initialize the points and lengths arrays
    points = Vector{SVector{2,Float64}}(undef, n)
    lengths = Vector{Float64}(undef, n)

    # Calculate the step size for parameter t based on whether the curve is closed or open
    step_size = curve.closed ? 1 / n : 1 / (n - 1)

    # Iterate through the range of t values and calculate points and lengths
    for i in 0:(n - 1)
        t = i * step_size
        points[i + 1] = curve.point(t)
        if i > 0
            lengths[i] = norm(points[i + 1] - points[i])
        end
    end

    # For a closed curve, set the length of the last segment
    if curve.closed
        lengths[end] = norm(points[end] - points[1])
    end

    return Segments(points, lengths)
end

end # module Curves
