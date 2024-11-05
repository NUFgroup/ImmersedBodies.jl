module Curves

export Curve, ClosedCurve, OpenCurve, Segments
export isclosed, arclength, partition
export TransformedCurve, Transformation, translate, rotate, scale
export LineSegment, Circle, ParameterizedCurve
export NacaParams, NACA, NACA4

using StaticArrays
using Interpolations
using FunctionWrappers: FunctionWrapper
using LinearAlgebra: dot, normalize, norm, I

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

"""
    AbstractSegments

A sequence of segments with associated [`points`](@ref) and [`lengths`](@ref).

See [`Segments`](@ref).
"""
abstract type AbstractSegments end

"""
    points(::AbstractSegments)

The matrix of `[x y]` points associated with each segment.
"""
function points end

"""
    lengths(::AbstractSegments)

The vector of lengths associated with each segment.
"""
function lengths end

"""
    Segments(points::AbstractMatrix, lengths::AbstractVector) :: AbstractSegments

A sequence of segments where the point `points[i, :]` is associated with length
`lengths[i]`.
"""
struct Segments <: AbstractSegments
    points::Matrix{Float64}
    lengths::Vector{Float64}
    function Segments(points, lengths)
        @assert axes(points, 1) == eachindex(lengths)
        return new(points, lengths)
    end
end

points(s::Segments) = s.points
lengths(s::Segments) = s.lengths

"""
    SplitSegments <: AbstractSegments

Segments created when [`partition`](@ref)ing a [`SplitCurve`](@ref). To retrieve the
individual parts, see [`parts`](@ref).
"""
struct SplitSegments{N} <: AbstractSegments
    points::Matrix{Float64}
    lengths::Vector{Float64}
    splits::NTuple{N,UnitRange{Int}}
    function SplitSegments(points, lengths, splits::NTuple{N}) where {N}
        @assert axes(points, 1) == eachindex(lengths)
        return new{N}(points, lengths, splits)
    end
end

points(s::SplitSegments) = s.points
lengths(s::SplitSegments) = s.lengths

"""
    separate(segments::SplitSegments)

Returns a tuple of [`Segments`](@ref) for each split part.
"""
function separate(s::SplitSegments)
    return map(s.splits) do r
        Segments(s.points[r, :], s.lengths[r])
    end
end

"""
    arclength(curve::Curve)

The total arclength of a curve.
"""
function arclength end

"""
    partition(curve::Curve, ds::AbstractFloat) :: AbstractSegments
    partition(curve::Curve, n::Integer) :: AbstractSegments

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
    Transformation

A 2-D transformation that is the combination of translation, rotation, and uniform scaling.

Constructed using [`translate`](@ref), [`rotate`](@ref), and [`scale`](@ref).
"""
struct Transformation
    offset::SVector{2,Float64} # Translation vector
    rotation::SMatrix{2,2,Float64,4} # Rotation matrix
    scale::Float64 # Uniform scale
end

(t::Transformation)(v) = t(SVector{2,Float64}(v))
function (t::Transformation)(v::AbstractVector{<:Real})
    return t.offset + t.rotation * t.scale * v
end
function (t1::Transformation)(t2::Transformation)
    offset = t1(t2.offset)
    rotation = t1.rotation * t2.rotation
    scale = t1.scale * t2.scale
    return Transformation(offset, rotation, scale)
end

const I2 = SMatrix{2,2,Float64}(I) # 2x2 identity
function _rot(θ::Real)
    c = cos(θ)
    s = sin(θ)
    return @SMatrix [c -s; s c]
end

"""
    translate(v) :: Transformation

A [`Transformation`](@ref) that translates by vector `v`.
"""
translate(v) = Transformation(v, I2, 1)

"""
    rotate(θ) :: Transformation

A [`Transformation`](@ref) that rotates counter-clockwise by angle `θ`.
"""
rotate(θ) = Transformation(zeros(SVector{2}), _rot(θ), 1)

"""
    scale(k) :: Transformation

A [`Transformation`](@ref) that uniformly scales by factor `k`.
"""
scale(k) = Transformation(zeros(SVector{2}), I2, k)

"""
    TransformedCurve

A curve that is transformed in space. Constructed by passing a [`Curve`](@ref) to a
[`Transformation`](@ref).
"""
struct TransformedCurve{C<:Curve} <: Curve
    base::C
    transform::Transformation
end

"""
    (t::Transformation)(curve::Curve)

Transform a curve, returning a [`TransformedCurve`](@ref).
"""
(t::Transformation)(curve::Curve) = TransformedCurve(curve, t)
(t::Transformation)(curve::TransformedCurve) = t(curve.transform)(curve.base)

isclosed(curve::TransformedCurve) = isclosed(curve.base)
arclength(curve::TransformedCurve) = abs(curve.transform.scale) * arclength(curve.base)
(curve::TransformedCurve)(t) = curve.transform(curve.base(t))

function partition(curve::TransformedCurve, n::Integer)
    segments = partition(curve.base, n)

    points = Curves.points(segments)
    lengths = Curves.lengths(segments)
    for i in eachindex(lengths)
        points[i, :] = curve.transform(@view points[i, :])
        lengths[i] *= curve.transform.scale
    end

    return segments
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
    @assert n > 0

    if n > 1
        xs, ys = (range(x1, x2, n) for (x1, x2) in zip(line.p1, line.p2))
        ds = hypot(xs[2] - xs[1], ys[2] - ys[1])
    else
        # TEMP: only uses the endpoint if there's one point
        xs, ys = (range(x1, x1, 1) for x1 in line.p2)
        ds = norm(line.p2 - line.p1)
    end

    points = [xs ys]
    lengths = fill(ds, n)

    return Segments(points, lengths)
end

struct SplitCurve{C<:Curve,N} <: Curve
    base::C
    splits::NTuple{N,LineSegment}
end

isclosed(curve::SplitCurve) = isclosed(curve.base)
arclength(curve::SplitCurve) = arclength(curve.base)

function partition(curve::SplitCurve, n::Integer)
    segments = partition(curve.base, n)
    if segments isa SplitSegments
        # TODO: Throw more useful error, or handle this case in a better way
        # If you split a curve twice, it's not clear which new split should correspond to
        # which part
        error("Cannot split curve twice")
    end

    points = Curves.points(segments)
    lengths = Curves.lengths(segments)
    ranges = cut_indices(points, curve.splits)
    return SplitSegments(points, lengths, ranges)
end

"""
    separate(curve::Curve, cut_at) :: SplitCurve
    separate(cut_at)(curve::Curve) :: SplitCurve

Mark cuts in `curve` such that a [`SplitSegments`](@ref) is returned when calling
[`partition`](@ref). `cut_at` is a tuple of [`LineSegment`](@ref)s that cut the curve.
"""
separate(curve::Curve, cut_at::Tuple{Vararg{LineSegment}}) = SplitCurve(curve, cut_at)
separate(cut_at::Tuple{Vararg{LineSegment}}) = curve::Curve -> SplitCurve(curve, cut_at)

"""
    separate(segments::Segments, cut_at::Tuple{Vararg{LineSegment}})

Separate `segments` into a tuple of [`Segments`](@ref). Separated in the same manner as
[`cut_ indices`](@ref).
"""
function separate(segments::Segments, cut_at::Tuple{Vararg{LineSegment}})
    points = Curves.points(segments)
    lengths = Curves.lengths(segments)
    ranges = cut_indices(points, cut_at)
    return map(ranges) do r
        @views Segments(points[r, :], lengths[r])
    end
end

"""
    cut_indices(points::AbstractMatrix, cut_at::Tuple{Vararg{LineSegment}}) -> ranges

Split the indices of points into a tuple of index ranges.

`points` is a matrix interpreted as a sequence of `[x y]` points. The points are split at
each sequential segment in `cut_at` such that `cut_at[i]` is between `ranges[i]` and
`ranges[i + 1]`. The range of points at index `i` is `points[ranges[i], :]`.
"""
function cut_indices(points::AbstractMatrix, cut_at::NTuple{N,LineSegment}) where {N}
    ranges = ntuple(_ -> Ref(1:2), N + 1)

    n = size(points, 1)
    n == 0 && return (1:0,)

    i1 = 1 # Index after last cut
    i2 = 1 # Current point index
    for j in eachindex(cut_at)
        i1 > n && error("Cannot find intersection with cut_at[$j]")

        line = cut_at[j]
        along, across = let v = line.p2 - line.p1, (x, y) = v
            (v / dot(v, v), SVector(-y, x))
        end

        p = @view points[i2, :]
        d1 = signbit(dot(p - line.p1, across))

        while true
            i2 += 1
            if i2 > n
                error("Cannot find intersection with cut_at[$j]")
            end

            p = @view points[i2, :]
            v = p - line.p1
            d2 = signbit(dot(v, across))

            # Dot product sign with `across` changes means we crossed the line
            # Dot product with `along` between 0 and 1 means the point is near the segment
            # If both, we found the intersection
            if d1 != d2 && 0 ≤ dot(v, along) ≤ 1
                break
            end
            d1 = d2
        end

        # Separate right before the point that crossed the line segment
        ranges[j][] = i1:(i2 - 1)
        i1 = i2
    end

    ranges[end][] = i1:n

    return map(i -> i[], ranges)
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

"""
    ParameterizedCurve(f; n_sample=100) :: Curve

A curve defined by the points `(x, y) = f(t)` for `0 ≤ t ≤ 1`.

`n_sample` points are used to determine an equal spacing when partitioning the curve.
"""
struct ParameterizedCurve{F,P1,P2} <: Curve
    f::F # parameter -> point
    param::P1 # normalized arclength -> parameter
    inv_param::P2 # parameter -> normalized arclength
    arclen::Float64
    closed::Bool # whether f(0) ≈ f(1)
    function ParameterizedCurve(f; n_sample=100)
        fv = SVector{2,Float64} ∘ f

        t_sample = range(0, 1, n_sample)

        s_sample = sample_arclengths(fv, t_sample)
        arclen = s_sample[end]

        s_sample ./= arclen # Normalize to [0, 1]

        inv_param = let s = interpolate(s_sample, BSpline(Linear()))
            Interpolations.scale(s, t_sample)
        end
        param = interpolate((s_sample,), t_sample, Gridded(Linear()))

        closed = norm(fv(1) - fv(0)) / arclen < 1e-8

        return let F = typeof(fv), P1 = typeof(param), P2 = typeof(inv_param)
            return new{F,P1,P2}(fv, param, inv_param, arclen, closed)
        end
    end
end

arclength(curve::ParameterizedCurve) = curve.arclen
isclosed(curve::ParameterizedCurve) = curve.closed
(curve::ParameterizedCurve)(t) = curve.f(curve.param(t))

function partition(curve::ParameterizedCurve, n::Integer)
    points = Matrix{Float64}(undef, n, 2)

    ts = if isclosed(curve)
        range(0, 1, n + 1)[1:n]
    else
        range(0, 1, n)
    end

    for (i, t) in enumerate(ts)
        points[i, :] .= curve(t)
    end

    lengths = point_array_lengths(points, isclosed(curve))
    return Segments(points, lengths)
end

function sample_arclengths(f, t)
    s = similar(t)

    p1 = f(t[1])
    s[1] = 0
    for i in eachindex(s)[2:end]
        p2 = f(t[i])
        s[i] = s[i - 1] + norm(p2 - p1)
        p1 = p2
    end

    return s
end

function point_array_lengths(points, closed::Bool)
    lengths = Vector{Float64}(undef, size(points, 1))
    return point_array_lengths!(lengths, points, closed)
end

function point_array_lengths!(lengths, points, closed::Bool)
    n = size(points, 1)
    @assert n > 1 "cannot determine length with single point"

    len(p, i, j) = norm(SVector{2}(p[i, :]) - SVector{2}(p[j, :]))

    lens = Vector{Float64}(undef, n - 1)
    for i in eachindex(lens)
        lens[i] = len(points, i, i + 1)
    end

    if closed
        endlen = len(points, n, 1)
        lengths[1] = (endlen + lens[1]) / 2
        lengths[n] = (lens[n - 1] + endlen) / 2
    else
        lengths[1] = lens[1]
        lengths[n] = lens[n - 1]
    end
    for i in 2:(n - 1)
        lengths[i] = (lens[i - 1] + lens[i]) / 2
    end

    return lengths
end

abstract type NacaParams end

function NacaParams(spec::AbstractString)
    if !all(isdigit, spec)
        throw(ArgumentError("NACA specification must only contain digits"))
    end

    n = length(spec)
    return if n == 4
        NACA4(spec)
    else
        throw(ArgumentError("$n-digit NACA airfoil not implemented"))
    end
end

macro naca_str(spec::AbstractString)
    params = NacaParams(spec)
    return :(NACA($params))
end

struct NACA{N<:NacaParams,C<:ParameterizedCurve} <: ClosedCurve
    params::N
    curve::C
    leading_edge::Float64 # Parameter of leading edge
    tailcut::Float64 # Portion of chord to flatten the trailing edge by
    function NACA(params::NacaParams; tailcut=0.0, open=false, n_sample=2000)
        curve = if open
            ParameterizedCurve(s -> _point(params, s, tailcut); n_sample)
        else
            p1 = _point(params, 0.0, tailcut)
            p2 = _point(params, 1.0, tailcut)
            cut_length = norm(p2 - p1)
            arclen = arclength(ParameterizedCurve(s -> _point(params, s, tailcut)))
            cut_s = cut_length / (cut_length + arclen)
            s1 = cut_s / 2
            s2 = 1 - cut_s / 2
            ParameterizedCurve(; n_sample) do s
                if s < s1
                    t = 0.5 + s / cut_s
                    cut = true
                elseif s > s2
                    t = (s - s2) / cut_s
                    cut = true
                else
                    cut = false
                end

                if cut
                    p2 + (p1 - p2) * t
                else
                    _point(params, (s - s1) / (s2 - s1), tailcut)
                end
            end
        end

        leading_edge = curve.inv_param(0.5)
        return new{typeof(params),typeof(curve)}(params, curve, leading_edge, tailcut)
    end
end

(curve::NACA)(s) = curve.curve(s)
arclength(curve::NACA) = arclength(curve.curve)
partition(curve::NACA, n::Integer) = partition(curve.curve, n)

leading_edge(airfoil::NACA) = airfoil.leading_edge
trailing_edge(airfoil::NACA) = 0.0

struct NACA4 <: NacaParams
    m::Float64 # Max camber
    p::Float64 # Location of max camber along chord
    t::Float64 # Max thickness
end

function NACA4(spec::AbstractString)
    @assert all(isdigit, spec)
    m = parse(Int8, spec[1]) / 100
    p = parse(Int8, spec[2]) / 10
    t = parse(Int8, spec[3:4]) / 100
    return NACA4(m, p, t)
end

# Equations from http://airfoiltools.com/airfoil/naca4digit

function naca_width(params::NACA4, x::Real)
    (; m, p, t) = params

    # NACA airfoil constants
    a0 = 0.2969
    a1 = -0.126
    a2 = -0.3516
    a3 = 0.2843
    a4 = -0.1036 # closed trailing edge

    5 * t * (a0 * sqrt(x) + a1 * x + a2 * x^2 + a3 * x^3 + a4 * x^4)
end

function _point(params::NACA4, s::Real, tailcut::Real=0.0)
    (; m, p, t) = params

    β = 2 * π * s
    x = 0.5 * (1 + cos(β)) * (1 - tailcut)

    if x < p
        yc = m / p^2 * (2 * p * x - x^2)
        dyc = 2 * m / p^2 * (p - x)
    else
        yc = m / (1 - p)^2 * (1 - 2 * p + 2 * p * x - x^2)
        dyc = 2 * m / (1 - p)^2 * (p - x)
    end

    yt = naca_width(params, x)

    θ = atan(dyc)
    z = SVector(x, yc)
    dz = yt * SVector(-sin(θ), cos(θ))

    return s < 0.5 ? z + dz : z - dz
end

end # module Curves
