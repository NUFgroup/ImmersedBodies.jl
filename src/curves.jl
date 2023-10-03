module Curves

using ..ImmersedBodies: Fluid, Panels, CartesianGrid, MultiDomainGrid

using Base: @kwdef
using LinearAlgebra: norm
using StaticArrays
using Interpolations
using LinearAlgebra: I

export Curve, ClosedCurve, OpenCurve, isclosed, arclength, partition
export Transformation, translate, rotate, scale, TransformedCurve
export LineSegment, Circle, ParameterizedCurve
export NacaParams, NACA, NACA4, NACA5, @naca_str, leading_edge, trailing_edge

"""
    Curve

A mathematical curve in 2D space. Calling `curve(t)` with `0 ≤ t ≤ 1` returns a point along
the curve.
"""
abstract type Curve end

"""
    ClosedCurve <: Curve

A curve where `curve(0.0) == curve(1.0)`.
"""
abstract type ClosedCurve <: Curve end

"""
    OpenCurve <: Curve

A curve where `curve(0.0) != curve(1.0)`.
"""
abstract type OpenCurve <: Curve end

"""
    isclosed(curve::Curve)

Return true if the `curve`'s starting and ending points are equal, false otherwise.
"""
isclosed(::ClosedCurve) = true
isclosed(::OpenCurve) = false

"""
    Panels(grid::CartesianGrid, curve::Curve)
    Panels(grid::MultiDomainGrid, curve::Curve)
    Panels(fluid::Fluid, curve::Curve)

Discretize the `curve` into panels that work well with the given `fluid` or `grid`.
"""
Panels(grid::CartesianGrid, curve::Curve) = partition(curve, 2 * grid.h)
Panels(grid::MultiDomainGrid, curve) = Panels(grid.base, curve)
Panels(fluid::Fluid, curve) = Panels(fluid.grid.base, curve)

"""
    arclength(curve::Curve)

The total arclength of a curve.
"""
function arclength end

"""
    partition(curve::Curve, ds::AbstractFloat) :: Panels
    partition(curve::Curve, n::Integer) :: Panels

Partition a curve into segments of approximately equal length. Either specify the target
segment length `ds` or the target segment count `n`. The curve's endpoints are preserved.
"""
function partition end

# Curves either need to define a method for ds or for n
# The other is satisfied by these defaults
function partition(curve::Curve, ds::AbstractFloat)
    n = round(Int, arclength(curve) / ds) + !isclosed(curve)
    partition(curve, n)
end
function partition(curve::Curve, n::Integer)
    ds = arclength(curve) / (n - !isclosed(curve))
    partition(curve, ds)
end

"""
    Transformation

A 2-D transformation that is the combination of translation, rotation, and uniform scaling.

Constructed using [`translate`](@ref), [`rotate`](@ref), and [`scale`](@ref).
"""
struct Transformation
    offset::SVector{2,Float64}  # Translation vector
    rotation::SMatrix{2,2,Float64,4}  # Rotation matrix
    scale::Float64  # Uniform scale
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

const I2 = SMatrix{2,2,Float64}(I)
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
    panels = partition(curve.base, n)

    for i in eachindex(panels.ds)
        pt = SVector{2}(@view panels.xb[i, :])
        panels.xb[i, :] = curve.transform(pt)
        panels.ds[i] *= curve.transform.scale
    end

    panels
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

arclength(line::LineSegment) = norm(line.p2 - line.p1)

function Base.show(io::IO, ::MIME"text/plain", line::LineSegment)
    print(io, "LineSegment: from ", line.p1, " to ", line.p2)
    nothing
end

function partition(line::LineSegment, n::Integer)
    x1, y1 = line.p1
    x2, y2 = line.p2
    xb = [range(x1, x2, n) range(y1, y2, n)]
    ds = fill(arclength(line) / (n - 1), n)
    Panels(; xb, ds)
end

"""
    Circle(; r=1.0, center=(0.0, 0.0)) :: Curve

A circle with radius `r` centered at `center`.
"""
@kwdef struct Circle <: ClosedCurve
    r::Float64 = 1.0
    center::SVector{2,Float64} = (0.0, 0.0)
end

function (circle::Circle)(t)
    s = 2 * pi * t
    circle.center + circle.r * SVector(cos(s), sin(s))
end

arclength(circle::Circle) = 2 * pi * circle.r

function Base.show(io::IO, ::MIME"text/plain", circle::Circle)
    print(io, "Circle: radius=", circle.r, ", center=", circle.center)
    nothing
end

function partition(circle::Circle, n::Integer)
    x0, y0 = circle.center
    r = circle.r

    t = 2pi / n * (0:(n - 1))

    xs = @. x0 + r * cos(t)
    ys = @. y0 + r * sin(t)
    ds = hypot(xs[2] - xs[1], ys[2] - ys[1])

    Panels(; xb=[xs ys], ds=fill(ds, n))
end

"""
    ParameterizedCurve(f; n_sample=100) :: Curve

A curve defined by the points `(x, y) = f(t)` for `0 ≤ t ≤ 1`. `n_sample` points are used to
determine an equal spacing when partitioning the curve.
"""
struct ParameterizedCurve{F,P1,P2} <: Curve
    f::F  # parameter -> point
    param::P1  # normalized arclength -> parameter
    inv_param::P2  # parameter -> normalized arclength
    arclen::Float64
    closed::Bool  # whether f(0) ≈ f(1)
    function ParameterizedCurve(f; n_sample=100)
        fv = SVector{2,Float64} ∘ f

        t_sample = range(0, 1, n_sample)

        s_sample = sample_arclengths(fv, t_sample)
        arclen = s_sample[end]

        s_sample ./= arclen  # Normalize to [0, 1]

        inv_param = let s = interpolate(s_sample, BSpline(Linear()))
            Interpolations.scale(s, t_sample)
        end
        param = interpolate((s_sample,), t_sample, Gridded(Linear()))

        closed = norm(fv(1) - fv(0)) / arclen < 1e-10

        let F = typeof(fv), P1 = typeof(param), P2 = typeof(inv_param)
            new{F,P1,P2}(fv, param, inv_param, arclen, closed)
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
    Panels(; xb=points, ds=lengths)
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

    s
end

function point_array_lengths(points, closed::Bool)
    lengths = Vector{Float64}(undef, size(points, 1))
    point_array_lengths!(lengths, points, closed)
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

    lengths
end

abstract type NacaParams end

NacaParams(spec::AbstractString) = NacaParams(String(spec))
function NacaParams(spec::String)
    if !all(isdigit, spec)
        throw(ArgumentError("NACA specification must only contain digits"))
    end
    NacaParams(Val(ncodeunits(spec)), spec)
end
function NacaParams(::Val{N}, ::String) where {N}
    throw(ArgumentError("$N-digit NACA airfoil not defined"))
end

macro naca_str(spec::AbstractString)
    params = NacaParams(spec)
    :(NACA($params))
end

struct NACA{N<:NacaParams,C<:ParameterizedCurve} <: ClosedCurve
    params::N
    curve::C
    leading_edge::Float64  # Parameter of leading edge
    function NACA(params::NacaParams)
        curve = ParameterizedCurve(s -> _point(params, s))
        leading_edge = curve.inv_param(0.5)
        new{typeof(params),typeof(curve)}(params, curve, leading_edge)
    end
end

(curve::NACA)(s) = curve.curve(s)
arclength(curve::NACA) = arclength(curve.curve)
partition(curve::NACA, n::Integer) = partition(curve.curve, n)

leading_edge(airfoil::NACA) = airfoil.leading_edge
trailing_edge(airfoil::NACA) = 0.0

struct NACA4 <: NacaParams
    m::Float64  # Max camber
    p::Float64  # Location of max camber along chord
    t::Float64  # Max thickness
end

function NacaParams(::Val{4}, spec::String)
    @assert all(isdigit, spec)
    m = parse(Int8, spec[1]) / 100
    p = parse(Int8, spec[2]) / 10
    t = parse(Int8, spec[3:4]) / 100
    NACA4(m, p, t)
end

function _point(params::NACA4, s::Real)
    # Equations from http://airfoiltools.com/airfoil/naca4digit
    (; m, p, t) = params

    # NACA airfoil constants
    a0 = 0.2969
    a1 = -0.126
    a2 = -0.3516
    a3 = 0.2843
    a4 = -0.1036  # closed trailing edge

    β = 2 * π * s
    x = 0.5 * (1 + cos(β))

    if x < p
        yc = m / p^2 * (2 * p * x - x^2)
        dyc = 2 * m / p^2 * (p - x)
    else
        yc = m / (1 - p)^2 * (1 - 2 * p + 2 * p * x - x^2)
        dyc = 2 * m / (1 - p)^2 * (p - x)
    end

    yt = 5 * t * (a0 * sqrt(x) + a1 * x + a2 * x^2 + a3 * x^3 + a4 * x^4)

    θ = atan(dyc)
    z = SVector(x, yc)
    dz = yt * SVector(-sin(θ), cos(θ))

    s < 0.5 ? z + dz : z - dz
end

struct NACA5 <: NacaParams
    l::Float64  #lift coeff . Multiply number by 3/20 to get Cl
    p::Float64  # Location of max camber along chord . Divide number by 20 to get Location
    q::Float64  # reflex or normal
    xx::Float64   #max thickness
end

function NacaParams(::Val{5}, spec::String)
    @assert all(isdigit, spec)
    l = parse(Int8, spec[1]) * (3 / 20)
    p = parse(Int8, spec[2])
    q = parse(Int8, spec[3])
    xx = parse(Int8, spec[4:5]) * 0.01
    NACA5(l, p, q, xx)
end

function _point(params::NACA5, s::Real)
    # Equations from http://airfoiltools.com/airfoil/naca5digit
    (; l, p, q, xx) = params

    β = 2 * π * s
    x = 0.5 * (1 + cos(β))

    #Getting r, k1,k2/k1 values
    if (10 * p + q) == 10
        r = 0.0580 * (l / 0.3)
        k1 = 361.4 * (l / 0.3)
    elseif (10 * p + q) == 20
        r = 0.126 * (l / 0.3)
        k1 = 51.64 * (l / 0.3)
    elseif (10 * p + q) == 30
        r = 0.2025 * (l / 0.3)
        k1 = 15.957 * (l / 0.3)
    elseif (10 * p + q) == 40
        r = 0.2900 * (l / 0.3)
        k1 = 6.643 * (l / 0.3)
    elseif (10 * p + q) == 50
        r = 0.3910 * (l / 0.3)
        k1 = 3.23 * (l / 0.3)
    elseif (10 * p + q) == 21
        r = 0.1300 * (l / 0.3)
        k1 = 51.99 * (l / 0.3)
        k2k1 = 0.000764 * (l / 0.3)
    elseif (10 * p + q) == 31
        r = 0.217 * (l / 0.3)
        k1 = 15.793 * (l / 0.3)
        k2k1 = 0.00677 * (l / 0.3)
    elseif (10 * p + q) == 41
        r = 0.3180 * (l / 0.3)
        k1 = 6.520 * (l / 0.3)
        k2k1 = 0.0303 * (l / 0.3)
    elseif (10 * p + q) == 51
        r = 0.4410 * (l / 0.3)
        k1 = 3.191 * (l / 0.3)
        k2k1 = 0.1355 * (l / 0.3)
    else
        println("Wrong airfoil")
    end

    # NACA airfoil constants
    a0 = 0.2969
    a1 = -0.126
    a2 = -0.3516
    a3 = 0.2843
    a4 = -0.1036  # closed trailing edge

    #Defining camber line
    if x < r && q == 0
        yc = (k1 / 6) * (x^3 - 3 * r * (x^2) + (r^2) * (3 - r) * x)
        dyc = (k1 / 6) * (3 * x^2 - 6 * r * x + (r^2) * (3 - r))
    elseif x < r && q == 1
        yc = (k1 / 6) * ((x - r)^3 - k2k1 * ((1 - r)^3) * x - (r^3) * x + r^3)
        dyc = (k1 / 6) * (3 * (x - r)^2 - k2k1 * (1 - r)^3 - r^3)
    elseif x > r && q == 0
        yc = (k1 / 6) * (r^3) * (1 - x)
        dyc = (-k1 / 6) * r^3
    elseif x > r && q == 1
        yc = (k1 / 6) * (k2k1 * (x - r)^3 - k2k1 * ((1 - r)^3) * x - (r^3) * x + r^3)
        dyc = (k1 / 6) * (3 * k2k1 * (x - r)^2 - k2k1 * (1 - r)^3 - r^3)
    end

    yt = 5 * xx * (a0 * sqrt(x) + a1 * x + a2 * x^2 + a3 * x^3 + a4 * x^4)

    θ = atan(dyc)
    z = SVector(x, yc)
    dz = yt * SVector(-sin(θ), cos(θ))

    s < 0.5 ? z + dz : z - dz
end

end # module
