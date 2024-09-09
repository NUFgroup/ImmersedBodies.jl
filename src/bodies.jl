module Bodies

using ..ImmersedBodies
using ..ImmersedBodies.Dynamics
using ..ImmersedBodies.Curves
using ..ImmersedBodies.PointSpacing
import ..ImmersedBodies: _show

using LinearAlgebra

using StaticArrays

export AbstractBody, AbstractPrescribedBody, BodyGroup, Panels, PanelView, npanels, bodypanels, deformation
export AbstractBodyPoint, BodyPointIndex, BodyPointParam
export RigidBody, PrescribedBody, DeformingBody, EulerBernoulliBeam, is_static
export SpringedMembrane, diatomic_phononic_material
export DeformingBodyBC, bc_point, ClampBC, PinBC
export reference_pos, n_variables, deforming
export StructureModel, LinearModel
export DeformationState, DeformationStateView

# typeof(@view matrix[i:j, :])
const MatrixRowView{T} = SubArray{
    T,2,Matrix{T},Tuple{UnitRange{Int},Base.Slice{Base.OneTo{Int}}},false
}

# typeof(@view vector[i:j])
const VectorView{T} = SubArray{T,1,Vector{T},Tuple{UnitRange{Int}},true}

"""
    AbstractBody

A structural body.
"""
abstract type AbstractBody end

abstract type AbstractPrescribedBody <: AbstractBody end

"""
    AbstractBodyPoint

A point on a body.
"""
abstract type AbstractBodyPoint end

"""
    BodyPointIndex(i) :: AbstractBodyPoint

The point on a body at index `i`.
"""
struct BodyPointIndex <: AbstractBodyPoint
    i::Int
end

"""
    BodyPointParam(t) :: AbstractBodyPoint

The point nearest to a portion of `t` along the arclength of a body.
"""
struct BodyPointParam <: AbstractBodyPoint
    t::Float64
    function BodyPointParam(t)
        if !(0 <= t <= 1)
            (throw ∘ DomainError)(t, "t must be between 0 and 1")
        end
        return new(t)
    end
end

Base.show(io::IO, ::MIME"text/plain", body::AbstractBody) = _show(io, body)

"""
    initial_pos!(xb::AbstractMatrix{Float64}, body::AbstractBody)

The initial `[x y]` panels on a body.
"""
function initial_pos! end

"""
    initial_lengths!(ds::AbstractVector{Float64}, body::AbstractBody)

The initial body segment lengths.
"""
function initial_lengths! end

"""
    npanels(body::AbstractBody)
    npanels(bodies::BodyGroup)
    npanels(panels::Panels)

The number of structural panels in a body or bodies.
"""
function npanels end

"""
    RigidBody(
        pos::AbstractMatrix{Float64},
        len::AbstractVector{Float64},
        frame::AbstractFrame=DiscretizationFrame()
    ) :: AbstractBody

A rigid body with optional prescribed motion.

# Arguments
- `pos`: Rows of `[x y]` points that define the body.
- `len`: Segment length at each point in `pos`.
- `motion`: Prescribed motion of the body.
"""
struct RigidBody{F<:AbstractFrame} <: AbstractPrescribedBody
    pos::Matrix{Float64} # panel positions
    len::Vector{Float64} # panel lengths
    frame::F # prescribed motion
end

RigidBody(pos, len) = RigidBody(pos, len, DiscretizationFrame())
initial_pos!(xb, body::RigidBody) = xb .= body.pos
initial_lengths!(ds, body::RigidBody) = ds .= body.len

"""
    RigidBody(segments::Segments, frame=DiscretizationFrame()) :: AbstractBody

# Arguments
- `segments::Segments`: Points and lengths on the body.
- `frame::AbstractFrame`: Prescribed motion of the body.
"""
function RigidBody(segments::Segments, frame=DiscretizationFrame())
    return RigidBody(Curves.points(segments), Curves.lengths(segments), frame)
end

npanels(body::RigidBody) = length(body.len)

is_static(body::RigidBody{DiscretizationFrame}, ::AbstractFrame) = true
is_static(body::RigidBody{F}, ::F) where {F<:BaseFrame} = true
is_static(body::RigidBody, ::AbstractFrame) = false

function _show(io::IO, body::RigidBody, prefix)
    print(io, prefix, "RigidBody:")
    if get(io, :compact, false)
        print(io, " in frame ", body.frame)
    else
        ioc = IOContext(io, :limit => true, :compact => true)

        print(ioc, '\n', prefix, "  points = ")
        summary(ioc, body.pos)

        print(ioc, '\n', prefix, "   frame = ")
        summary(ioc, body.frame)
    end

    return nothing
end

struct PrescribedBody{F1,F2,F3} <: AbstractPrescribedBody
    n_panels::Int
    x::F1
    v::F2
    ds::F3
    auto_ds::Bool
end

PrescribedBody(n_panels, x, v, ds) = PrescribedBody(n_panels, x, v, ds, false)

initial_pos!(xb, body::PrescribedBody) = body.x(xb, 0)
initial_lengths!(ds, body::PrescribedBody) = body.ds(ds, 0)

npanels(body::PrescribedBody) = body.n_panels

is_static(::PrescribedBody, _) = false

abstract type DeformingBody <: AbstractBody end

abstract type StructureModel end

Base.@kwdef struct LinearModel <: StructureModel
    m::Vector{Float64}
    kb::Vector{Float64}
end

abstract type DeformingBodyBC{P<:AbstractBodyPoint} end

# TODO: Only shorten the type name when displaying a value of this type
#       println(bc) should shorten the type name
#       println(typeof(bc)) should keep the entire type name
function Base.show(io::IO, T::Type{<:DeformingBodyBC})
    return print(io, nameof(T)) # Omit type parameters by default
end

"""
    ClampBC(point::AbstractBodyPoint) :: DeformingBodyBC

Fixes the position and rotation of a point on a deforming body.
"""
struct ClampBC{P} <: DeformingBodyBC{P}
    point::P
end

bc_point(bc::ClampBC) = bc.point
set_point(::ClampBC, p::AbstractBodyPoint) = ClampBC(p)

"""
    PinBC(point::AbstractBodyPoint) :: DeformingBodyBC

Fixes the position of a point on a deforming body.
"""
struct PinBC{P} <: DeformingBodyBC{P}
    point::P
end

bc_point(bc::PinBC) = bc.point
set_point(::PinBC, p::AbstractBodyPoint) = PinBC(p)

struct EulerBernoulliBeam{M<:StructureModel,B<:DeformingBodyBC{BodyPointIndex}} <:
       DeformingBody
    model::M
    xref::Matrix{Float64} # Reference locations about which displacements are determined
    ds0::Vector{Float64} # Line segment lengths on body in undeformed configuration
    bcs::Vector{B} # Boundary conditions
end

reference_pos(body::EulerBernoulliBeam) = body.xref

# 2 displacements per poin
n_variables(body::EulerBernoulliBeam{LinearModel}) = 2 * npanels(body)

initial_pos!(xb, body::EulerBernoulliBeam) = xb .= body.xref
initial_lengths!(ds, body::EulerBernoulliBeam) = ds .= body.ds0

npanels(body::EulerBernoulliBeam) = length(body.ds0)

function EulerBernoulliBeam(
    ::Type{LinearModel},
    segments::Segments,
    bcs::AbstractVector{<:ClampBC};
    m::Float64,
    kb::Float64,
)
    xref = segments.points
    ds0 = segments.lengths

    nel = size(xref, 1) - 1
    model = LinearModel(; m=fill(m, nel), kb=fill(kb, nel))

    ts = cumsum(ds0)
    @. ts = (ts - ts[1]) / (ts[end] - ts[1])

    # Make sure bcs are based on indices and not parameters along the curve
    bc_indices = [to_index_bc(ts, bc) for bc in bcs]

    return EulerBernoulliBeam(model, xref, ds0, bc_indices)
end

"""
    search_closest(xs, x)

Return the index of the closest element in `xs` to `x` assuming `xs` is sorted.
"""
function search_closest(xs, x)
    i2 = searchsortedfirst(xs, x)

    i2 > lastindex(xs) && return lastindex(xs)
    i2 == firstindex(xs) && return firstindex(xs)

    i1 = i2 - 1
    return xs[i2] - x > x - xs[i1] ? i1 : i2
end

# Ensure the boundary condition is based on an index along the body
function to_index_bc(::AbstractVector{Float64}, bc::DeformingBodyBC{BodyPointIndex})
    return bc
end
function to_index_bc(ts::AbstractVector{Float64}, bc::DeformingBodyBC{BodyPointParam})
    point = bc_point(bc)
    i = search_closest(ts, point.t)
    return set_point(bc, BodyPointIndex(i))
end

function _show(io::IO, body::EulerBernoulliBeam, prefix)
    print(io, prefix)
    summary(io, body)
    print(io, ':')
    if get(io, :compact, false)
        print(io, " with bcs ", body.bcs)
    else
        ioc = IOContext(io, :limit => true, :compact => true)

        print(ioc, '\n', prefix, "  reference points = ")
        summary(ioc, body.xref)

        print(ioc, '\n', prefix, "    boundary conds = ", body.bcs)
    end

    return nothing
end

struct SpringedMembrane{F} <: DeformingBody
    deformed::F
    s_undef_end::Float64
    xref::Matrix{Float64} # Reference position of membrane
    ds0::Vector{Float64} # Lengths at each membrane reference position
    normals::Matrix{Float64} # Normals at each point in membrane
    deform_weights::Matrix{Float64} # Weights for deforming membrane
    spring_normal::SVector{2,Float64} # Direction that the spring can respond in
    m::Vector{Float64}
    k::Vector{Float64}
    kg::Vector{Float64} # Grounded spring constants
end

reference_pos(body::SpringedMembrane) = body.xref

n_variables(body::SpringedMembrane) = length(body.m)

initial_pos!(xb, body::SpringedMembrane) = xb .= body.xref
initial_lengths!(ds, body::SpringedMembrane) = ds .= body.ds0

npanels(body::SpringedMembrane) = length(body.ds0)

function SpringedMembrane(
    curve;
    ds,
    m::AbstractVector{Float64},
    k::AbstractVector{Float64},
    kg::AbstractVector{Float64},
    align_normal,
)
    segments = partition(ParameterizedCurve(curve), ds)

    xref = segments.points
    ds0 = segments.lengths

    nb = size(xref, 1)

    @assert axes(m) == axes(k)

    s = cumsum(ds0)

    # -1 to 1 across compliant section
    r = @. 2 * (s - s[1]) / (s[end] - s[1]) - 1

    function secant(i1, i2)
        dx = xref[i2, 1] - xref[i1, 1]
        dy = xref[i2, 2] - xref[i1, 2]
        return SVector(dx, dy) / hypot(dx, dy)
    end

    function rotate(v)
        vx, vy = v
        return SVector(-vy, vx)
    end

    i1 = nb ÷ 2
    i2 = i1 + 1
    flipnormal = dot(rotate(secant(i1, i2)), align_normal) < 0

    function normal(i1, i2)
        v = rotate(secant(i1, i2))
        return flipnormal ? -v : v
    end

    normals = zeros(nb, 2)
    normals[1, :] .= normal(1, 2)
    normals[end, :] .= normal(nb - 1, nb)
    for i in 2:(nb - 1)
        normals[i, :] .= normal(i - 1, i + 1)
    end

    # Take the normal at the middle point to be the spring direction
    spring_normal = SVector{2}(@view normals[end ÷ 2, :])

    deform_weights = zeros(nb, 2)

    s_undef_end = PointSpacing.arclength(curve; sample=500)
    function deformed(k)
        function(t)
            SVector(curve(t)) + k * spring_normal * membrane_distribution(2t-1)
        end
    end

    return SpringedMembrane(deformed, s_undef_end, xref, ds0, normals, deform_weights, spring_normal, m, k, kg)
end

membrane_distribution(x) = exp(-(3x)^2 / 2)

function diatomic_phononic_material(
    segments::Segments;
    n_cell::Int,
    m::NTuple{2,Float64},
    k::NTuple{2,Float64},
    kg::NTuple{2,Float64}=(0.0, 0.0),
    first_spring::Bool,
    kw...,
)
    n_spring = 2 * n_cell

    mvec = Vector{Float64}(undef, n_spring)
    kvec = Vector{Float64}(undef, n_spring)
    kgvec = Vector{Float64}(undef, n_spring)

    m1, m2 = m
    k1, k2 = k
    kg1, kg2 = kg
    for i1 in 1:2:n_spring
        i2 = i1 + 1
        mvec[i1] = m1
        mvec[i2] = m2
        kvec[i1] = k1
        kvec[i2] = k2
        kgvec[i1] = kg1
        kgvec[i2] = kg2
    end

    if !first_spring
        kvec[1] = 0.0
    end

    return SpringedMembrane(segments; m=mvec, k=kvec, kg=kgvec, kw...)
end

struct DeformingBodyIndex
    i_body::Int
    i_panel::Int
    i_deform_panel::Int
end

"""
    BodyGroup(bodies::Vector{<:AbstractBody})

A collection of bodies.
"""
struct BodyGroup{B<:AbstractBody} <: AbstractVector{B}
    bodies::Vector{B} # all bodies
    npanel::Int
    npanel_deform::Int
    deforming::Vector{DeformingBodyIndex}
    index_to_deform::Dict{Int,Int}
    function BodyGroup(bodies::Vector{B}) where {B<:AbstractBody}
        deforming = DeformingBodyIndex[]
        index_to_deform = Dict{Int,Int}()

        n_panel = 0
        n_panel_deform = 0
        for (i_body, body) in enumerate(bodies)
            n = npanels(body)

            if body isa DeformingBody
                push!(deforming, DeformingBodyIndex(i_body, n_panel, n_panel_deform))
                index_to_deform[i_body] = lastindex(deforming)
                n_panel_deform += n
            end

            n_panel += n
        end

        return new{B}(bodies, n_panel, n_panel_deform, deforming, index_to_deform)
    end
end

Base.size(bodies::BodyGroup) = size(bodies.bodies)
Base.getindex(bodies::BodyGroup, i) = bodies.bodies[i]
Base.IndexStyle(::BodyGroup) = IndexLinear()

npanels(bodies::BodyGroup) = bodies.npanel

deforming(bodies::BodyGroup) = (bodies[i.i_body] for i in bodies.deforming)

Base.show(io::IO, ::MIME"text/plain", bodies::BodyGroup) = _show(io, bodies)

function _show(io::IO, bodies::BodyGroup, prefix)
    print(io, prefix)
    summary(io, bodies)
    print(io, ":\n")

    indent = prefix * "  "
    for body in bodies
        _show(io, body, indent)
        println(io)
    end
    return nothing
end

"""
    PanelView

View into sequences of panels in [`Panels`](@ref).
"""
struct PanelView
    i_panel::Int
    pos::MatrixRowView{Float64}
    vel::MatrixRowView{Float64}
    len::VectorView{Float64}
    traction::MatrixRowView{Float64}
end

npanels(panels::PanelView) = length(panels.len)

"""
    Panels

A group of structural panels for multiple bodies.

Coordinates are given in the frame that the fluid is discretized in.
"""
struct Panels
    pos::Matrix{Float64} # panel positions
    vel::Matrix{Float64} # panel velocities
    len::Vector{Float64} # panel lengths
    traction::Matrix{Float64} # traction on each panel
    indices::Vector{UnitRange{Int}} # index range for each body
    perbody::Vector{PanelView} # panels grouped by body
end

function Panels(bodies::BodyGroup)
    n = npanels(bodies)

    pos = zeros(n, 2)
    vel = zeros(n, 2)
    len = zeros(n)
    traction = zeros(n, 2)

    indices = Vector{UnitRange}(undef, length(bodies))
    perbody = Vector{PanelView}(undef, length(bodies))
    i_panel = 0
    for (i, body) in enumerate(bodies)
        n_panel = npanels(body)
        r = i_panel .+ (1:n_panel)
        indices[i] = r
        perbody[i] = @views PanelView(i_panel, pos[r, :], vel[r, :], len[r], traction[r, :])
        i_panel += n_panel
    end

    for (body, panels) in zip(bodies, perbody)
        initial_pos!(panels.pos, body)
        initial_lengths!(panels.len, body)
    end

    return Panels(pos, vel, len, traction, indices, perbody)
end

npanels(p::Panels) = length(p.len)

"""
    bodypanels(state::AbstractState) :: Panels

Return the [`Panels`](@ref) of a state.
"""
function bodypanels end

"""
    deformation(state::AbstractState) :: DeformationState

Return the [`DeformationState`](@ref) of a state.
"""
function deformation end

function prescribe_motion!(::PanelView, ::AbstractFrame, ::RigidBody, t::Float64)
    # TODO: Implement prescribed for other reference frames
    throw(ArgumentError("unsupported reference frame"))
end

function prescribe_motion!(
    ::PanelView, ::F1, ::RigidBody{F2}, ::Float64
) where {F1<:BaseFrame,F2<:Union{F1,DiscretizationFrame}}
    return nothing
end

function prescribe_motion!(
    panels::PanelView, ::F1, body::RigidBody{OffsetFrame{F2}}, t::Float64
) where {F1<:BaseFrame,F2<:Union{F1,DiscretizationFrame}}
    f = body.frame(t) # evaluate frame at time t
    r0 = f.r
    v0 = f.v
    Ω = f.Ω
    c = f.cθ # cos(θ)
    s = f.sθ # sin(θ)

    Rx = @SMatrix [c -s; s c]
    Rv = Ω * @SMatrix [-s -c; c -s]

    for (r, v, rb) in zip(eachrow.((panels.pos, panels.vel, body.pos))...)
        r .= r0 + Rx * rb
        v .= v0 + Rv * rb
    end

    return nothing
end

function prescribe_motion!(
    panels::PanelView, ::GlobalFrame,
    body::PrescribedBody, t::Float64
)
    body.x(panels.pos, t)
    body.v(panels.vel, t)
    body.ds(panels.len, t)
    if body.auto_ds
        ImmersedBodies.Solvers.recalculate_lengths!(panels)
    end

    return nothing
end

"""
    DeformationStateView

A view into [`DeformationState`](@ref).
"""
struct DeformationStateView
    χ::VectorView{Float64} # Structural displacements
    ζ::VectorView{Float64} # Structural velocities
    ζdot::VectorView{Float64} # Structural accels
end

"""
    DeformationState

The state of deformation of all deforming bodies.
"""
struct DeformationState
    χ::Vector{Float64} # Structural displacements
    ζ::Vector{Float64} # Structural velocities
    ζdot::Vector{Float64} # Structural accels
    perbody::Vector{DeformationStateView}
end

function DeformationState(vars_per_body)
    n_bodies = length(vars_per_body)
    n_vars = sum(vars_per_body; init=0)

    χ = zeros(n_vars)
    ζ = zeros(n_vars)
    ζdot = zeros(n_vars)

    perbody = Vector{DeformationStateView}(undef, n_bodies)
    i = 0
    for (i_body, n) in enumerate(vars_per_body)
        r = i .+ (1:n)
        perbody[i_body] = @views DeformationStateView(χ[r], ζ[r], ζdot[r])
        i += n
    end

    return DeformationState(χ, ζ, ζdot, perbody)
end

function DeformationState(bodies::BodyGroup)
    return DeformationState(n_variables(body) for body in deforming(bodies))
end

function Base.similar(deform::DeformationState)
    vars_per_body = (length(def.χ) for def in deform.perbody)
    return DeformationState(vars_per_body)
end

function Base.copy!(dst::DeformationState, src::DeformationState)
    # Crude check that both states cover the same bodies
    @assert length(dst.perbody) == length(src.perbody)

    dst.χ .= src.χ
    dst.ζ .= src.ζ
    dst.ζdot .= src.ζdot

    return dst
end

end # module Bodies
