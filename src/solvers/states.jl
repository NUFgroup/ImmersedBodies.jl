# using Debugger

"""
    PsiOmegaState

A state that uses the streamfunction-vorticity formulation.
"""
abstract type PsiOmegaGridState <: AbstractState end

"""
    PsiOmegaGridQuantities

Quantities common to [`PsiOmegaState`](@ref) states.
"""
Base.@kwdef struct PsiOmegaGridQuantities
    q::Matrix{Float64} # Flux
    q0::Matrix{Float64} # Base flux
    Γ::Matrix{Float64} # Circulation
    ψ::Matrix{Float64} # Streamfunction
    F̃b::Vector{Float64} # Body forces * dt
    panels::Panels # Structural panels
    u::MVector{2,Float64} # [ux, uy] freestream velocity
end

function PsiOmegaGridQuantities(prob::Problem{<:PsiOmegaFluidGrid})
    (; nq, nΓ) = prob.fluid.gridindex

    nlevel = (nlevels ∘ discretized)(prob.fluid)

    q = zeros(nq, nlevel)
    q0 = zeros(nq, nlevel)
    Γ = zeros(nΓ, nlevel)
    ψ = zeros(nΓ, nlevel)
    panels = Panels(prob.bodies)
    F̃b = zeros(2 * npanels(panels))
    u = zeros(MVector{2,Float64})

    return PsiOmegaGridQuantities(; q, q0, Γ, ψ, F̃b, panels, u)
end

bodypanels(state::PsiOmegaGridState) = quantities(state).panels

function update_traction!(
    qty::PsiOmegaGridQuantities, prob::Problem{<:PsiOmegaFluidGrid{CNAB}}
)
    # Update the panels field to reflect F̃b.

    h = (gridstep ∘ baselevel ∘ discretized)(prob.fluid)
    dt = timestep(prob)
    ds = qty.panels.len
    F̃b = reshape(qty.F̃b, :, 2)

    @. qty.panels.traction = F̃b * h / (dt * ds)

    return nothing
end

"""
    StatePsiOmegaGridCNAB <: PsiOmegaState

The state of a streamfcn-vorticity formulation simulation with the
Crank-Nicolson/Adams-Bashforth scheme.
"""
mutable struct StatePsiOmegaGridCNAB <: PsiOmegaGridState
    qty::PsiOmegaGridQuantities
    nonlin::Vector{Matrix{Float64}} # Memory of nonlinear terms
    cfl::Float64
    t0::Float64 # Starting time
    i::Int # Index of the current timestep
    t::Float64 # Time of the current timestep
end

timevalue(state::StatePsiOmegaGridCNAB) = state.t
timeindex(state::StatePsiOmegaGridCNAB) = state.i

function StatePsiOmegaGridCNAB(prob::Problem{<:PsiOmegaFluidGrid{CNAB}}, t::Float64)
    nΓ = prob.fluid.gridindex.nΓ
    nlevel = (nlevels ∘ discretized)(prob.fluid)
    nstep = length(timestep_scheme(prob).β)

    qty = PsiOmegaGridQuantities(prob)
    nonlin = [zeros(nΓ, nlevel) for _ in 1:nstep]
    return StatePsiOmegaGridCNAB(qty, nonlin, 0.0, t, 0, t)
end

statetype(::Problem{<:PsiOmegaFluidGrid{CNAB}}) = StatePsiOmegaGridCNAB
quantities(s::StatePsiOmegaGridCNAB) = s.qty

all_levels(prob::Problem{<:PsiOmegaFluidGrid{CNAB}}) = 1:nlevels(discretized(prob.fluid))
grid_quantity(f, coords::AbstractArray{<:Tuple,0}) = GridQuantity(f, coords[])
grid_quantity(f, coords::AbstractVector{<:Tuple}) = MultiLevelGridQuantity(f, coords)

"""
    flow_velocity(direction::Direction, prob::Problem{<:PsiOmegaFluidGrid{CNAB}}; [level])

Return the flow velocity along `direction` on the the `level` grid sublevels.
"""
function Quantities.flow_velocity(
    ::XAxis{DiscretizationFrame},
    prob::Problem{<:PsiOmegaFluidGrid{CNAB}};
    level=all_levels(prob),
)
    grids = discretized(prob.fluid)
    u_coords = [xflux_ranges(sublevel(grids, lev)) for lev in level]
    h = [gridstep(grids, lev) for _ in 1:1, _ in 1:1, lev in level]

    return grid_quantity(u_coords) do state::StatePsiOmegaGridCNAB
        qty = quantities(state)
        u0, _ = split_flux(qty.q0, prob.fluid.gridindex, level)
        u, _ = split_flux(qty.q, prob.fluid.gridindex, level)
        return @. (u0 + u) / h
    end
end

function Quantities.flow_velocity(
    ::YAxis{DiscretizationFrame},
    prob::Problem{<:PsiOmegaFluidGrid{CNAB}};
    level=all_levels(prob),
)
    grids = discretized(prob.fluid)
    v_coords = [yflux_ranges(sublevel(grids, lev)) for lev in level]
    h = [gridstep(grids, lev) for _ in 1:1, _ in 1:1, lev in level]

    return grid_quantity(v_coords) do state::StatePsiOmegaGridCNAB
        qty = quantities(state)
        _, v0 = split_flux(qty.q0, prob.fluid.gridindex, level)
        _, v = split_flux(qty.q, prob.fluid.gridindex, level)
        return @. (v0 + v) / h
    end
end

"""
    streamfunction(prob::Problem{<:PsiOmegaFluidGrid{CNAB}}; [level])

Return the streamfunction on the `level` grid sublevels.
"""
function Quantities.streamfunction(
    prob::Problem{PsiOmegaFluidGrid{CNAB,GlobalFrame}}; level=all_levels(prob)
)
    grids = discretized(prob.fluid)
    coords = [circ_ranges(sublevel(grids, lev)) for lev in level]

    return grid_quantity(coords) do state::StatePsiOmegaGridCNAB
        qty = quantities(state)
        # @bp
        v, u = qty.u
        ψ = unflatten_circ(qty.ψ, prob.fluid.gridindex, level)
        for i in axes(ψ, 3)
            xs, ys = coords[i]
            @views ψ[:, :, i] .+= (u * y - v * x for (x, y) in Iterators.product(xs, ys))
        end

        return ψ
    end
end

function Quantities.streamfunction(
    prob::Problem{PsiOmegaFluidGrid{CNAB,OffsetFrame{GlobalFrame}}}; level=all_levels(prob)
)
    # TODO: Implement for moving domain
    # return error("streamfunction not implemented for moving domains")

    # Implementation for non-rotating motion:
    fluid = prob.fluid
    frame = fluid.frame

    grids = discretized(fluid)
    coords = [circ_ranges(sublevel(grids, lev)) for lev in level]
    # @bp

    return grid_quantity(coords) do state::StatePsiOmegaGridCNAB
        qty = quantities(state)
        v0, u0 = qty.u
        instant = frame.f(state.t)
        # @bp
        vframe, uframe = instant.v
        u0 -= uframe
        v0 -= vframe
        cθ = instant.cθ
        sθ = instant.sθ
        Ω = instant.Ω
        # Ox, Oy = instant.r

        u = u0 * cθ - v0 * sθ
        v = u0 * sθ + v0 * cθ

        ψ = unflatten_circ(qty.ψ, fluid.gridindex, level)

        for i in axes(ψ, 3)
            xs, ys = coords[i]
            @views ψ[:, :, i] .+= ((u - Ω * x) * y - (v + Ω * y) * x for (x, y) in Iterators.product(xs, ys))
        end
        # @bp
        return ψ
    end
end

"""
    vorticity(prob::Problem{<:PsiOmegaFluidGrid{CNAB}}; [level])

Return the vorticity on the `level` grid sublevels.
"""
function Quantities.vorticity(
    prob::Problem{<:PsiOmegaFluidGrid{CNAB}}; level=all_levels(prob)
)
    grids = discretized(prob.fluid)
    coords = [circ_ranges(sublevel(grids, lev)) for lev in level]
    h = [gridstep(grids, lev) for _ in 1:1, _ in 1:1, lev in level]

    return grid_quantity(coords) do state::StatePsiOmegaGridCNAB
        Γ_flat = quantities(state).Γ
        Γ = unflatten_circ(Γ_flat, prob.fluid.gridindex, level)
        return @. Γ / h^2
    end
end
