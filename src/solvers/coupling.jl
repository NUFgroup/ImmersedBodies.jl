struct Reg
    basegrid::UniformGrid
    gridindex::PsiOmegaGridIndexing
    body_idx::Matrix{Int}
    supp_idx::UnitRange{Int}
    weight::Array{Float64,4}
    idx_offsets::Vector{Int} # Cumulative index of each body
end

struct RegT
    reg::Reg
end

function E_linearmap(reg::Reg)
    regT = RegT(reg)

    nf = 2 * size(reg.body_idx, 1)
    nq = reg.gridindex.nq

    return LinearMap(regT, reg, nf, nq)
end

function Reg(prob::Problem{<:PsiOmegaFluidGrid}, panels::Panels)
    nf = 2 * npanels(panels)
    body_idx = zeros(Int, size(panels.pos))

    # TODO: Add external option for supp
    supp = 6
    supp_idx = (-supp):supp
    weight = zeros(nf, 2, 2 * supp + 1, 2 * supp + 1)

    basegrid = baselevel(discretized(prob.fluid))
    gridindex = prob.fluid.gridindex

    nbodies = length(panels.perbody)
    idx_offsets = zeros(Int, nbodies)
    for i in 1:(nbodies - 1)
        idx_offsets[i + 1] = idx_offsets[i] + npanels(panels.perbody[i])
    end

    reg = Reg(basegrid, gridindex, body_idx, supp_idx, weight, idx_offsets)
    for (i, point) in enumerate(eachrow(panels.pos))
        update!(reg, point, i)
    end

    return reg
end

function update!(reg::Reg, panels::PanelView, bodyindex::Int)
    offset = reg.idx_offsets[bodyindex]
    indices = offset .+ (1:npanels(panels))
    for (i, point) in zip(indices, eachrow(panels.pos))
        update!(reg, point, i)
    end
end

function update!(reg::Reg, bodypoint::AbstractVector, index::Int)
    h = gridstep(reg.basegrid)
    nx = reg.gridindex.nx
    ny = reg.gridindex.ny

    x0, y0 = minimum.(xycoords(reg.basegrid))
    px, py = bodypoint

    # Nearest indices of body relative to grid
    ibody = floor(Int, (px - x0) / h)
    jbody = floor(Int, (py - y0) / h)

    if (
        !all(in(1:nx), ibody .+ extrema(reg.supp_idx)) ||
        !all(in(1:ny), jbody .+ extrema(reg.supp_idx))
    )
        error("Body outside innermost fluid grid")
    end

    reg.body_idx[index, :] .= (ibody, jbody)

    # Get regularized weight near IB points (u-vel points)
    x = @. x0 + h * (ibody - 1 + reg.supp_idx)
    y = permutedims(@. y0 + h * (jbody - 1 + reg.supp_idx))
    @. reg.weight[index, 1, :, :] = δh(x, px, h) * δh(y + h / 2, py, h)
    @. reg.weight[index, 2, :, :] = δh(x + h / 2, px, h) * δh(y, py, h)

    return reg
end

function (reg::Reg)(q_flat, fb_flat)
    # Matrix E'

    nb = size(reg.body_idx, 1)

    q_flat .= 0
    fb = reshape(fb_flat, nb, 2)
    qx, qy = split_flux(q_flat, reg.gridindex)

    for k in 1:nb
        i = reg.body_idx[k, 1] .+ reg.supp_idx
        j = reg.body_idx[k, 2] .+ reg.supp_idx
        @views @. qx[i, j] += reg.weight[k, 1, :, :] * fb[k, 1]
        @views @. qy[i, j] += reg.weight[k, 2, :, :] * fb[k, 2]

        # TODO: Throw proper exception or remove
        if !isfinite(sum(x -> x^2, qx[i, j]))
            error("infinite flux")
        end
    end

    return q_flat
end

function (regT::RegT)(fb_flat, q_flat)
    # Matrix E
    reg = regT.reg

    nb = size(reg.body_idx, 1)

    fb_flat .= 0
    fb = reshape(fb_flat, nb, 2)
    qx, qy = split_flux(q_flat, reg.gridindex)

    for k in 1:nb
        i = reg.body_idx[k, 1] .+ reg.supp_idx
        j = reg.body_idx[k, 2] .+ reg.supp_idx
        fb[k, 1] += @views dot(qx[i, j], reg.weight[k, 1, :, :])
        fb[k, 2] += @views dot(qy[i, j], reg.weight[k, 2, :, :])
    end

    return fb_flat
end

function δh(rf, rb, dr)
    # Discrete delta function used to relate flow to structure quantities

    # Take points on the flow domain (r) that are within the support (supp) of the IB points
    # (rb), and evaluate delta( abs(r - rb) )

    # Currently uses the Yang3 smooth delta function (see Yang et al, JCP, 2009), which has
    # a support of 6*h (3*h on each side)

    # Note that this gives slightly different answers than Fortran at around 1e-4,
    # apparently due to slight differences in the floating point arithmetic.  As far as I
    # can tell, this is what sets the bound on agreement between the two implementations.
    # It's possible this might be improved with arbitrary precision arithmetic (i.e.
    # BigFloats), but at least it doesn't seem to be a bug.

    # Note: the result is delta * h

    r = abs(rf - rb)
    r1 = r / dr
    r2 = r1 * r1
    r3 = r2 * r1
    r4 = r3 * r1

    return if (r1 <= 1.0)
        a5 = asin((1.0 / 2.0) * sqrt(3.0) * (2.0 * r1 - 1.0))
        a8 = sqrt(1.0 - 12.0 * r2 + 12.0 * r1)

        4.166666667e-2 * r4 +
        (-0.1388888889 + 3.472222222e-2 * a8) * r3 +
        (-7.121664902e-2 - 5.208333333e-2 * a8 + 0.2405626122 * a5) * r2 +
        (-0.2405626122 * a5 - 0.3792313933 + 0.1012731481 * a8) * r1 +
        8.0187537413e-2 * a5 - 4.195601852e-2 * a8 + 0.6485698427

    elseif (r1 <= 2.0)
        a6 = asin((1.0 / 2.0) * sqrt(3.0) * (-3.0 + 2.0 * r1))
        a9 = sqrt(-23.0 + 36.0 * r1 - 12.0 * r2)

        -6.250000000e-2 * r4 +
        (0.4861111111 - 1.736111111e-2 * a9) .* r3 +
        (-1.143175026 + 7.812500000e-2 * a9 - 0.1202813061 * a6) * r2 +
        (0.8751991178 + 0.3608439183 * a6 - 0.1548032407 * a9) * r1 - 0.2806563809 * a6 +
        8.22848104e-3 +
        0.1150173611 * a9

    elseif (r1 <= 3.0)
        a1 = asin((1.0 / 2.0 * (2.0 * r1 - 5.0)) * sqrt(3.0))
        a7 = sqrt(-71.0 - 12.0 * r2 + 60.0 * r1)

        2.083333333e-2 * r4 +
        (3.472222222e-3 * a7 - 0.2638888889) * r3 +
        (1.214391675 - 2.604166667e-2 * a7 + 2.405626122e-2 * a1) * r2 +
        (-0.1202813061 * a1 - 2.449273192 + 7.262731481e-2 * a7) * r1 +
        0.1523563211 * a1 +
        1.843201677 - 7.306134259e-2 * a7
    else
        0.0
    end
end

struct RedistributionWeights{ME}
    E::ME
    weights::Vector{Float64}
    qtmp::Vector{Float64}
    function RedistributionWeights(E, weights, qtmp)
        redist = new{typeof(E)}(E, weights, qtmp)
        update!(redist)
        return redist
    end
end

function RedistributionWeights(; E, qtmp)
    weights = similar(qtmp)
    return RedistributionWeights(E, weights, qtmp)
end

function W_linearmap(redist::RedistributionWeights)
    n = size(redist.E, 1)
    return LinearMap(redist, n)
end

function update!(redist::RedistributionWeights)
    (; weights, E) = redist

    nf = size(E, 1)
    mul!(weights, E', ones(nf))
    for (i, w) in pairs(weights)
        weights[i] = w < 1e-10 ? 0.0 : 1.0 / w
    end

    return nothing
end

function (redist::RedistributionWeights)(qout, qin)
    (; weights, qtmp, E) = redist

    mul!(qtmp, E', qin)
    qtmp .*= weights
    mul!(qout, E, qtmp)

    return nothing
end

function deforming_coord_unfilter(bodies::BodyGroup)
    nb_all = npanels(bodies)
    nb_deform = bodies.npanel_deform

    # z_all and z_deform are flattened [xs ys] matrices
    function deform_to_all(z_all, z_deform)
        z_all .= 0

        for index in bodies.deforming
            n = npanels(bodies[index.i_body])
            for j in 1:n
                i_deform = index.i_deform_panel + j
                i_all = index.i_panel + j

                z_all[i_all] = z_deform[i_deform]
                z_all[nb_all + i_all] = z_deform[nb_deform + i_deform]
            end
        end

        return z_all
    end

    return LinearMap(deform_to_all, 2 * nb_all, 2 * nb_deform)
end

function deforming_coord_filter(bodies::BodyGroup)
    nb_all = npanels(bodies)
    nb_deform = bodies.npanel_deform

    # z_all and z_deform are flattened [xs ys] matrices
    function all_to_deform(z_deform, z_all)
        for index in bodies.deforming
            n = npanels(bodies[index.i_body])
            for j in 1:n
                i_deform = index.i_deform_panel + j
                i_all = index.i_panel + j

                z_deform[i_deform] = z_all[i_all]
                z_deform[nb_deform + i_deform] = z_all[nb_all + i_all]
            end
        end

        return z_deform
    end

    return LinearMap(all_to_deform, 2 * nb_deform, 2 * nb_all)
end

abstract type SurfaceCoupler end

# Solve the Poisson equation (25) in Colonius & Taira (2008).
struct RigidSurfaceCoupler{MBinv,ME} <: SurfaceCoupler
    Binv::MBinv
    E::ME
    Ftmp::Vector{Float64}
    Q::Vector{Float64}
    h::Float64
end

function RigidSurfaceCoupler(; basegrid::UniformGrid, Binv, E, Ftmp, Q)
    h = gridstep(basegrid)
    return RigidSurfaceCoupler(Binv, E, Ftmp, Q, h)
end

function (coupler::RigidSurfaceCoupler)(state::StatePsiOmegaGridCNAB, qs)
    # Bodies moving in the grid frame
    # Solve the Poisson problem for bc2 = 0 (???) with nonzero boundary velocity ub
    # Bf̃ = Eq - ub
    #    = ECψ - ub

    (; F̃b, q0, panels) = quantities(state)
    (; Binv, E, Ftmp, Q, h) = coupler

    @views @. Q = qs[:, 1] + q0[:, 1]

    mul!(Ftmp, E, Q) # E*(qs .+ state.q0)

    ub = vec(panels.vel) # Flattened velocities

    # TODO: Is it worth dispatching to avoid this calculation for static bodies (ub = 0)
    @. Ftmp -= ub * h # Enforce no-slip conditions

    mul!(F̃b, Binv, Ftmp)

    return nothing
end

struct DeformingSurfaceCoupler{
    P<:Problem,O<:StackedStructureOps,R<:RedistributionWeights,ME,MB,MW,MF,MU,MS
} <: SurfaceCoupler
    prob::P
    qtot::Vector{Float64}
    ops::O
    reg::Reg
    redist::R
    E::ME
    B::MB
    W::MW
    filter_deform::MF
    unfilter_deform::MU
    force_to_traction::MS
    deform_k::DeformationState
end

function DeformingSurfaceCoupler(; prob, state, qtot, qtmp, reg, E, B)
    panels = bodypanels(state)

    ops = StackedStructureOps(prob, panels)
    init!(ops, prob.bodies, panels, state.qty.deform)

    redist = RedistributionWeights(; E, qtmp)

    W = W_linearmap(redist)
    filter_deform = deforming_coord_filter(prob.bodies)
    unfilter_deform = deforming_coord_unfilter(prob.bodies)

    nf = 2 * npanels(panels)
    h = (gridstep ∘ baselevel ∘ discretized)(prob.fluid)
    dt = timestep(prob)

    # Un-scales F̃b to be a physical traction
    force_to_traction = LinearMap(nf) do f_flat, F̃b_flat
        f = reshape(f_flat, :, 2)
        F̃b = reshape(F̃b_flat, :, 2)
        @. f = F̃b * h / (dt * panels.len)
        return f_flat
    end

    deform_k = similar(state.qty.deform)

    return DeformingSurfaceCoupler(
        prob,
        qtot,
        ops,
        reg,
        redist,
        E,
        B,
        W,
        filter_deform,
        unfilter_deform,
        force_to_traction,
        deform_k,
    )
end

function (coupler::DeformingSurfaceCoupler)(state::StatePsiOmegaGridCNAB, qs)
    (; F̃b, f_rdst, q0, panels, deform) = quantities(state)
    (;
        prob,
        qtot,
        ops,
        reg,
        redist,
        E,
        B,
        W,
        filter_deform,
        unfilter_deform,
        force_to_traction,
        deform_k,
    ) = coupler
    (; M, K, Khat_inv, stru_to_fluid_offset, fluid_to_stru_force) = ops
    # stru_to_fluid_offset: Equivalent of Itilde in previous versions
    # fluid_to_stru_force: Equivalent of Q * Itilde' in pervious versions

    dt = timestep(prob)
    h = (gridstep ∘ baselevel ∘ discretized)(prob.fluid)

    # Total flux; doesn't change over an FSI loop
    @views @. qtot = qs[:, 1] + q0[:, 1]

    # TODO: Add external option for tolerance
    tol_fsi = 1.e-5

    copy!(deform_k, deform)

    (; χ, ζ, ζdot) = deform
    χ_k = deform_k.χ
    ζ_k = deform_k.ζ
    ζdot_k = deform_k.ζdot

    # TODO: Define `W` and `force_to_traction` to only act on the deforming bodies so that
    # computation is not wasted when filtering the deforming bodies

    # References the operations in ops, so these are updated when ops is updated
    # This is QWx in the Fortran code
    B2 =
        unfilter_deform *
        stru_to_fluid_offset *
        Khat_inv *
        fluid_to_stru_force *
        filter_deform *
        W *
        force_to_traction *
        2 / dt * h
    Binv = Binv_linearmap(prob, B, B2)

    err_fsi = Inf
    while err_fsi > tol_fsi
        # Update structural operations (matrices) with current deformation
        update!(ops, prob.bodies, panels, deform_k)

        # Update the regularization (E matrix) for each deforming body
        for (; i_body) in prob.bodies.deforming
            update!(reg, panels.perbody[i_body], i_body)
        end

        # Update the redistribution weights
        update!(redist)

        # Develop RHS for linearized system for stresses
        F̃_kp1 = similar(F̃b)
        mul!(F̃_kp1, E, qtot)  # E*(qs + q0)... using fb here as working array

        # Enforce no-slip conditions for rigid bodies
        F̃_kp1_mat = reshape(F̃_kp1, :, 2)
        for (i, body) in enumerate(prob.bodies)
            if body isa RigidBody
                p = panels.perbody[i]
                n = size(p.vel, 1)
                @views @. F̃_kp1_mat[p.i_panel .+ (1:n), :] -= h * p.vel
            end
        end

        r_c = 2 / dt * (χ_k - χ) - ζ
        r_ζ = M * (ζdot + 4 / dt * ζ + 4 / dt^2 * (χ - χ_k)) - K * χ_k

        r_ζ = Khat_inv * r_ζ
        F_bg = -h * (2 / dt * r_ζ + r_c)
        F_sm = (unfilter_deform * stru_to_fluid_offset) * F_bg
        rhsf = F_sm + F̃_kp1

        mul!(F̃b, Binv, rhsf) # Solve for the stresses

        # Redistribute
        mul!(f_rdst, W * force_to_traction, F̃b)

        Δχ1 = Khat_inv * fluid_to_stru_force * filter_deform * f_rdst
        Δχ = r_ζ + Δχ1

        max_χ = maximum(abs, χ_k)
        max_Δχ = maximum(abs, Δχ)
        err_fsi = max_χ > 1e-13 ? max_Δχ / max_χ : max_Δχ

        # Update all structural quantities
        @. χ_k = χ_k + Δχ
        @. ζ_k = -ζ + 2 / dt * (χ_k - χ)
        @. ζdot_k = 4 / dt^2 * (χ_k - χ) - 4 / dt * ζ - ζdot

        update!(ops, prob.bodies, panels, deform_k)

        # Update the body panels that interface with the fluid
        for (op, i) in zip(ops.perbody, prob.bodies.deforming)
            body = prob.bodies[i.i_body]
            bodypanels = panels.perbody[i.i_body]

            structure_to_fluid_position(op)(vec(bodypanels.pos), χ_k)
            bodypanels.pos .+= reference_pos(body)

            # Recalculate segment lengths based on positions
            recalculate_lengths!(bodypanels)

            mul!(vec(bodypanels.vel), structure_to_fluid_offset(op), ζ_k)
        end
    end

    copy!(deform, deform_k)

    return nothing
end

function recalculate_lengths!(panels::PanelView)
    # TODO: Use a better way to recalculate the lengths.

    for i in 1:(npanels(panels) - 1)
        dx = panels.pos[i + 1, 1] - panels.pos[i, 1]
        dy = panels.pos[i + 1, 2] - panels.pos[i, 2]
        panels.len[i] = hypot(dx, dy)
    end

    # Keep the previous length, unless it is zero.
    if abs(panels.len[end]) < eps(eltype(panels.len))
        panels.len[end] = panels.len[end-1]
    end

    return nothing
end
