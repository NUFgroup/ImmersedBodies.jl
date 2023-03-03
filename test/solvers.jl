module TestSolvers

using ImmersedBodies
using ImmersedBodies.Solvers
using Test

@testset "solvers" begin
    dx = 0.02
    xspan = (-1.0, 3.0)
    yspan = (-2.0, 2.0)
    basegrid = UniformGrid(dx, xspan, yspan)
    grids = MultiLevelGrid(basegrid, 5)

    @test gridstep(basegrid) ≈ dx
    @test gridstep(baselevel(grids)) ≈ dx
    @test gridstep(grids, 3) ≈ 4 * dx

    dt = 0.004
    scheme = CNAB(dt)

    @test timestep(scheme) ≈ dt

    freestream(t) = (1.0, -0.1)
    Re = 100.0

    function max_vel(state)
        qty = state.qty
        nq, nlev = size(qty.q)
        return maximum(1:nlev) do lev
            h = gridstep(grids, lev)
            maximum(1:nq) do i
                abs((qty.q[i, lev] + qty.q0[i, lev]) / h)
            end
        end
    end

    @test_throws "invalid reference frame" PsiOmegaFluidGrid(
        scheme, grids, freestream, DiscretizationFrame(); Re
    )

    let
        fluid = PsiOmegaFluidGrid(scheme, grids, freestream; Re)

        # Body outside of innermost fluid grid
        curve = Curves.LineSegment((0.0, 1.5), (0.0, 2.1))
        body = RigidBody(fluid, curve)
        bodies = BodyGroup([body])

        prob = Problem(fluid, bodies)

        @test_throws "outside innermost fluid grid" solve(prob, (0.0, 5 * dt))
    end

    @testset "static rigid bodies" begin
        fluid = PsiOmegaFluidGrid(scheme, grids, freestream; Re)

        curve = Curves.Circle(0.5)
        body = RigidBody(fluid, curve)
        bodies = BodyGroup([body])

        prob = Problem(fluid, bodies)
        state = initstate(prob)

        solve!(state, prob, 10 * dt)
        @test max_vel(state) < 20
    end

    @testset "moving grid" begin
        frame = OffsetFrame(GlobalFrame()) do t
            r = 0.1 .* (cos(t), sin(t))
            v = 0.1 .* (-sin(t), cos(t))
            θ = 0.2 * sin(t)
            Ω = 0.2 * cos(t)
            return OffsetFrameInstant(r, v, θ, Ω)
        end

        fluid = PsiOmegaFluidGrid(scheme, grids, freestream, frame; Re)

        curve = Curves.Circle(0.5)
        body = RigidBody(fluid, curve)
        bodies = BodyGroup([body])

        prob = Problem(fluid, bodies)
        state = initstate(prob)

        solve!(state, prob, 10 * dt)
        @test max_vel(state) < 20
    end

    @testset "moving rigid bodies" begin
        fluid = PsiOmegaFluidGrid(scheme, grids, freestream; Re)

        function offset(t)
            r = 0.1 .* (cos(t), sin(t))
            v = 0.1 .* (-sin(t), cos(t))
            θ = 0.2 * sin(t)
            Ω = 0.2 * cos(t)
            return OffsetFrameInstant(r, v, θ, Ω)
        end
        body1 = RigidBody(
            fluid,
            Curves.LineSegment((-0.3, 0.5), (0.3, 0.5)),
            OffsetFrame(offset, DiscretizationFrame()),
        )
        body2 = RigidBody(
            fluid,
            Curves.LineSegment((-0.3, -0.5), (0.3, -0.5)),
            OffsetFrame(offset, GlobalFrame()),
        )
        bodies = BodyGroup([body1, body2])

        prob = Problem(fluid, bodies)
        state = initstate(prob)

        solve!(state, prob, 10 * dt)
        @test max_vel(state) < 20
    end
end

end # module