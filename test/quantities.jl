module TestQuantities

using ImmersedBodies.Quantities
using ImmersedBodies
# using Debugger
using Test
# using Plots

@testset "quantities" begin
    @testset "moving grid quantities" begin
        dx = 0.02
        xspan = (-1.0, 3.0)
        yspan = (-2.0, 2.0)
        basegrid = UniformGrid(dx, xspan, yspan)
        grids = MultiLevelGrid(basegrid, 5)

        scheme = default_scheme(basegrid; Umax=hypot(1, 2), cfl=0.2, safety=2.0)
        dt = timestep(scheme)

        flow_vec = (1.0, 2.0)

        Re = 100.0
        flow = FreestreamFlow(t -> flow_vec; Re=Re)
        flow_0 = FreestreamFlow(t -> (0.0, 0.0); Re=Re)
        horizontalflow = FreestreamFlow(t -> (hypot(2.0, 1.0), 0); Re=Re)

        # rotatingflow needs to account for the effects of changing velocity at all points.
        rotatingflow = FreestreamFlow(t -> hypot(2.0, 1.0) .* (cos(t), -sin(t)); Re=Re)

        movingframe = OffsetFrame(GlobalFrame()) do t
            r = (0.0, 0.0)
            v = .- flow_vec
            θ = 0.0
            Ω = 0.0
            return OffsetFrameInstant(r, v, θ, Ω)
        end

        rotatedframe = OffsetFrame(GlobalFrame()) do t
            r = (0.0, 0.0)
            v = (0.0, 0.0)
            θ = -atan(flow_vec[2] / flow_vec[1])
            Ω = 0.0
            return OffsetFrameInstant(r, v, θ, Ω)
        end

        rotatingframe = OffsetFrame(GlobalFrame()) do t
            r = (0.0, 0.0)
            v = (0.0, 0.0)
            θ = 0.0
            Ω = 1.0
            return OffsetFrameInstant(r, v, θ + Ω * t, Ω)
        end
        
        flowingfluid = PsiOmegaFluidGrid(flow, grids; scheme)
        nonflowingfluid = PsiOmegaFluidGrid(flow_0, grids; scheme, frame=movingframe)
        rotatedfluid = PsiOmegaFluidGrid(horizontalflow, grids; scheme, frame=rotatedframe)
        horizontalfluid = PsiOmegaFluidGrid(horizontalflow, grids; scheme, frame=rotatingframe)
        rotatingfluid = PsiOmegaFluidGrid(rotatingflow, grids; scheme)

        curve = Curves.Circle(0.5)
        stationarybody = RigidBody(partition(curve, flowingfluid))
        stationarybodies = BodyGroup([stationarybody])
        movingbody = RigidBody(partition(curve, nonflowingfluid))
        movingbodies = BodyGroup([movingbody])
        rotatedbody = RigidBody(partition(curve, rotatedfluid))
        rotatedbodies = BodyGroup([rotatedbody])
        rotatingbody = RigidBody(partition(curve, horizontalfluid))
        rotatingbodies = BodyGroup([rotatingbody])
        centerbody = RigidBody(partition(curve, rotatingfluid))
        centerbodies = BodyGroup([centerbody])

        flowprob = Problem(flowingfluid, stationarybodies)
        moveprob = Problem(nonflowingfluid, movingbodies)
        rotatedprob = Problem(rotatedfluid, rotatedbodies)
        rotatingprob = Problem(horizontalfluid, rotatingbodies)
        baselinerotatingprob = Problem(rotatingfluid, centerbodies)

        flowstate = initstate(flowprob)
        movestate = initstate(moveprob)
        rotatedstate = initstate(rotatedprob)
        rotatingstate = initstate(rotatingprob)
        baselinerotatingstate = initstate(baselinerotatingprob)

        solve!(flowstate, flowprob, 10 * dt)
        solve!(movestate, moveprob, 10 * dt)
        solve!(rotatedstate, rotatedprob, 10 * dt)
        solve!(rotatingstate, rotatingprob, 10 * dt)
        solve!(baselinerotatingstate, baselinerotatingprob, 10 * dt)

        basesf = Quantities.streamfunction(flowprob)(flowstate)

        # Test for unrotated motion
        @test basesf ≈ Quantities.streamfunction(moveprob)(movestate)
        # @run Quantities.streamfunction(rotatingprob)(rotatingstate)
        # Test for a stationary object at a given angular displacement
        @test Quantities.streamfunction(rotatedprob)(rotatedstate) ≈ basesf

        # Test for a stationary rotating object
        # This test will not pass. The rotating reference frame is treated as a rigid body, so the velocities of different points vary,
        # but the rotating flow has a constant velocity at all points. This is a feature, not a bug, and we need to make better tests.

        # TODO: Create a better way to test the streamfunction in a moving reference frame.
        @test Quantities.streamfunction(rotatingprob)(rotatingstate) ≈ Quantities.streamfunction(baselinerotatingprob)(baselinerotatingstate)
    end
end

end # module
