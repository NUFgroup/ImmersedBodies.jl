module TestQuantities

using ImmersedBodies.Quantities
using ImmersedBodies.Bodies
using ImmersedBodies
using Test
# using Plots
# Plots was used here for other tests, but those tests broke the repository.
# May be worth looking into fixing the Plots library for the repo.

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

        movingframe = OffsetFrame(GlobalFrame()) do t
            r = (0.0, 0.0)
            v = .-flow_vec
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
        horizontalfluid = PsiOmegaFluidGrid(
            horizontalflow, grids; scheme, frame=rotatingframe
        )

        curve = Curves.Circle(0.5)
        stationarybody = RigidBody(partition(curve, flowingfluid))
        stationarybodies = BodyGroup([stationarybody])
        movingbody = RigidBody(partition(curve, nonflowingfluid))
        movingbodies = BodyGroup([movingbody])
        rotatedbody = RigidBody(partition(curve, rotatedfluid))
        rotatedbodies = BodyGroup([rotatedbody])
        rotatingbody = RigidBody(partition(curve, horizontalfluid))
        rotatingbodies = BodyGroup([rotatingbody])

        flowprob = Problem(flowingfluid, stationarybodies)
        moveprob = Problem(nonflowingfluid, movingbodies)
        rotatedprob = Problem(rotatedfluid, rotatedbodies)
        rotatingprob = Problem(horizontalfluid, rotatingbodies)

        flowstate = initstate(flowprob)
        movestate = initstate(moveprob)
        rotatedstate = initstate(rotatedprob)
        rotatingstate = initstate(rotatingprob)

        solve!(flowstate, flowprob, 10 * dt)
        solve!(movestate, moveprob, 10 * dt)
        solve!(rotatedstate, rotatedprob, 10 * dt)
        solve!(rotatingstate, rotatingprob, 10 * dt)

        basesf = Quantities.streamfunction(flowprob)(flowstate)

        # Test for unrotated motion
        @test basesf ≈ Quantities.streamfunction(moveprob)(movestate)

        # Test for a stationary object at a given angular displacement
        @test Quantities.streamfunction(rotatedprob)(rotatedstate) ≈ basesf

        smallcylinder = Curves.Circle(0.02)

        # Offset the frame by 10 * dt radians so that it is at 0 radians in 10 timesteps.
        # This makes it easy to calculate the expected stream function.
        offsetrotatingframe = OffsetFrame(GlobalFrame()) do t
            r = (0.0, 0.0)
            v = (0.0, 0.0)
            θ = -10.0 * dt
            Ω = 1.0
            return OffsetFrameInstant(r, v, θ + t * Ω, Ω)
        end

        relativeflow = (hypot(2.0, 1.0) + 1.98, 0.98)

        # Expected stream function at (-0.98, 1.98)
        expectedsf = -relativeflow[1] * -0.98 + relativeflow[2] * 1.98

        thisfluid = PsiOmegaFluidGrid(
            horizontalflow, grids; scheme, frame=offsetrotatingframe
        )
        smallcylinderbody = RigidBody(partition(smallcylinder, thisfluid))
        smallcylindergroup = BodyGroup([smallcylinderbody])
        smallcylinderprob = Problem(thisfluid, smallcylindergroup)

        smallstate = initstate(smallcylinderprob)
        solve!(smallstate, smallcylinderprob, 10 * dt)

        # Access streamfunction at (-0.98, 1.98) on the lowest grid level.
        sf = Quantities.streamfunction(smallcylinderprob)(smallstate)[1, 199, 1]

        @test isapprox(sf, expectedsf, rtol=0.01)
    end
end

end # module
