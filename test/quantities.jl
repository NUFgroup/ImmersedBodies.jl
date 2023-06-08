module TestQuantities

using ImmersedBodies.Quantities
using ImmersedBodies
using Debugger
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
        flow_0 = FreestreamFlow(t -> (0.0, 0.0), Re=Re)

        frame = OffsetFrame(GlobalFrame()) do t
            r = (0.0, 0.0)
            v = .- flow_vec
            θ = 0.0
            Ω = 0.0
            return OffsetFrameInstant(r, v, θ, Ω)
        end
        
        flowingfluid = PsiOmegaFluidGrid(flow, grids; scheme)
        nonflowingfluid = PsiOmegaFluidGrid(flow_0, grids; scheme, frame)

        curve = Curves.Circle(0.5)
        stationarybody = RigidBody(partition(curve, flowingfluid))
        stationarybodies = BodyGroup([stationarybody])
        movingbody = RigidBody(partition(curve, nonflowingfluid))
        movingbodies = BodyGroup([movingbody])

        flowprob = Problem(flowingfluid, stationarybodies)
        moveprob = Problem(nonflowingfluid, movingbodies)

        flowstate = initstate(flowprob)
        movestate = initstate(moveprob)

        solve!(flowstate, flowprob, 10 * dt)
        solve!(movestate, moveprob, 10 * dt)

        # @run @test Quantities.streamfunction(flowprob)(flowstate) ≈ Quantities.streamfunction(moveprob)(movestate)
        @run Quantities.streamfunction(flowprob)(flowstate)
        @run Quantities.streamfunction(moveprob)(movestate)
    end
end

end # module
