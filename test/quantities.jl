module TestQuantities

using ImmersedBodies.Quantities
using ImmersedBodies
using Debugger
using Test
# using Plots

@testset "quantities" begin

    # @testset "plotting" begin
    #     let xs = LinRange(0, 1, 3), ys = LinRange(0, 2, 4)
    #         z = xs .+ ys'
    #         val = GridValue(z, (xs, ys))

    #         p, = plot(val)
    #         @test p[1][:seriestype] == :heatmap # Make sure the recipe was applied
    #     end

    #     let coords = [(LinRange(-t, t, 6), LinRange(-0.5t, 0.5t, 3)) for t in 2.0 .^ (0:2)]
    #         z = cat((xs .+ ys' for (xs, ys) in coords)...; dims=3)
    #         val = MultiLevelGridValue(z, coords)

    #         p, = plot(val)
    #         @test p[1][:seriestype] == :heatmap # Make sure the recipe was applied
    #         @test length(p.series_list) == size(z, 3) # One heatmap per sublevel
    #     end
    # end

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
