struct ParameterizedCurve{F} <: Curve
    point::F
    arclength::Float64
    closed::Bool
end

function ParameterizedCurve(f)
    # f(t) -> (x(t), y(t))
    arclength = 0 # get arclength
    itp = Interpolation() # define interpolation

    function point(t)
        s = itp(t)
        return SVector(f(t))
    end

    closed = true # is first point approx equal last point?

    return ParameterizedCurve(point, arclength, closed)
end

isclosed(curve::ParameterizedCurve) = curve.closed
arclength(curve::ParameterizedCurve) = curve.arclength

(curve::ParameterizedCurve)(t) = curve.point(t)

struct Segments
    points::Vector{SVector{2,Float64}}
    lengths::Vector{Float64}
end

function partition(curve::ParameterizedCurve, n::Integer)
    # n is the number of points you should return

    # Initialize the points and lengths arrays
    points = Vector{SVector{2,Float64}}(undef, n)
    lengths = Vector{Float64}(undef, n)

    # Calculate the step size for parameter t based on whether the curve is closed or open
    step_size = curve.closed ? 1 / n : 1 / (n - 1)

    # Iterate through the range of t values and calculate points and lengths
    for i in 0:(n - 1)
        t = i * step_size
        points[i + 1] = curve.point(t)
        if i > 0
            lengths[i] = norm(points[i + 1] - points[i])
        end
    end

    # For a closed curve, set the length of the last segment
    if curve.closed
        lengths[end] = norm(points[end] - points[1])
    end

    return Segments(points, lengths)
end
