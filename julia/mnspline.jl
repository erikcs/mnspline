# mnspline.jl

type Spline
    x::Vector{Float64}
    y::Vector{Float64}
    n::Int
    y2::Vector{Float64}
end

function mnspline(x::AbstractVector, y::AbstractVector)
    n = length(x)
    length(y) == n || error("First two arguments must be same length")

    xin = convert(Vector{Float64}, x)
    yin = convert(Vector{Float64}, y)

    y2 = Array(Float64, n)
        
    r = ccall( (:spline, "mnspline.so"), Int,
        (Ptr{Cdouble}, Ptr{Cdouble}, Csize_t, Ptr{Cdouble}),
        xin, yin, n, y2) == 0 || error("Terminated on malloc failure")

    return Spline(xin, yin, n, y2)
end

function evaluate(spline::Spline, X::AbstractVector)
    m = length(X)

    Xin = convert(Vector{Float64}, X)
    
    yint = Array(Float64, m)

    r = ccall( (:splint, "mnspline.so"), Int,
        (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Csize_t,
        Ptr{Cdouble}, Ptr{Cdouble}, Csize_t),
        spline.x, spline.y, spline.y2, spline.n,
        Xin, yint, m)
    
    return yint
end

Base.call(spline::Spline, X::AbstractVector) = evaluate(spline, X)

