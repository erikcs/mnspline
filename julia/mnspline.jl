# mnspline.jl

nixlib = "mnspline.so"
winlib = "mnspline.dll"
const lib = @unix? nixlib : winlib

immutable Spline
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
        
    ccall( (:spline, lib), Int,
        (Ptr{Cdouble}, Ptr{Cdouble}, Csize_t, Ptr{Cdouble}),
        xin, yin, n, y2) == 0 || error("Terminated on malloc failure")

    return Spline(xin, yin, n, y2)
end

function evaluate(spline::Spline, X::AbstractVector, blookup::Bool=false)
    m = length(X)

    Xin = convert(Vector{Float64}, X)
    
    yint = Array(Float64, m)

    ccall( (:splint, lib), Int,
        (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Csize_t,
        Ptr{Cdouble}, Ptr{Cdouble}, Csize_t, Cint),
        spline.x, spline.y, spline.y2, spline.n,
        Xin, yint, m, blookup)
    
    return yint
end

Base.call(spline::Spline, X::AbstractVector; blookup::Bool=false) =
                                                evaluate(spline, X, blookup)

## INPLACE
function spline!(yint::Array{Float64}, y2::Array{Float64}, x::Array{Float64},
                 y::Array{Float64}, X::Array{Float64}, blookup::Int)
        n = length(x)
        m = length(X)

        ccall( (:spline, lib), Int,
        (Ptr{Cdouble}, Ptr{Cdouble}, Csize_t, Ptr{Cdouble}),
        x, y, n, y2)

        r = ccall( (:splint, lib), Int,
        (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Csize_t,
        Ptr{Cdouble}, Ptr{Cdouble}, Csize_t, Cint),
        x, y, y2, n,
        X, yint, m, blookup)
 
end

