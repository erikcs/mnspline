Note: this was written when Julia was in version 0.4, and has not been
updated since (so will most likely not work with the latest Julia version)
###### Installation
```
$git clone https://github.com/erikcs/mnspline
$cd julia
$make
julia> include("mnspline.jl")
...
julia> spline = mnspline(x, y)
julia> interpolated_values = spline(Xqueries)
```
