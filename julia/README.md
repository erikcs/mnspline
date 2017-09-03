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
