###### Installation
```
$git clone https://github.com/nuffe/mnspline
$cd julia
$make
julia> include("mnspline.jl")
...
julia> spline = mnspline(x, y)
julia> interpolated_values = spline(Xqueries)
```
