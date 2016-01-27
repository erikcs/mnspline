using Base.Test

include("mnspline.jl")

des = [0.952391, 0.607689 , -0.352613, -0.973527, -0.703281,
       0.213946 , 0.936819, 0.788606, -0.0478781, -0.34354]

x = 1:10.0
y = sin(x)
X = collect(1.5:10.5)
X[10] = 9.8

spl = mnspline(x, y)
res = spl(X)
@test_approx_eq_eps res des 1e-6
res = spl(X, blookup=true)
@test_approx_eq_eps res des 1e-6
println("Test OK")

