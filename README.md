# mnspline
multithreaded (natural) cubic spline interpolation

C implementation based on *Numerical Recipes in C*, with wrappers for ```Julia, Matlab, Python (todo)```

# Examples

###### gridsize  = 100 000
###### julia
```julia
n = 100000
nloops = div(10000000, n)
x = Float64[1:n;]
y = sin(x)
X = Float64[1.5:0.5:(n/2.0 + 1.0);]

splD = Dierckx.Spline1D(x, y; k=2)
spl = mnspline(x, y)

@time for i=1:nloops
    splD(X)
end

@time for i=1:nloops
    spl(X)
end

```

```
  0.255939 seconds (401 allocations: 76.311 MB, 2.95% gc time)
  0.174147 seconds (5.87 k allocations: 76.547 MB, 3.71% gc time)
```

###### matlab
```matlab
tic
for i=1:nloops
    interp1(x, y, X, 'spline');
end
toc

tic
for i=1:nloops
    mnspline(x, y, X);
end
toc

```

```
Elapsed time is 1.012905 seconds.
Elapsed time is 0.313663 seconds.
```

```
Platform Info:
  System: Darwin (x86_64-apple-darwin13.4.0)
  CPU: Intel(R) Core(TM) i7-5650U CPU @ 2.20GHz
```
