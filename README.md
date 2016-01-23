# mnspline
multithreaded (natural) cubic spline interpolation

C implementation based on *Numerical Recipes in C*, with wrappers for ```Julia, Matlab, Python (todo)```

# Examples

###### gridsize  = 100 000 (monotonically increasing)
###### julia
called with linear probing as the lookup method (```blookup=0```)
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
  0.303440 seconds (401 allocations: 76.311 MB, 13.35% gc time)
  0.080551 seconds (301 allocations: 76.303 MB, 8.49% gc time)
```
###### gridsize  = 100 000 (unordered)
called with bisection as the lookup method (```blookup=1```)

```julia
...
shuffle!(X)
@time for i=1:nloops
    splD(X)
end

@time for i=1:nloops
    spl(X, blookup=true)
end
```
```
111.561258 seconds (401 allocations: 76.311 MB, 0.01% gc time)
  0.744946 seconds (401 allocations: 76.311 MB, 0.88% gc time)
```





###### matlab (gridsize  = 100 000, monotonically increasing, ```blookup=0```)
```matlab
tic
for i=1:nloops
    interp1(x, y, X, 'spline');
end
toc

tic
for i=1:nloops
    mnspline(x, y, X, 0);
end
toc

```

```
Elapsed time is 1.120193 seconds.
Elapsed time is 0.183580 seconds.
```

```
Platform Info:
  System: Darwin (x86_64-apple-darwin13.4.0)
  CPU: Intel(R) Core(TM) i7-5650U CPU @ 2.20GHz
```
