# mnspline
multithreaded (natural) cubic spline interpolation

C implementation with wrappers for ```Julia, Matlab, Python```

# Examples

###### gridsize  = 100 000 (increasing)
###### julia
called with linear probing as the lookup method (```blookup=0```, default)
```julia
n = 100000
nloops = div(10000000, n)
x = Float64[1:n;]
y = sin(x)
X = Float64[1.5:0.5:(n/2.0 + 1.0);]

splD = Dierckx.Spline1D(x, y)
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
the bulk of the time is spent solving the tridiagonal system (which is still serial):

```julia
@time for i=1:nloops
    splD = Dierckx.Spline1D(x, y)
    splD(X)
end

@time for i=1:nloops
    spl = mnspline(x, y)
    spl(X, blookup=false)
end
```

```
1.934707 seconds (1.80 k allocations: 1.826 GB, 5.68% gc time)
0.215632 seconds (701 allocations: 152.615 MB, 8.22% gc time)
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

###### matlab (gridsize  = 100 000, increasing)
(```blookup=0```)
```matlab
n = 100000;
nloops = floor(10000000 / n);
x = (1:n)';
y = sin(x);
X = (1.5:0.5:(n/2 + 1))';

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

###### python (gridsize  = 100 000, increasing)
```python
n = 100000
nloops = 10000000 / n
x = np.arange(1.0, n + 1)
y = np.sin(x)
X = np.arange(1.5, (n/2.0 + 1.0) + 0.5, 0.5)

s = time.time()
for n in range(nloops):
    spl = Spline(x, y)
    spl(X)
print 'elapsed time (secs): ', time.time() - s
```

```
elapsed time (secs):  0.178623199463
```

```
Platform Info:
  System: Darwin (x86_64-apple-darwin13.4.0)
  CPU: Intel(R) Core(TM) i7-5650U CPU @ 2.20GHz
```
