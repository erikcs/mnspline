###### Installation (OSX)
```
$git clone https://github.com/nuffe/mnspline
$cd matlab
$./mnspline_mexbuild.sh

matlab> blookup = 0; % =1 for bisection
matlab> interpolated_values = mnspline(x, y, Xqueries, blookup)
```
