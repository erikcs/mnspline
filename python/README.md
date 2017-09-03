###### Installation (OSX, with `gcc` from `brew`)
```
$git clone https://github.com/erikcs/mnspline
$cd python
$python setup.py build_ext --inplace # in current directory, or:
$python setup install # install on system

>>> from pymnspline import Spline
>>> itp = Spline(x, y)
>>> interpolated_values = itp(Xqueries, blookup=True/False)
```
