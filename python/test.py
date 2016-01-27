from pymnspline import Spline
import numpy as np
from numpy.testing import assert_almost_equal

des = np.array([0.952391, 0.607689 , -0.352613, -0.973527, -0.703281,
                0.213946 , 0.936819, 0.788606, -0.0478781, -0.34354])

x = np.arange(1.0, 11)
X = np.arange(1.5, 11)
X[-1] = 9.8
y = np.sin(x)

spl = Spline(x, y)

res = spl(X, False)
assert_almost_equal(res, des, decimal=6)
res = spl(X, True)
assert_almost_equal(res, des, decimal=6)

print "Test OK"
