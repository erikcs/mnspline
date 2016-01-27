import numpy as np
cimport numpy as np

cdef extern from "../src/mnspline.h":
    int spline(const double *px, const double *py, size_t n, double *py2);

    int splint(const double *pxa, const double *pya, const double *py2a,
           size_t n, const double* px, double *py, size_t nx, int blookup);

cdef wspline(np.ndarray[np.double_t, ndim=1] x,
             np.ndarray[np.double_t, ndim=1] y):
    n = len(x)
    assert n == len(y)
    cdef np.ndarray[np.double_t, ndim=1] y2 = np.empty((n,), dtype=np.double)

    r = spline(<double*> x.data, <double*> y.data,
               <size_t> n, <double*> y2.data)
    if r != 0:
        raise RuntimeError("spline: malloc error")

    return y2

cdef wsplint(np.ndarray[np.double_t, ndim=1] x,
             np.ndarray[np.double_t, ndim=1] y,
             np.ndarray[np.double_t, ndim=1] y2,
             np.ndarray[np.double_t, ndim=1] X,
             int blookup):
    n = len(x)
    N = len(X)
    assert n == len(y) == len(y2)
    assert blookup == 0 or blookup == 1
    cdef np.ndarray[np.double_t, ndim=1] Y = np.empty((N,), dtype=np.double)

    splint(<double*> x.data, <double*> y.data, <double*> y2.data,
           <size_t> n, <double*> X.data, <double*> Y.data, <size_t> N,
           <int> blookup)

    return Y

cdef class Spline:
    """Estimates a (natural) cubic spline

    Parameters
    -----------
    x : np.ndarray
        the function evaluation points, x1 < x2 < ... < xN
    y : np.ndarray
        the function evaluated at above points x1, x2, ..., xN

    Methods
    -----------
    __call__(X, blookup=False)
    Return the interpolated function at query points X (ndarray)
    with lookup method linear (default: blookup=False),
    or bisection (blookup=True)
    """
    cdef np.ndarray x
    cdef np.ndarray y
    cdef np.ndarray y2

    def __init__(self, np.ndarray x, np.ndarray y):
        self.x = x
        self.y = y
        self.y2 = wspline(self.x, self.y)

    def __call__(self, np.ndarray X, int blookup=False):
        return wsplint(self.x, self.y, self.y2, X, blookup)

