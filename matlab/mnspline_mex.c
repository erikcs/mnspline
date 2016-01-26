/*==========================================================
 * mnspline.c - (natural) cubic spline interpolation
 *
 * Interpolates a function y: Nx1
 * evaluated at x: Nx1 with x(1) < x(2) <... x(N) 
 *  (this condition is not checked in input validation)
 * at query points X: Mx1
 * and outputs Y: Mx1
 * with lookup method `blookup`: = 0 linear probe
 *                               = 1 bisection
 *
 * The calling syntax is:
 *
 *		Y = mnspline(x, y, X, blookup)
 *
 *========================================================*/

#include "mex.h"
#include "../src/mnspline.h"
#include <stdlib.h>

/*  Inefficient if the matlab function is called with
 *  the same (x, y, .) arguments several times */
int mwrapper(const double *x, const double *y, size_t n,
             const double *X, double *Y, size_t N, int blookup)
{  
    double *y2;
    if ( (y2 = (double*) malloc(n * sizeof(double))) == NULL )
        return -1;

    spline(x, y, n, y2);
    splint(x, y, y2, n,
           X, Y, N, blookup);

    free(y2); 

    return 0;
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    /* check for proper number of arguments */
    if ( nrhs!=4 )
        mexErrMsgIdAndTxt("mnspline:nrhs", "Four inputs required.");

    /* get dimension of the inputs */
    int n = mxGetM(prhs[0]);
    int N = mxGetM(prhs[2]);
    if ( n != mxGetM(prhs[1]) )
        mexErrMsgIdAndTxt("mnspline:xxx", 
                "First two inputs needs to be of same dimension");

    /* create output array and get a ptr to the data in the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize) N, 1, mxREAL);
    double *Y = mxGetPr(plhs[0]);

    /* create a pointer to the real data in the input matrix  */
    double *x  = mxGetPr(prhs[0]);
    double *y  = mxGetPr(prhs[1]);
    double *X  = mxGetPr(prhs[2]);

    /* get the scalar input value */
    int blookup = mxGetScalar(prhs[3]);

    int r =  mwrapper(x, y, n,
                      X, Y, N, blookup);

    if (r != 0)
        mexErrMsgIdAndTxt("mnspline:xxx", "Terminated on malloc failure");
}

