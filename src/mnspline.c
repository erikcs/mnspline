/* mnspline.c - created January 2016
 * Natural cubic spline interpolation
 * 
 * based on "3.3 Cubic Spline Interpolation" in 
 * Numerical Recipes in C: The Art of Scientific Computing, 2nd edition.
 *
 * With additions: 
 *  performs the lookup on the array in parallel (with OpenMP)
 *  caches lookup results, with either a linear probe or bisection
 */

#include <stdlib.h>
#include <stdbool.h>
#include <omp.h>

/*#define DEBUG 1*/
#ifdef DEBUG
#include <stdio.h>
#endif

inline size_t
b_search(const double *pxa, double x, size_t idxlow, size_t idxhigh);

inline bool
lin_search(const double *pxa, double x, size_t idxlow, size_t idxhigh,
           size_t *index);

/* spline - calculate the 2nd. derivatives of the interpolating function
 * only needs to be called once - gives input to splint
 * the boundary condition at x1 and xN is zero.
 * Does no input validation (can be done in wrapper)
 *
 * Returns -1 on malloc failure.  
 *
 * px: Ptr to array of function evaluation points, x[1] < x[2] < ... < x[n] 
 * py: Ptr to array of function evaluated at above points
 * n: Size of array px and py
 * py2: Ptr to array, returns the 2nd derivatives of the interpolating function
 */
int
spline(const double *px, const double *py, size_t n, double *py2)
{
    double qn = 0.0; /* Upper boundary condition set to be 'natural' */
    double un = 0.0;  /* */
    double p;
    double sig;
    double *pu;

    if ( (pu = (double*) malloc((n-1) * sizeof(double))) == NULL )
        return -1;

    py2[0] = 0.0; /* Lower boundary condition set to be 'natural' */
    pu[0] = 0.0;  /* */

    /* The Tridiagonal algorithm */
    for (size_t i = 1; i < n - 1; i++)
    {
        sig     = (px[i] - px[i-1]) / (px[i+1] - px[i-1]);
        p       = sig * py2[i-1] + 2.0;
        py2[i]  = (sig - 1.0) / p;
        pu[i]   = (py[i+1] - py[i]) / (px[i+1] - px[i]) - 
                  (py[i] - py[i-1]) / (px[i] - px[i-1]);
        pu[i]   = (6.0 * pu[i] / (px[i+1] - px[i-1]) -
                  sig * pu[i-1]) / p; 
    }

    py2[n-1] = (un - qn * pu[n-2]) / (qn * py2[n-2] + 1.0);

    for (size_t k = n - 1; k-- > 0; )
        py2[k] = py2[k] * py2[k+1] + pu[k];

    free(pu);
    return 0;
}

/* splint - perform the interpolation
 * returns 0 - currently does no error checking (can be done in wrapper)
 *
 * pxa: Same input as to spline - px
 * pya: Same input as to spline - py
 * py2a: The output from spline - py2
 * n: Same input as to spline - n
 * px: Ptr to array of function evaluation points to be interpolated
 * py: Ptr to array, returns the interpolated function values evaluated at px
 * nx: Size of array px and py
 * blookup: 0=does linear probe on interpolation table,
 *          1=does binary search on interpolation table
 */
int
splint(const double *pxa, const double *pya, const double *py2a,
       size_t n, const double* px, double *py, size_t nx, int blookup)
{
    size_t klo;
    size_t khi;
    size_t prev_idx = 0;

    double h;
    double b;
    double a;
    double x;

    #pragma omp parallel for \
    firstprivate(prev_idx) \
    private(klo, khi, h, b, a, x)
    for (size_t i = 0; i < nx; i++)
    {
        x = px[i];

        /* LINEAR */
        if (blookup == 0) 
        {
            if ( lin_search(pxa, x, prev_idx,  n - 1, &klo) )
            {
                khi = klo + 1;
                prev_idx = klo;
            }
            else /* The linear probe did not find our x-val */
            {
                klo = b_search(pxa, x, 0, n - 1);
                khi = klo + 1;
                prev_idx = klo;
            }
        }
       /* BISECTION */
       else if (blookup == 1) 
       {
           if ( x >= pxa[prev_idx + 1] ) 
               prev_idx = b_search(pxa, x, prev_idx, n - 1);

           else if ( x < pxa[prev_idx] )
               prev_idx = b_search(pxa, x, 0, prev_idx);

        klo = prev_idx;
        khi = klo + 1;
       }

       /* COMPUTE */
        h     = pxa[khi] - pxa[klo];
        a     = (pxa[khi] - px[i]) / h;
        b     = (px[i] - pxa[klo]) / h;
        py[i] = a * pya[klo] + b * pya[khi] + 
            ((a*a*a -a) * py2a[klo] + (b*b*b - b) * py2a[khi]) * h*h / 6.0;
    }

    return 0;
}

/* Returns an index i s.t. pxa[i] <= x < pxa[i+1] */
inline size_t
b_search(const double *pxa, double x, size_t idxlow, size_t idxhigh)
{
    size_t ilo = idxlow;
    size_t ihi = idxhigh;
    size_t mid;
    while (ihi > ilo + 1)
    {
        mid = ilo + ((ihi - ilo) / 2);
        if (pxa[mid] > x)
            ihi = mid;
        else
            ilo = mid;
    }

    return ilo;
}

/* Linear probing for an i st. pxa[i] <= x < pxa[i+1] */
inline bool
lin_search(const double *pxa, double x, size_t idxlow, size_t idxhigh,
           size_t *index)
{
    bool found = false;

    for (size_t i = idxlow; i < idxhigh; i++)
    {
        if (pxa[i] <= x && pxa[i+1] > x)
        {
            found = true;
            *index = i;
            break;
        }
    }

    return found; 
}

#ifdef DEBUG
#include <stdio.h>
int main() {

    double x[10] = {1.0, 2.0, 3.0, 4.0, 5.0,
                    6.0, 7.0, 8.0, 9.0, 10.0};
    double xint[10] = {1.5, 2.5, 3.5, 4.5, 5.5,
                       6.5, 7.5, 8.5, 9.5, 9.8};
    double y[10] = {
        0.841471, 0.909297, 0.14112, -0.756802, -0.958924,
        -0.279415, 0.656987, 0.989358, 0.412118, -0.544021
    };

    double y2desired[10] = {
        0.0, -1.23211, -0.087569, 0.803919, 1.0467,
        0.299074, -0.701635, -1.11672, -0.289171, 0.0
    };
    double yintdesired[10] = {
        0.952391, 0.607689 , -0.352613, -0.973527, -0.703281,
        0.213946 , 0.936819, 0.788606, -0.0478781, -0.34354
    };

    double y2[10];
    double yint[10];

    spline(&x[0], &y[0], 10, &y2[0]);
    splint(&x[0], &y[0], &y2[0], 10,
            &xint[0], &yint[0], 10, 0);

    printf("%s", "linear probe\n");
    printf("y2\ty2(tst)\tyint\tyint(tst)\n");
    for (int i=0; i < 10; i++) 
        printf("%.3f\t%.3f\t%.3f\t%.3f\n",
                y2[i], y2desired[i], yint[i], yintdesired[i]);

    splint(&x[0], &y[0], &y2[0], 10,
            &xint[0], &yint[0], 10, 1);

    printf("%s", "bisection\n");
    printf("y2\ty2(tst)\tyint\tyint(tst)\n");
    for (int i=0; i < 10; i++) 
        printf("%.3f\t%.3f\t%.3f\t%.3f\n",
                y2[i], y2desired[i], yint[i], yintdesired[i]);
}
#endif

