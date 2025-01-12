#include "matrix.hpp"
#include "cross.hpp"
#include <stdio.h>
#include <sys/time.h>
#include <math.h>

#ifndef CALCNORMS
#define CALCNORMS 0
#endif

double f(unsigned i, unsigned j, double h)
{
    return 1.0 / (h * (double)(i + 1) + h * (double)(j + 1));
    //return pow(h * (double)(i + 1) * h * (double)(j), 1.0 / 6.0);
    /*if(i >= j) {
        return exp(-3.0 * h * double(i));
    } else {
        return 0.0;
    }*/
    /* double tmp = (i + 1) / ((double)(j) + 1);
    return pow(tmp, 1.0 / 3) + pow(tmp, - 1.0 / 3); */
}

/* calc differnce in microseconds of timevalues */
double time_diff(struct timeval *end, struct timeval *start)
{
    return (end->tv_sec - start->tv_sec) + 
        (end->tv_usec - start->tv_usec) * 1e-6;
}

int main()
{
    unsigned r, m = 10, n = 10;
#if CALCNORMS
    unsigned i, j;
    double tmp, fn, cn;
#endif
    double h = 0.01;
    struct timeval tv1, tv2;
    FunctionMatrix A(m, n, h, f);
    Cross A_approx(A);
    double eps;
    /*for(i = 0; i < m; ++i) {
        for(j = 0; j < n; ++j) {
            printf("[%d][%d]=%f ", i, j, A.elem(i, j));
        }
        putchar('\n');
    }*/
    for(eps = 1; eps >= 1e-6; eps /= 10) {
        printf("eps: %g\n", eps);
        gettimeofday(&tv1, NULL);
        A_approx.Approximate(eps);
        gettimeofday(&tv2, NULL);
        r = A_approx.GetRank();
        printf("approx-rank -> %d\ntime to calc: %f\n", r,
                time_diff(&tv2, &tv1));
#if CALCNORMS
        fn = 0;
        cn = 0;
        for(i = 0; i < m; ++i) {
            for(j = 0; j < n; ++j) {
                tmp = fabs(A.elem(i,j) - A_approx.elem(i, j));
                fn += tmp * tmp;
                if(tmp > cn) {
                    cn = tmp;
                }
            }
        }
        fn = sqrt(fn);
        printf("F: %g\nC: %g\n", fn, cn);
#endif
    }
    return 0;
}
