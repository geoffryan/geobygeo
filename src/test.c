#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "test.h"

void oscillator(double t, double *x, void *args, double *xdot)
{
    //Harmonic Oscillator with unit frequency.

    xdot[0] = x[1];
    xdot[1] = -x[0];
}

void test_oscillator(double (*step)(double, double*, double *, int, void *,
                        void (*)(double,double*,void*,double*)), 
                        double t0, double t1, int n, char *filename)
{
    double t, dt, dt0;
    double x0[2], x1[2], x[2];
    x0[0] = cos(t0);    x0[1] = -sin(t0);
    x1[0] = cos(t1);    x1[1] = -sin(t1);

    dt = (t1-t0) / n;
    t = t0;
    x[0] = x0[0];
    x[1] = x0[1];

    FILE *f;
    if(filename != NULL)
    {
        f = fopen(filename, "w");
        fprintf(f,"%.12lg %.12lg %.12lg\n", t, x[0], x[1]);
    }

    while(t < t1)
    {
        dt0 = step(t, x, &dt, 2, NULL,  &oscillator);
        t += dt0;
        if(filename != NULL)
            fprintf(f,"%.12lg %.12lg %.12lg\n", t, x[0], x[1]);
    }
    if(filename != NULL)
        fclose(f);
    

    printf("Num:   %.12lf %.12lf\n", x[0], x[1]);
    printf("Exact: %.12lf %.12lf\n", x1[0], x1[1]);
    printf("Err: %.12lg\n", fabs(x[0]-x1[0])+fabs(x[1]-x1[1]));
}
