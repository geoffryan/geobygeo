#include <stdlib.h>
#include "ode.h"
#include "test.h"

void forward_euler(double t, double *x, double *dt, int n, void *args,
                    void (*xdot)(double,double*,void*,double*))
{
    // Simple forward euler integration scheme.  Updates 'x' in place.

    int i;
    double *k = (double *) malloc(n * sizeof(double));
    
    xdot(t, x, args, k);

    for(i=0; i<n; i++)
        x[i] += (*dt)*k[i];

    free(k);
}

void rk2(double t, double *x, double *dt, int n, void *args,
            void (*xdot)(double,double*,void*,double*))
{
    // Basic Runge-Kutta 2nd order integration scheme.  Updates 'x' in place.

    int i;
    double *k1 = (double *) malloc(n * sizeof(double));
    double *k2 = (double *) malloc(n * sizeof(double));
    double *x1 = (double *) malloc(n * sizeof(double));
    
    xdot(t, x, args, k1);

    for(i=0; i<n; i++)
        x1[i] = x[i] + 0.5*(*dt)*k1[i];
    xdot(t, x1, args, k2);

    for(i=0; i<n; i++)
        x[i] += (*dt)*k2[i];

    free(k1);
    free(k2);
    free(x1);
}

void rk4(double t, double *x, double *dt, int n, void *args,
            void (*xdot)(double,double*,void*,double*))
{
    // Basic Runge-Kutta 4th order integration scheme.  Updates 'x' in place.

    int i;
    double *k1 = (double *) malloc(n * sizeof(double));
    double *k2 = (double *) malloc(n * sizeof(double));
    double *k3 = (double *) malloc(n * sizeof(double));
    double *k4 = (double *) malloc(n * sizeof(double));
    double *x1 = (double *) malloc(n * sizeof(double));
    double *x2 = (double *) malloc(n * sizeof(double));
    double *x3 = (double *) malloc(n * sizeof(double));
    
    xdot(t, x, args, k1);

    for(i=0; i<n; i++)
        x1[i] = x[i] + 0.5*(*dt)*k1[i];
    xdot(t, x1, args, k2);

    for(i=0; i<n; i++)
        x2[i] = x[i] + 0.5*(*dt)*k2[i];
    xdot(t, x2, args, k3);

    for(i=0; i<n; i++)
        x3[i] = x[i] + (*dt)*k3[i];
    xdot(t, x3, args, k4);

    for(i=0; i<n; i++)
        x[i] += (*dt)*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;

    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(x1);
    free(x2);
    free(x3);
}
