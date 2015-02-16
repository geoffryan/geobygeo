#include <math.h>
#include "metric.h"

// Schwarzschild metric in spherical schwarzschild coordinates.
// x = (t, r, theta, phi)
// args = (G*M/c^2)

void  metric_ig_schw_sc(double *g, double *x, void *args)
{
    int i;
    double sint = sin(x[2]);
    double M = ((double *)args)[0];

    for(i=0; i<16; i++)
        g[i] = 0.0;

    g[0] = 1.0 / (-1.0 + 2*M/x[1]);
    g[5] = 1.0 - 2.0*M/x[1];
    g[10] = 1.0/(x[1]*x[1]);
    g[15] = 1.0/(x[1]*x[1]*sint*sint);
}

void metric_dig_schw_sc(double *dg, double *x, void *args)
{
    int i;
    for(i=0; i<64; i++)
        dg[i] = 0.0;
    
    double sint = sin(x[2]);
    double cost = cos(x[2]);
    double M = ((double *)args)[0];

    dg[16*1 +  0] = 2.0*M / ((2*M-x[1])*(2*M-x[1]));
    dg[16*1 +  5] = 2.0*M / (x[1]*x[1]);
    dg[16*1 + 10] = -2.0/(x[1]*x[1]*x[1]*sint*sint);
    dg[16*1 + 15] = -2.0/(x[1]*x[1]*x[1]*sint*sint);
    dg[16*2 + 15] = -2.0*cost/(x[1]*x[1]*sint*sint*sint);
}
