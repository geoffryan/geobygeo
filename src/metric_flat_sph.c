#include <math.h>
#include "metric.h"

// Flat space metric in spherical polar coordinates.
// x = (t, r, theta, phi)

void  metric_ig_flat_sph(double *g, double *x, void *args)
{
    int i;
    double sint = sin(x[2]);
    for(i=0; i<16; i++)
        g[i] = 0.0;
    g[0] = -1.0;
    g[5] = 1.0;
    g[10] = 1.0/(x[1]*x[1]);
    g[15] = 1.0/(x[1]*x[1]*sint*sint);
}

void metric_dig_flat_sph(double *dg, double *x, void *args)
{
    int i;
    for(i=0; i<64; i++)
        dg[i] = 0.0;
    
    double sint = sin(x[2]);
    double cost = cos(x[2]);

    dg[16*1 + 10] = -2.0/(x[1]*x[1]*x[1]);
    dg[16*1 + 15] = -2.0/(x[1]*x[1]*x[1]*sint*sint);
    dg[16*2 + 15] = -2.0*cost/(x[1]*x[1]*sint*sint*sint);
}
