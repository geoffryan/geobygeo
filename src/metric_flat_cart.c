#include "metric.h"

// Flat space metric in cartesian coordinates.

void  metric_ig_flat_cart(double *g, double *x, void *args)
{
    int i;
    for(i=0; i<16; i++)
        g[i] = 0.0;
    g[0] = -1.0;
    g[5] = 1.0;
    g[10] = 1.0;
    g[15] = 1.0;
}

void metric_dig_flat_cart(double *dg, double *x, void *args)
{
    int i;
    for(i=0; i<64; i++)
        dg[i] = 0.0;
}
