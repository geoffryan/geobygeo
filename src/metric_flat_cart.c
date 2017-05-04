#include "par.h"
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

void metric_cart2coord_flat_cart(double *xc, double *x, void *args)
{
    int i;
    for(i=0; i<4; i++)
        x[i] = xc[i];
}

void metric_vec2coordb_flat_cart(double *x, double *uc, double *u, void *args)
{
    int i;

    u[0] = -uc[0];
    for(i=1; i<4; i++)
        u[i] = uc[i];
}

int metric_shadow_flat_cart(double *x, void *args)
{
    return 0;
}

int metric_fix_domain_flat_cart(double *x, double *u, void *args)
{
    return 0;
}

