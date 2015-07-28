#include <math.h>
#include <stdio.h>
#include "metric.h"

// Schwarzschild metric in spherical kerr-schild coordinates.
// x = (t, r, theta, phi)
// args = (G*M/c^2)

void  metric_ig_schw_ks(double *g, double *x, void *args)
{
    int i;
    double sint = sin(x[2]);
    double M = ((double *)args)[0];

    for(i=0; i<16; i++)
        g[i] = 0.0;

    g[0] = -1.0 - 2*M/x[1];
    g[1] = 2*M/x[1];
    g[4] = 2*M/x[1];
    g[5] = 1.0 - 2.0*M/x[1];
    g[10] = 1.0/(x[1]*x[1]);
    g[15] = 1.0/(x[1]*x[1]*sint*sint);
}

void metric_dig_schw_ks(double *dg, double *x, void *args)
{
    int i;
    for(i=0; i<64; i++)
        dg[i] = 0.0;
    
    double sint = sin(x[2]);
    double cost = cos(x[2]);
    double M = ((double *)args)[0];

    dg[16*1 +  0] = 2.0*M / (x[1]*x[1]);
    dg[16*1 +  1] = -2.0*M / (x[1]*x[1]);
    dg[16*1 +  4] = -2.0*M / (x[1]*x[1]);
    dg[16*1 +  5] = 2.0*M / (x[1]*x[1]);
    dg[16*1 + 10] = -2.0/(x[1]*x[1]*x[1]);
    dg[16*1 + 15] = -2.0/(x[1]*x[1]*x[1]*sint*sint);
    dg[16*2 + 15] = -2.0*cost/(x[1]*x[1]*sint*sint*sint);
}

void metric_cart2coord_schw_ks(double *xc, double *x, void *args)
{
    // No factors of M show up, may not be right.
    
    x[0] = xc[0];
    x[1] = sqrt(xc[1]*xc[1] + xc[2]*xc[2] + xc[3]*xc[3]);
    x[2] = acos(xc[3]/x[1]);
    x[3] = atan2(xc[2], xc[1]);
}

void metric_vec2coordb_schw_ks(double *x, double *uc, double *u, void *args)
{
    double nr[3], nt[3], np[3], uo[4];
    double M = ((double *)args)[0];    
    double r = x[1];
    double cost = cos(x[2]);
    double sint = sin(x[2]);
    double cosp = cos(x[3]);
    double sinp = sin(x[3]);

    //Spherical Unit Vectors
    nr[0] = sint*cosp;
    nr[1] = sint*sinp;
    nr[2] = cost;
    nt[0] = cost*cosp;
    nt[1] = cost*sinp;
    nt[2] = -sint;
    np[0] = -sinp;
    np[1] = cosp;
    np[2] = 0.0;

    //Transform cartesian uc[] to orthonormal spherical basis.
    uo[0] = uc[0];
    uo[1] = nr[0]*uc[1] + nr[1]*uc[2] + nr[2]*uc[3];
    uo[2] = nt[0]*uc[1] + nt[1]*uc[2] + nt[2]*uc[3];
    uo[3] = np[0]*uc[1] + np[1]*uc[2] + np[2]*uc[3];

    // Transform from orthonormal to coordinate basis. u[] is covariant.
    u[0] = (-uo[0] + 2*M/r*uo[1]) / sqrt(1.0+2.0*M/r);
    u[1] = sqrt(1.0+2.0*M/r) * uo[1];
    u[2] = r * uo[2];
    u[3] = r*sint * uo[3];

    //double u2 = u[0]*u[0]*(-1-2*M/r) + 4*M/r*u[0]*u[1] + u[1]*u[1]*(1-2*M/r)
    //            + u[2]*u[2]/(r*r) + u[3]*u[3]/(r*r*sint*sint);
    //printf("norm: %g\n", u2);
}

int metric_shadow_schw_ks(double *x, void *args)
{
    double M = ((double*)args)[0];

    if(x[1] < 2*M + 1.0e-6)
        return 1;

    return 0;
}
