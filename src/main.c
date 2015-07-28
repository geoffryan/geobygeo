#include <stdio.h>
#include <math.h>
#include "test.h"
#include "ode.h"
#include "metric.h"
#include "geodesic.h"
#include "imageDisc.h"

int main(int argc, char *argv[])
{

    
    metric_ig  = &metric_ig_schw_ks;
    metric_dig = &metric_dig_schw_ks;
    metric_cart2coord = &metric_cart2coord_schw_ks;
    metric_vec2coordb = &metric_vec2coordb_schw_ks;
    metric_shadow = &metric_shadow_schw_ks;
    
/*   
    metric_ig  = &metric_ig_flat_sph;
    metric_dig = &metric_dig_flat_sph;
    metric_cart2coord = &metric_cart2coord_flat_sph;
    metric_vec2coordb = &metric_vec2coordb_flat_sph;
    metric_shadow = &metric_shadow_schw_ks;
*/


    
    double c[4], n[3];

    double inclination = M_PI/3.0;

    c[0] = 0.0;
    c[1] = 0.0;
    c[2] = 1.4e5 * cos(inclination);
    c[3] = 1.4e5 * sin(inclination);
    n[0] = 0.0;
    n[1] = 1.0 * cos(inclination);
    n[2] = 1.0 * sin(inclination);

    double args[1];
    args[0] = 1.0;

    imageDisc(c, n, 80.0, 80.0, 800, 800, args, "disc_im_60_schw_ks_800_800.txt");

    /*
    int n = 200;
    double x0[4], u0[4], x[4], u[4];
    double t0, t1;

    double r0 = 3.0;
    double r1 = 10.0;

    int i,j;
    int nr, np;
    nr = 20;
    np = 16;

    for(i=0; i<nr; i++)
    {
        double r = r0 + i*((r1-r0)/(nr-1));
        double ut = -0.1;
        double up = 0.2;
        double ur = 0.0;
        double u00 = sqrt(ur*ur + (ut*ut+up*up)/(r*r));
        for(j=0; j<np; j++)
        {
            x0[0] = 0.0;
            x0[1] = r;
            x0[2] = 0.5*M_PI;
            x0[3] = j*(2.0*M_PI/np);
            u0[0] = u00;
            u0[1] = ur;
            u0[2] = ut;
            u0[3] = up;

            char outname[128];
            sprintf(outname, "out_schw_%d_%d.txt", i, j);
            printf("Running %s...\n", outname);
            geo_integrate_generic(x0, u0, x, u, 0.0, 1000.0, 0.1, args, 
                                    &dopr54, outname);
        }

        
    }
    */

    return 0;
}
