#include <stdio.h>
#include <math.h>
#include "par.h"
#include "grid.h"
#include "test.h"
#include "ode.h"
#include "metric.h"
#include "geodesic.h"
#include "imageDisc.h"

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        printf("\nusage: geobygeo parfile.par\n");
        printf("exitting\n\n");
        return 0;
    }

    struct parList pars = PAR_DEFAULT;
    struct Grid grid = GRID_DEFAULT;
    read_pars(&pars, argv[1]);
    print_pars(&pars, NULL);

    int err = 0;
    err += setup_metric(&pars);
    err += setup_grid(&grid, &pars);

    if(err != 0)
    {
        printf("\nbad setup\n");
        printf("exitting\n\n");
        free_grid(&grid);
        return 0;
    }

    FILE *f = fopen("disc_im_grid.txt", "w");
    fclose(f);
    print_pars(&pars, "disc_im_grid.txt");

    double args[1];
    args[0] = 1.0;
    imageDiscGrid(&grid, args, "disc_im_grid.txt");

    free_grid(&grid);

    /*
    double c[4], n[3];
    double inclination = 75.0*M_PI/180.0;

    c[0] = 0.0;
    c[1] = 0.0;
    c[2] = 1.4e5 * cos(inclination);
    c[3] = 1.4e5 * sin(inclination);
    n[0] = 0.0;
    n[1] = 1.0 * cos(inclination);
    n[2] = 1.0 * sin(inclination);

    imageDisc(c, n, 100.0, 100.0, 10, 10, args, 
                "disc_im_90_schw_ks_10_10.txt");
    */
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
