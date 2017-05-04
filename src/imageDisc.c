#ifdef USEOMP
#include <omp.h>
#endif
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "par.h"
#include "grid.h"
#include "geodesic.h"
#include "metric.h"
#include "ode.h"
#include "imageDisc.h"

void imageDisc(double center[], double n[], double width, double height, 
                int NX, int NY, void *args, char filename[])
{
    /*
     * Tabulate the positions and 4-velocities of geodesics orginating
     * on the equatorial plane and terminating on a plane at position center[],
     * with normal n[], width, and height.  Geodesics are integrated backwards
     * from an nx by ny grid on the image plane.
     *
     * It is implicitly assumed the image plane is located in near-flat space 
     * and is static in the current coordinate system.  All geometrical
     * calculations at the image plane are performed in a cartesian sense.
     */

    int i, j, k;

    double *dat = (double *)malloc(16 * NX * NY * sizeof(double));

    double nx[3];
    double ny[3];
    double nz[3];
    
    double norm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);

    nz[0] = n[0]/norm;
    nz[1] = n[1]/norm;
    nz[2] = n[2]/norm;

    nx[0] = 1.0;
    nx[1] = nz[1] != 0.0 ? -nx[0]*nz[0]/nz[1] : 0.0;
    nx[2] = 0.0;
    norm = sqrt(nx[0]*nx[0] + nx[1]*nx[1] + nx[2]*nx[2]);
    nx[0] /= norm;
    nx[1] /= norm;
    nx[2] /= norm;

    ny[0] = nz[1]*nx[2] - nz[2]*nx[1];
    ny[1] = nz[2]*nx[0] - nz[0]*nx[2];
    ny[2] = nz[0]*nx[1] - nz[1]*nx[0];

    for(i=0; i<NX; i++)
    {
        double xi = (i*width)/(NX-1) - 0.5*width;

        for(j=0; j<NY; j++)
        {
            double yi = (j*height)/(NY-1) - 0.5*height;

            double xc[4], x[4], x1[4], u[4], uc[4], u1[4];
            
            xc[0] = center[0];
            uc[0] = 1.0;
            for (k=1; k<4; k++)
            {
                xc[k] = center[k] + xi*nx[k-1] + yi*ny[k-1];
                uc[k] = nz[k-1];
            }

            metric_cart2coord(xc, x, args);
            metric_vec2coordb(x, uc, u, args);

            double l = fabs(center[1])+fabs(center[2])+fabs(center[3]);
            double t0 = 0.0;
            double t1 = -10.0 * l;
            double dt0 = -l/1000.0;

            //printf("xc = %g %g %g %g\n", xc[0], xc[1], xc[2], xc[3]);
            //printf("x  = %g %g %g %g\n", x[0], x[1], x[2], x[3]);
            //printf("uc = %g %g %g %g\n", uc[0], uc[1], uc[2], uc[3]);
            //printf("u  = %g %g %g %g\n", u[0], u[1], u[2], u[3]);
            printf("%d %d\n", i, j);

            char rayname[80];
            sprintf(rayname, "ray_%d_%d.txt", i, j);
            geo_integrate_surface(x, u, x1, u1, t0, t1, dt0, 2, 0.5*M_PI, args,
                                    &dopr54, NULL);

            int ind = 16 * (NY*i + j);
            for(k=0; k<4; k++)
            {
                dat[ind+k] = x[k];
                dat[ind+k+4] = u[k];
                dat[ind+k+8] = x1[k];
                dat[ind+k+12] = u1[k];
            }
        }
    }

    FILE *f = fopen(filename, "w");
    for(i=0; i<NX; i++)
        for(j=0; j<NY; j++)
        {
            double xi = (i*width)/(NX-1) - 0.5*width;
            double yi = (j*height)/(NY-1) - 0.5*height;

            int ind = 16 * (NY*i + j);

            fprintf(f, "%d %d %.12lg %.12lg", i, j, xi, yi);
            for(k=0; k<16; k++)
                fprintf(f, " %.12lg", dat[ind+k]);
            fprintf(f, "\n");
        }
    fclose(f);

    free(dat);
}

void imageDiscGrid(struct Grid *g, void *args, char filename[])
{
    /*
     * Tabulate the positions and 4-velocities of geodesics orginating
     * on the equatorial plane and terminating on a plane at position center[],
     * with normal n[], width, and height.  Geodesics are integrated backwards
     * from an nx by ny grid on the image plane.
     *
     * It is implicitly assumed the image plane is located in near-flat space 
     * and is static in the current coordinate system.  All geometrical
     * calculations at the image plane are performed in a cartesian sense.
     */

    int N = g->N;

    double *dat = (double *)malloc(16 * N * sizeof(double));

    double R = g->distance;

    int i;
#ifdef USEOMP
    #pragma omp parallel for
#endif
    for(i=0; i<N; i++)
    {
        int k;
        double xc[4], x[4], x1[4], u[4], uc[4], u1[4];
        
        for(k=0; k<4; k++)
        {
            xc[k] = g->xc[4*i+k];
            uc[k] = g->nc[4*i+k];
        }

        metric_cart2coord(xc, x, args);
        metric_vec2coordb(x, uc, u, args);

        double t0 = 0.0;
        double t1 = -10.0 * R;
        double dt0 = -R/1000.0;

        printf("%d\n", i);

        char rayname[80];
        sprintf(rayname, "ray_%d_%d.txt", i, i);
        geo_integrate_surface(x, u, x1, u1, t0, t1, dt0, 2, 0.5*M_PI, args,
                                &dopr54, NULL);

        int ind = 16 * i;
        for(k=0; k<4; k++)
        {
            dat[ind+k] = x[k];
            dat[ind+k+4] = u[k];
            dat[ind+k+8] = x1[k];
            dat[ind+k+12] = u1[k];
        }
    }

    int N1 = g->N1;
    int N2 = g->N2;

    FILE *f = fopen(filename, "a");
    for(i=0; i<N; i++)
    {
        int k;
        double xi = g->xi[2*i+0];
        double yi = g->xi[2*i+1];

        int ind = 16 * i;

        fprintf(f, "%d %d %.12lg %.12lg", i/N2, i%N2, xi, yi);
        for(k=0; k<16; k++)
            fprintf(f, " %.12lg", dat[ind+k]);
        fprintf(f, "\n");
    }
    fclose(f);

    free(dat);
}
