#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "par.h"
#include "metric.h"
#include "grid.h"

int setup_grid(struct Grid *g, struct parList *par)
{
    int choice = par->grid;
    int err = 0;

    if(choice == 0)
    {
        grid_initialize(g, par);
        grid_generate_cart(g);
    }
    else if(choice == 1)
    {
        grid_initialize(g, par);
        grid_generate_elps(g);
    }
    else
    {
        printf("Bad Grid choice: %d\n", choice);
        err = 1;
    }

    return err;
}

void grid_initialize(struct Grid *g, struct parList *par)
{
    int N1 = par->N1;
    int N2 = par->N2;

    double X1a = par->X1a;
    double X1b = par->X1b;
    double X2a = par->X2a;
    double X2b = par->X2b;

    double distance = par->distance;
    double inclination = M_PI * par->inclination / 180.0;
    double azimuth = M_PI * par->azimuth / 180.0;

    if(N1 < 1)
        N1 = 1;
    if(N2 < 1)
        N2 = 1;

    int N = N1*N2;

    g->xi = (double *) malloc(N*2*sizeof(double));
    g->xc = (double *) malloc(N*4*sizeof(double));
    g->nc = (double *) malloc(N*4*sizeof(double));

    g->mode = par->grid;
    g->N1 = N1;
    g->N2 = N2;
    g->N = N;

    g->X1a = X1a;
    g->X1b = X1b;
    g->X2a = X2a;
    g->X2b = X2b;

    g->distance = distance;
    g->inclination = inclination;
    g->azimuth = azimuth;
}

void grid_generate_cart(struct Grid *g)
{
    int Nx = g->N1;
    int Ny = g->N2;

    double R = g->distance;
    double inc = g->inclination;
    double azi = g->azimuth;

    double xa = g->X1a;
    double xb = g->X1b;
    double ya = g->X2a;
    double yb = g->X2b;

    double c[3], nx[3], ny[3], nz[3];
    double cosi = cos(inc);
    double sini = sin(inc);
    double cosa = cos(azi);
    double sina = sin(azi);

    c[0] = R * sin(inc) * cos(azi);
    c[1] = R * sin(inc) * sin(azi);
    c[2] = R * cos(inc);

    nz[0] = sini * cosa;
    nz[1] = sini * sina;
    nz[2] = cosi;

    nx[0] = -sina;
    nx[1] = cosa;
    nx[2] = 0.0;
    
    ny[0] = -cosi * cosa;
    ny[1] = -cosi * sina;
    ny[2] = sini;

    double dx = Nx > 1 ? (xb-xa) / (Nx - 1) : 0.0;
    double dy = Ny > 1 ? (yb-ya) / (Ny - 1) : 0.0;

    double *xi = g->xi;
    double *xc = g->xc;
    double *nc = g->nc;

    int i,j;
    for(i=0; i<Nx; i++)
    {
        for(j=0; j<Ny; j++)
        {
            double X = xa + i*dx;
            double Y = ya + j*dy;
            xi[2*(Ny*i+j) + 0] = X;
            xi[2*(Ny*i+j) + 1] = Y;

            printf("%g, %g\n", X, Y);

            xc[4*(Ny*i+j)+0] = 0.0;
            xc[4*(Ny*i+j)+1] = c[0] + X*nx[0] + Y*ny[0];
            xc[4*(Ny*i+j)+2] = c[1] + X*nx[1] + Y*ny[1];
            xc[4*(Ny*i+j)+3] = c[2] + X*nx[2] + Y*ny[2];
            
            nc[4*(Ny*i+j)+0] = 1.0;
            nc[4*(Ny*i+j)+1] = nz[0];
            nc[4*(Ny*i+j)+2] = nz[1];
            nc[4*(Ny*i+j)+3] = nz[2];
        }
    }
}

void grid_generate_elps(struct Grid *g)
{
    int Nr = g->N1;
    int Np = g->N2;

    double R = g->distance;
    double inc = g->inclination;
    double azi = g->azimuth;

    double ra = g->X1a;
    double rb = g->X1b;

    double c[3], nx[3], ny[3], nz[3];
    double cosi = cos(inc);
    double sini = sin(inc);
    double cosa = cos(azi);
    double sina = sin(azi);

    c[0] = R * sin(inc) * cos(azi);
    c[1] = R * sin(inc) * sin(azi);
    c[2] = R * cos(inc);

    nz[0] = sini * cosa;
    nz[1] = sini * sina;
    nz[2] = cosi;

    nx[0] = -sina;
    nx[1] = cosa;
    nx[2] = 0.0;
    
    ny[0] = -cosi * cosa;
    ny[1] = -cosi * sina;
    ny[2] = sini;

    double dphi = 2*M_PI / Np;

    double *xi = g->xi;
    double *xc = g->xc;
    double *nc = g->nc;

    int i,j;
    double r = ra;

    for(i=0; i<Nr; i++)
    {
        if(Nr > 1)
            r = ra * pow(rb/ra, ((double) i) / ((double)(Nr-1)));

        for(j=0; j<Np; j++)
        {
            double phi = j*dphi;
            double X = r * cos(phi);
            double Y = r * sin(phi) * cosi;

            xi[2*(Np*i+j) + 0] = r;
            xi[2*(Np*i+j) + 1] = phi;

            //printf("%g, %g, %g, %g\n", X, Y, r, phi);

            xc[4*(Np*i+j)+0] = 0.0;
            xc[4*(Np*i+j)+1] = c[0] + X*nx[0] + Y*ny[0];
            xc[4*(Np*i+j)+2] = c[1] + X*nx[1] + Y*ny[1];
            xc[4*(Np*i+j)+3] = c[2] + X*nx[2] + Y*ny[2];
            
            nc[4*(Np*i+j)+0] = 1.0;
            nc[4*(Np*i+j)+1] = nz[0];
            nc[4*(Np*i+j)+2] = nz[1];
            nc[4*(Np*i+j)+3] = nz[2];
        }
    }
}

void free_grid(struct Grid *g)
{
    if(g->xi != NULL)
        free(g->xi);
    if(g->xc != NULL)
        free(g->xc);
    if(g->nc != NULL)
        free(g->nc);
}
