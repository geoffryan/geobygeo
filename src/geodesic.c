#include <stdio.h>
#include <stdlib.h>
#include "par.h"
#include "metric.h"
#include "ode.h"
#include "output.h"
#include "geodesic.h"

/*
 * Routines for numerically solving the Geodesic Equation.  
 * Coordinates are all expressed using contravariant values: x[mu] == x^\mu. 
 * Four-velocities are all contravariant, since these are the conjugate
 * momenta to x^mu in the Hamiltonian picture of geodesic motion.
 * u[mu] == u_\mu.
 *
 * The organization is based on the Hamiltonian picture, where the free motion
 * of a particle in curved spacetime is determined by the Hamiltonian:
 * H = 1/2 * p_\mu p^\mu.
 */

void geo_dHdP(double t, double *x, double *u, void *args, double *xdot)
{
    // The first part of Hamilton's equations: dq/dt = dH/dp.  Simplifies
    // to dx^\mu/dt = u^\mu.  Fills 'xdot' with the current values of dx/dt.

    double g[16];
    metric_ig(g, x, args);

    int i,j;

    for(i=0; i<4; i++)
    {
        xdot[i] = 0.0;
        for(j=0; j<4; j++)
            xdot[i] += g[4*i+j]*u[j];
    }
}

void geo_dHdQ(double t, double *x, double *u, void *args, double *udot)
{
    // The second part of Hamilton's equations: dp/dt = -dH/dq.  Simplifies
    // to d u_\mu/dt = -1/2 * u_\nu * u_\rho * d_\mu g^{\nu\rho}.  
    // Fills 'udot' with the current values of du/dt.

    double dg[64];
    metric_dig(dg, x, args);

    int i,j,k;

    for(i=0; i<4; i++)
    {
        udot[i] = 0.0;
        for(j=0; j<4; j++)
            for(k=0; k<4; k++)
                udot[i] += dg[16*i+4*j+k]*u[j]*u[k];
        udot[i] *= -0.5;
    }
}

void geo_xudot(double t, double *xu, void *args, double *xudot)
{
    // Fills 'xudot' with the current values of dx/dt and du/dt.

    double *x = &(xu[0]);
    double *u = &(xu[4]);
    double *xdot = &(xudot[0]);
    double *udot = &(xudot[4]);

    geo_dHdP(t, x, u, args, xdot);
    geo_dHdQ(t, x, u, args, udot);
}

void geo_integrate_generic(double *x0, double *u0, double *x, double *u,
                           double t0, double t1, double dt0, void *args,
                           double (*step)(double, double*, double*, int, void*,
                             void (*)(double,double*,void*,double*)),
                           char filename[])
{
    // Integrate a geodesic with initial data x^mu=x0, u_mu=u0 from t0 to t1.
    // Begin with timestep dt0, and use 'step' integrator.
    
    int i;
    double t, dt;
    double *xu = (double *) malloc(8 * sizeof(double));

    for(i=0; i<4; i++)
    {
        xu[i] = x0[i];
        xu[4+i] = u0[i];
    }

    t = t0;
    dt = dt0;
   
    //TODO: set n_args?
    if(filename != NULL)
    {
        output_print_init((double *)args, 0, filename);
        output_print_step(t, xu, 8, filename);
    }

    while(t < t1)
    {
        dt0 = step(t, xu, &dt, 8, args, &geo_xudot);
        t += dt0;

        if(filename != NULL)
            output_print_step(t, xu, 8, filename);
        for(i=0; i<8; i++)
            if(xu[i] != xu[i])
                t = t1;
    }

    for(i=0; i<4; i++)
    {
        x[i] = xu[i];
        u[i] = xu[4+i];
    }

    free(xu);
}

void geo_integrate_surface(double *x0, double *u0, double *x, double *u,
                           double t0, double t1, double dt0, int mu, double C,
                           void *args,
                           double (*step)(double, double*, double*, int, void*,
                             void (*)(double,double*,void*,double*)),
                           char filename[])
{
    // Integrate a geodesic with initial data x^mu=x0, u_mu=u0 until it hits
    // surface x[mu] = C. Ends also if t > t1.
    // Begin with timestep dt0, and use 'step' integrator.
    
    int i, j;
    double t, dt;
    double xu[8];

    for(i=0; i<4; i++)
    {
        xu[i] = x0[i];
        xu[4+i] = u0[i];
    }

    t = t0;
    dt = dt0;
   
    //TODO: set n_args?
    if(filename != NULL)
    {
        output_print_init((double *)args, 0, filename);
        output_print_step(t, xu, 8, filename);
    }

    double xu0[8];

    while(((x0[mu] < C && xu[mu] < C) || (x0[mu] > C && xu[mu] > C)) && t > t1
            && !metric_shadow(xu, args))
    {
        for(i=0; i<8; i++)
            xu0[i] = xu[i];

        dt0 = step(t, xu, &dt, 8, args, &geo_xudot);
        metric_fix_domain(&(xu[0]), &(xu[4]), args);
        t += dt0;

        if(filename != NULL)
            output_print_step(t, xu, 8, filename);
        for(i=0; i<8; i++)
            if(xu[i] != xu[i])
                t = t1;

        double g[16];
        metric_ig(g, xu, args);
        double sum = 0.0;
        //double sum2 = xu[4]*xu[4]*g[0];
        for(i=1; i<4; i++)
            for(j=1; j<4; j++)
                sum += g[4*i+j]*xu[4+i]*xu[4+j];
        for(i=1; i<4; i++)
            sum += 2*g[i]*xu[4]*xu[4+i];

    //    printf("%g %lg %lg %lg\n", t, xu[1], sum+sum2, sum/sum2);

//        printf("%g %g %g %g %g %g %g %g %g %g\n", t, dt0, xu[0], xu[1], xu[2], xu[3],
 //xu[4], xu[5], xu[6], xu[7]);
    }

    if(t > t1 && !metric_shadow(xu, args))
    {
        double w = (C-xu0[mu]) / (xu[mu]-xu0[mu]);

        for(i=0; i<4; i++)
        {
            x[i] = (1.0-w)*xu0[i] + w*xu[i];
            u[i] = (1.0-w)*xu0[4+i] + w*xu[4+i];
        }
    }
    else
    {
        for(i=0; i<4; i++)
        {
            x[i] = 0.0;
            u[i] = 0.0;
        }
    }
}
