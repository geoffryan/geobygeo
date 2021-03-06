#ifndef GEO_GEO
#define GEO_GEO

void geo_dHdP(double t, double *x, double *u, void *args, double *xdot);
void geo_dHdQ(double t, double *x, double *u, void *args, double *udot);
void geo_xudot(double t, double *xu, void *args, double *xudot);

void geo_integrate_generic(double *x0, double *u0, double *x, double *u,
                           double t0, double t1, double dt0, void *args,
                           double (*step)(double, double*, double*, int, void*,
                                void (*)(double,double*,void*,double*)),
                           char filename[]);
void geo_integrate_surface(double *x0, double *u0, double *x, double *u,
                           double t0, double t1, double dt0, int mu, double C,
                           void *args,
                           double (*step)(double, double*, double*, int, void*,
                             void (*)(double,double*,void*,double*)),
                           char filename[]);

#endif
