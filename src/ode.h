#ifndef GEO_ODE
#define GEO_ODE

void forward_euler(double t, double *x, double *dx, int n, void *args,
                    void (*xdot)(double,double*,void*,double*));
void rk2(double t, double *x, double *dx, int n, void *args,
                    void (*xdot)(double,double*,void*,double*));
void rk4(double t, double *x, double *dx, int n, void *args,
                    void (*xdot)(double,double*,void*,double*));

#endif
