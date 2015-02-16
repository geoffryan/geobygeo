#ifndef GEO_TEST
#define GEO_TEST

void oscillator(double t, double *x, void *args, double *xdot);
void test_oscillator(double (*step)(double, double*, double *, int, void *,
                        void (*)(double,double*,void*,double*)), double t0,
                        double t1, int n, char *filename);

#endif
