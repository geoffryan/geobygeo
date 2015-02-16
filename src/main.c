#include <stdio.h>
#include <math.h>
#include "test.h"
#include "ode.h"
#include "metric.h"
#include "geodesic.h"

int main(int argc, char *argv[])
{
    int n = 200;

    metric_ig  = &metric_ig_flat_sph;
    metric_dig = &metric_dig_flat_sph;

    double x0[4];
    double u0[4];

    double r = 0.2;
    double b = 0.1;
    double v = 0.01;

    x0[0] = 0.0;
    x0[1] = 1.0;
    x0[2] = asin(b/r);
    x0[3] = 0.0;

    u0[0] = -1.0/sqrt(1.0-v*v);
    u0[1] = 0.0;
    u0[2] = 0.0;
    u0[3] = b*v;

    geo_integrate_generic(x0, u0, 0.0, 4.0, 0.01, NULL, &forward_euler, 
                            "fe.txt");
    geo_integrate_generic(x0, u0, 0.0, 4.0, 0.01, NULL, &rk2, 
                            "rk2.txt");
    geo_integrate_generic(x0, u0, 0.0, 4.0, 0.01, NULL, &rk4, 
                            "rk4.txt");

    //test_oscillator(&forward_euler, 0.0, 4.25*M_PI, n);
    //test_oscillator(&rk2, 0.0, 4.25*M_PI, n);
    //test_oscillator(&rk4, 0.0, 4.25*M_PI, n);

    return 0;
}
