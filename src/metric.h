#ifndef GEO_METRIC
#define GEO_METRIC

/*
 * metric_ig(g, x, args):   Fills the instantiated array g with the values of
 *                          g^{\mu\nu} at position x^\mu and parameters 'args'
 *                          It is assumed 'g' has at least length 16, and x
 *                          has length at least 4.
 *
 * metric_dig(dg, x, args): Fills the instantiated array 'dg' with the values
 *                          of \partial_\mu g^{\nu\rho} at position x^\mu and 
 *                          parameters 'args'.  It is assumed 'dg' has length
 *                          at least 64, and x has length at least 4.  'dg' is 
 *                          filled with the derivative as the outer-most index.
 *                          ie. d_i g^[jk] == dg[16*i+4*j+k]
 */

void  (*metric_ig)(double *, double *, void *);
void (*metric_dig)(double *, double *, void *);
void (*metric_cart2coord)(double *, double *, void *);
void (*metric_vec2coordb)(double *, double *, double *, void *);
int (*metric_shadow)(double *, void *);

int setup_metric(struct parList *pars);

void  metric_ig_flat_cart(double *g, double *x, void *args);
void metric_dig_flat_cart(double *dg, double *x, void *args);
void metric_cart2coord_flat_cart(double *xc, double *x, void *args);
void metric_vec2coordb_flat_cart(double *x, double *uc, double *u, void *args);
int metric_shadow_flat_cart(double *x, void *args);

void  metric_ig_flat_sph(double *g, double *x, void *args);
void metric_dig_flat_sph(double *dg, double *x, void *args);
void metric_cart2coord_flat_sph(double *xc, double *x, void *args);
void metric_vec2coordb_flat_sph(double *x, double *uc, double *u, void *args);
int metric_shadow_flat_sph(double *x, void *args);

void  metric_ig_schw_sc(double *g, double *x, void *args);
void metric_dig_schw_sc(double *dg, double *x, void *args);
void metric_cart2coord_schw_sc(double *xc, double *x, void *args);
void metric_vec2coordb_schw_sc(double *x, double *uc, double *u, void *args);
int metric_shadow_schw_sc(double *x, void *args);

void  metric_ig_schw_ks(double *g, double *x, void *args);
void metric_dig_schw_ks(double *dg, double *x, void *args);
void metric_cart2coord_schw_ks(double *xc, double *x, void *args);
void metric_vec2coordb_schw_ks(double *x, double *uc, double *u, void *args);
int metric_shadow_schw_ks(double *x, void *args);

#endif
