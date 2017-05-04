#include <stdio.h>
#include "par.h"
#include "metric.h"

int setup_metric(struct parList *pars)
{
    int choice = pars->metric;
    int err = 0;

    if(choice == 0)
    {
        metric_ig  = &metric_ig_flat_cart;
        metric_dig = &metric_dig_flat_cart;
        metric_cart2coord = &metric_cart2coord_flat_cart;
        metric_vec2coordb = &metric_vec2coordb_flat_cart;
        metric_shadow = &metric_shadow_flat_cart;
        metric_fix_domain = &metric_fix_domain_flat_cart;
    }
    else if(choice == 1)
    {
        metric_ig  = &metric_ig_flat_sph;
        metric_dig = &metric_dig_flat_sph;
        metric_cart2coord = &metric_cart2coord_flat_sph;
        metric_vec2coordb = &metric_vec2coordb_flat_sph;
        metric_shadow = &metric_shadow_flat_sph;
        metric_fix_domain = &metric_fix_domain_flat_sph;
    }
    else if(choice == 2)
    {
        metric_ig  = &metric_ig_schw_sc;
        metric_dig = &metric_dig_schw_sc;
        metric_cart2coord = &metric_cart2coord_schw_sc;
        metric_vec2coordb = &metric_vec2coordb_schw_sc;
        metric_shadow = &metric_shadow_schw_sc;
        metric_fix_domain = &metric_fix_domain_schw_sc;
    }
    else if(choice == 3)
    {
        metric_ig  = &metric_ig_schw_ks;
        metric_dig = &metric_dig_schw_ks;
        metric_cart2coord = &metric_cart2coord_schw_ks;
        metric_vec2coordb = &metric_vec2coordb_schw_ks;
        metric_shadow = &metric_shadow_schw_ks;
        metric_fix_domain = &metric_fix_domain_schw_ks;
    }
    else
    {
        printf("Bad Metric choice: %d\n", choice);
        err = 1;
    }

    return err;
}
