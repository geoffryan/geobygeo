#ifndef GEO_GRID
#define GEO_GRID

struct Grid
{
    int mode;//Type of grid
    int N;   //Total number of grid points
    int N1;  //Number of grid points in X1 dir
    int N2;  //Number of grid points in X2 dir

    double distance;        //Distance to origin
    double inclination;     //Angle (in rad) from z-axis
    double azimuth;         //Azimuthal angle (in rad)

    double X1a;
    double X1b;
    double X2a;
    double X2b;

    double *xi;     //Image plane coordinates of each grid point
    double *xc;     //Cartesian bulk-space coordinates of each grid point
    double *nc;     //Cartesian normal vector to image plane at each grid
                    // point.
};

const static struct Grid GRID_DEFAULT = {
    .mode = 0,
    .N = 1,
    .N1 = 1,
    .N2 = 1,
    .distance = 1.0e3,
    .inclination = 0.0,
    .azimuth = 0.0,
    .X1a = 0.0,
    .X1b = 0.0,
    .X2a = 0.0,
    .X2b = 0.0,
    .xi = NULL,
    .xc = NULL,
    .nc = NULL
};

int setup_grid(struct Grid *g, struct parList *par);
void grid_initialize(struct Grid *g, struct parList *par);
void grid_generate_cart(struct Grid *g);
void grid_generate_elps(struct Grid *g);
void free_grid(struct Grid *g);

#endif
