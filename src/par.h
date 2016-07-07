#ifndef GEO_PAR
#define GEO_PAR

enum{VAR_DBL, VAR_INT, VAR_LON, VAR_STR};

struct parList
{
    int metric;
    int grid;
    int N1;
    int N2;
    double X1a;
    double X1b;
    double X2a;
    double X2b;
    double distance;
    double inclination;
    double azimuth;
};

const static struct parList PAR_DEFAULT = {
    .metric = 0,
    .grid = 0,
    .N1 = 1,
    .N2 = 1,
    .X1a = 0.0,
    .X1b = 0.0,
    .X2a = 0.0,
    .X2b = 0.0,
    .inclination = 0.0,
    .azimuth = 0.0
};

int readvar(char filename[], char key[], int vtype, void *ptr);
void read_pars(struct parList *theParList, char filename[]);
void print_pars(struct parList *theParList, char filename[]);

#endif
