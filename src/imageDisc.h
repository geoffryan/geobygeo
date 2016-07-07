#ifndef GEO_IMGDISC
#define GEO_IMGDISC

void imageDisc(double center[], double n[], double width, double height, 
                int nx, int ny, void *args, char filename[]);
void imageDiscGrid(struct Grid *g, void *args, char filename[]);

#endif
