#include <stdio.h>

void output_print_init(double *args, int n, char *filename)
{
    FILE *f = fopen(filename, "w");
    int i;
    
    for(i=0; i<n; i++)
        fprintf(f, "%.12lg ", args[i]);
    fprintf(f, "\n");

    fclose(f);
}

void output_print_step(double t, double *x, int n, char *filename)
{
    FILE *f = fopen(filename, "a");
    int i;
    fprintf(f, "%.12lg", t);
    for(i=0; i<n; i++)
        fprintf(f, " %.12lg", x[i]);
    fprintf(f, "\n");

    fclose(f);
}
