#define gv_message(msg) fprintf(stdout, "[geovis] %s: %s\n", __FUNCTION__, msg);
#define gv_error(msg) fprintf(stderr, "[ERROR] %s: %s\n", __FUNCTION__, msg);

//Callback Functions
void IdleFunc();
void ReshapeFunc(int width, int height);
void DisplayFunc();
void KeyboardFunc(unsigned char key, int x, int y);
void SpecialFunc(int key, int x, int y);

//Everything Else
int gv_create_window(char *title, int x, int y, int width, int height);
void gv_init(int argc, char *argv[]);
void gv_init_gl(int width, int height);
void gv_exit();
void gv_build_scene();
void gv_build_cutplanes();
void gv_build_cutplane(int ax);
void gv_draw_cutplanes();
void gv_build_geodesics();
void gv_build_geodesic(int p);
void gv_draw_geodesics();
void gv_draw_geodesic(int p);
void gv_wirecube(float side);
void gv_cmap(double val, int cmap, GLfloat *rrr, GLfloat *ggg, GLfloat *bbb);


struct gv_cutplane
{
    GLuint tex;
    int ind;
    int cmap;
    GLfloat X0[3];
    GLfloat X1[3];
    GLfloat X2[3];
    GLfloat X3[3];
};
typedef struct gv_cutplane gv_cutplane;

struct gv_geodesic
{
    int n;
    GLfloat *lambda;
    GLfloat *t;
    GLfloat *x;
    GLfloat *y;
    GLfloat *z;
};
typedef struct gv_geodesic gv_geodesic;

const struct gv_cutplane cutplane_default = {0, 0, 0, {0.0f,0.0f,0.0f}, {0.0f,0.0f,0.0f},
    {0.0f,0.0f,0.0f}, {0.0f,0.0f,0.0f}};
const struct gv_geodesic geodesic_default = {0, NULL, NULL, NULL, NULL, NULL};

