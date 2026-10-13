#ifndef SPOOKYLIGHTS3D_H
#define SPOOKYLIGHTS3D_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265354

struct image {
 void *rgbdata;
 int sx;
 int sy;
};

struct point3D {
 double px, py, pz, pw;
};

struct ray3D {
 struct point3D p0;
 struct point3D d;
};

struct albedosPhong {
 double ra; // ambient
 double rd; // diffuse
 double rs; // specular
 double rg; // recursive reflection
};

struct colourRGB {
 double R, G, B;
};

struct object3D {
 struct albedosPhong alb;
 struct colourRGB col;
 double T[4][4];
 double Tinv[4][4];
 double shinyness;
 // each primitive supplies its own ray intersection function
 void (*intersect)(struct object3D *obj, struct ray3D *ray, double *lambda,
                   struct point3D *p, struct point3D *n);
 struct object3D *next;
};

struct pointLS {
 struct colourRGB col;
 struct point3D p0;
 struct pointLS *next;
};

struct view {
 struct point3D e, u, v, w;
 double f, wl, wt, wsize;
 double W2C[4][4];
 double C2W[4][4];
};

extern struct object3D *object_list;
extern struct pointLS *light_list;
extern int MAX_DEPTH;

void buildScene(void);
void rayTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os);
void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os,
                  struct object3D **obj, struct point3D *p, struct point3D *n);
void rtShade(struct object3D *obj, struct point3D *p, struct point3D *n,
             struct ray3D *ray, int depth, struct colourRGB *col);

#endif
