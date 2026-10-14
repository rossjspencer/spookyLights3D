#ifndef SPOOKYLIGHTS3D_UTILS_H
#define SPOOKYLIGHTS3D_UTILS_H

#include "spookyLights3D.h"
extern double TOL;
extern double eye4x4[4][4];

// multiply transformation matrices, leaving the result in b
static inline void matMult(double A[4][4], double B[4][4])
{
 double C[4][4];
 for (int i = 0; i < 4; i++)
  for (int j = 0; j < 4; j++)
   C[i][j] = A[i][0]*B[0][j] + A[i][1]*B[1][j] + A[i][2]*B[2][j] + A[i][3]*B[3][j];
 memcpy(B, C, sizeof C);
}

static inline void matVecMult(double A[4][4], struct point3D *pt)
{
 struct point3D pr;
 pr.px = A[0][0]*pt->px + A[0][1]*pt->py + A[0][2]*pt->pz + A[0][3]*pt->pw;
 pr.py = A[1][0]*pt->px + A[1][1]*pt->py + A[1][2]*pt->pz + A[1][3]*pt->pw;
 pr.pz = A[2][0]*pt->px + A[2][1]*pt->py + A[2][2]*pt->pz + A[2][3]*pt->pw;
 pr.pw = A[3][0]*pt->px + A[3][1]*pt->py + A[3][2]*pt->pz + A[3][3]*pt->pw;
 *pt = pr;
}

static inline void normalize(struct point3D *v)
{
 double l = sqrt(v->px*v->px + v->py*v->py + v->pz*v->pz);
 v->px /= l;
 v->py /= l;
 v->pz /= l;
}

static inline double dot(struct point3D *u, struct point3D *v)
{
 return u->px*v->px + u->py*v->py + u->pz*v->pz;
}

static inline struct point3D *cross(struct point3D *u, struct point3D *v)
{
 struct point3D *cp = (struct point3D *)calloc(1, sizeof(struct point3D));
 if (cp == NULL) return NULL;
 cp->px = u->py*v->pz - v->py*u->pz;
 cp->py = u->pz*v->px - v->pz*u->px;
 cp->pz = u->px*v->py - v->px*u->py;
 cp->pw = 1;
 return cp;
}

static inline void initRay(struct ray3D *ray, struct point3D *p0, struct point3D *d)
{
 ray->p0 = *p0;
 ray->d = *d;
}

struct pointLS *newPLS(struct point3D *p0, double r, double g, double b);
void rayTransform(struct ray3D *ray_orig, struct ray3D *ray_transformed, struct object3D *obj);
void normalTransform(struct point3D *n_orig, struct point3D *n_transformed, struct object3D *obj);
void insertObject(struct object3D *o, struct object3D **list);
void insertPLS(struct pointLS *l, struct pointLS **list);

// each geometric primitive supplies its own intersection function
struct object3D *newPlane(double ra, double rd, double rs, double rg, double r, double g, double b, double shiny);
struct object3D *newSphere(double ra, double rd, double rs, double rg, double r, double g, double b, double shiny);
struct object3D *newCyl(double ra, double rd, double rs, double rg, double r, double g, double b, double shiny);
struct object3D *newCone(double ra, double rd, double rs, double rg, double r, double g, double b, double shiny);
struct object3D *newCube(double ra, double rd, double rs, double rg, double r, double g, double b, double shiny);
struct object3D *newIcosahedron(double ra, double rd, double rs, double rg, double r, double g, double b, double shiny);
void planeIntersect(struct object3D *plane, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n);
void sphereIntersect(struct object3D *sphere, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n);
void cylIntersect(struct object3D *cylinder, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n);
void coneIntersect(struct object3D *cone, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n);
void cubeIntersect(struct object3D *cube, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n);
void icosahedronIntersect(struct object3D *icosahedron, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n);

void invert(double *T, double *Tinv);
void RotateXMat(double T[4][4], double theta);
void RotateYMat(double T[4][4], double theta);
void RotateZMat(double T[4][4], double theta);
void TranslateMat(double T[4][4], double tx, double ty, double tz);
void ScaleMat(double T[4][4], double sx, double sy, double sz);
void RotateX(struct object3D *o, double theta);
void RotateY(struct object3D *o, double theta);
void RotateZ(struct object3D *o, double theta);
void Translate(struct object3D *o, double tx, double ty, double tz);
void Scale(struct object3D *o, double sx, double sy, double sz);
void printmatrix(double mat[4][4]);
struct view *setupView(struct point3D *e, struct point3D *g, struct point3D *up,
                       double f, double wl, double wt, double wsize);
struct image *newImage(int size_x, int size_y);
void imageOutput(struct image *im, const char *filename);
void deleteImage(struct image *im);
void cleanup(struct object3D *o_list, struct pointLS *l_list);

#endif
