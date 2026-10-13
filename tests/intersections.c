#include "utils.h"

static int near(double a, double b)
{
 return fabs(a - b) < 1e-6;
}

static int checkHit(struct object3D *obj, struct point3D origin, struct point3D direction,
                    double expectedLambda, struct point3D expectedPoint)
{
 struct ray3D ray;
 struct point3D point, normal;
 double lambda = -1;
 initRay(&ray, &origin, &direction);
 obj->intersect(obj, &ray, &lambda, &point, &normal);
 return near(lambda, expectedLambda) && near(point.px, expectedPoint.px) &&
        near(point.py, expectedPoint.py) && near(point.pz, expectedPoint.pz);
}

int main(void)
{
 struct object3D *sphere = newSphere(.1, .8, .2, 0, 1, 1, 1, 8);
 struct object3D *plane = newPlane(.1, .8, .2, 0, 1, 1, 1, 8);
 struct object3D *cylinder = newCyl(.1, .8, .2, 0, 1, 1, 1, 8);
 if (sphere == NULL || plane == NULL || cylinder == NULL) return 1;

 // verify hits against the canonical surfaces
 if (!checkHit(sphere, (struct point3D){0,0,-3,1}, (struct point3D){0,0,1,0},
               2, (struct point3D){0,0,-1,1})) return 2;
 if (!checkHit(plane, (struct point3D){0,0,-3,1}, (struct point3D){0,0,1,0},
               3, (struct point3D){0,0,0,1})) return 3;
 if (!checkHit(cylinder, (struct point3D){2,0,.5,1}, (struct point3D){-1,0,0,0},
               1, (struct point3D){1,0,.5,1})) return 4;
 if (!checkHit(cylinder, (struct point3D){0,0,2,1}, (struct point3D){0,0,-1,0},
               1, (struct point3D){0,0,1,1})) return 5;

 // verify an intersection still works after non-uniform scaling
 Scale(sphere, 2, 1, 1);
 invert(&sphere->T[0][0], &sphere->Tinv[0][0]);
 if (!checkHit(sphere, (struct point3D){3,0,0,1}, (struct point3D){-1,0,0,0},
               1, (struct point3D){2,0,0,1})) return 6;

 // verify a miss reports a negative distance
 struct ray3D ray;
 struct point3D point, normal;
 double lambda = -1;
 struct point3D origin = {0,0,-3,1};
 struct point3D direction = {0,1,0,0};
 initRay(&ray, &origin, &direction);
 sphere->intersect(sphere, &ray, &lambda, &point, &normal);
 if (lambda >= 0) return 7;

 free(sphere);
 free(plane);
 free(cylinder);
 puts("plane, sphere, and cylinder intersection checks passed");
 return 0;
}
