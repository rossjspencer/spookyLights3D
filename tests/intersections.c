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

static int checkNormal(struct object3D *obj, struct point3D origin, struct point3D direction,
                       struct point3D expectedNormal)
{
 struct ray3D ray;
 struct point3D point, normal;
 double lambda = -1;
 initRay(&ray, &origin, &direction);
 obj->intersect(obj, &ray, &lambda, &point, &normal);
 return lambda > 0 && near(normal.px, expectedNormal.px) &&
        near(normal.py, expectedNormal.py) && near(normal.pz, expectedNormal.pz);
}

static int testNewPrimitives(void)
{
 struct object3D *cone = newCone(.1, .8, .2, 0, 1, 1, 1, 8);
 struct object3D *cube = newCube(.1, .8, .2, 0, 1, 1, 1, 8);
 struct object3D *icosahedron = newIcosahedron(.1, .8, .2, 0, 1, 1, 1, 8);
 if (!cone || !cube || !icosahedron) return 10;

 // test the cone's curved wall, circular base, and inward origin
 if (!checkHit(cone, (struct point3D){2,0,.5,1}, (struct point3D){-1,0,0,0},
               1.5, (struct point3D){.5,0,.5,1})) return 11;
 double d = sqrt(.5);
 if (!checkNormal(cone, (struct point3D){2,0,.5,1}, (struct point3D){-1,0,0,0},
                  (struct point3D){d,0,d,0})) return 12;
 if (!checkHit(cone, (struct point3D){0,0,-2,1}, (struct point3D){0,0,1,0},
               2, (struct point3D){0,0,0,1})) return 13;
 if (!checkNormal(cone, (struct point3D){0,0,-2,1}, (struct point3D){0,0,1,0},
                  (struct point3D){0,0,-1,0})) return 14;
 if (!checkHit(cone, (struct point3D){0,0,.5,1}, (struct point3D){1,0,0,0},
               .5, (struct point3D){.5,0,.5,1})) return 15;

 // test entry and exit faces of the canonical cube
 if (!checkHit(cube, (struct point3D){2,.5,.5,1}, (struct point3D){-1,0,0,0},
               1, (struct point3D){1,.5,.5,1})) return 16;
 if (!checkNormal(cube, (struct point3D){2,.5,.5,1}, (struct point3D){-1,0,0,0},
                  (struct point3D){1,0,0,0})) return 17;
 if (!checkHit(cube, (struct point3D){.5,.5,.5,1}, (struct point3D){0,0,-1,0},
               .5, (struct point3D){.5,.5,0,1})) return 18;
 if (!checkNormal(cube, (struct point3D){.5,.5,.5,1}, (struct point3D){0,0,-1,0},
                  (struct point3D){0,0,-1,0})) return 19;

 // verify that all twenty faces of the centered icosahedron can be intersected from inside
 double phi = (1 + sqrt(5.0)) / 2;
 double vertices[12][3] = {
  {0, .5, phi/2}, {0, -.5, phi/2}, {0, .5, -phi/2}, {0, -.5, -phi/2},
  {.5, phi/2, 0}, {-.5, phi/2, 0}, {.5, -phi/2, 0}, {-.5, -phi/2, 0},
  {phi/2, 0, .5}, {-phi/2, 0, .5}, {phi/2, 0, -.5}, {-phi/2, 0, -.5}
 };
 int faces[20][3] = {
  {0, 1, 8}, {0, 4, 5}, {0, 5, 9}, {0, 8, 4}, {0, 9, 1},
  {1, 6, 8}, {5, 4, 2}, {7, 1, 9}, {7, 3, 6}, {7, 6, 1},
  {10, 2, 4}, {10, 3, 2}, {10, 4, 8}, {10, 6, 3}, {10, 8, 6},
  {11, 2, 3}, {11, 3, 7}, {11, 5, 2}, {11, 7, 9}, {11, 9, 5}
 };
 for (int i = 0; i < 20; i++)
 {
  struct point3D center = {0,0,0,0};
  for (int j = 0; j < 3; j++)
  {
   center.px += vertices[faces[i][j]][0]/3;
   center.py += vertices[faces[i][j]][1]/3;
   center.pz += vertices[faces[i][j]][2]/3;
  }
  struct point3D origin = {0,0,0,1};
  struct point3D expected = center;
  expected.pw = 1;
  struct ray3D ray;
  struct point3D hit, normal;
  double lambda = -1;
  initRay(&ray, &origin, &center);
  icosahedron->intersect(icosahedron, &ray, &lambda, &hit, &normal);
  if (!near(lambda, 1) || !near(hit.px, expected.px) ||
      !near(hit.py, expected.py) || !near(hit.pz, expected.pz) ||
      dot(&normal, &center) <= 0) return 20;
 }

 // verify scaled and translated intersections for the three additional shapes
 Scale(cone, 2, 2, 2);
 Translate(cone, 3, 0, 0);
 invert(&cone->T[0][0], &cone->Tinv[0][0]);
 if (!checkHit(cone, (struct point3D){7,0,1,1}, (struct point3D){-1,0,0,0},
               3, (struct point3D){4,0,1,1})) return 21;
 Scale(cube, 2, 1, 1);
 Translate(cube, -2, 0, 0);
 invert(&cube->T[0][0], &cube->Tinv[0][0]);
 if (!checkHit(cube, (struct point3D){1,.5,.5,1}, (struct point3D){-1,0,0,0},
               1, (struct point3D){0,.5,.5,1})) return 22;
 Scale(icosahedron, 2, 1, .5);
 Translate(icosahedron, 2, 0, 0);
 invert(&icosahedron->T[0][0], &icosahedron->Tinv[0][0]);
 struct point3D center = {0,0,0,0};
 for (int j = 0; j < 3; j++)
 {
  center.px += vertices[faces[0][j]][0]/3;
  center.py += vertices[faces[0][j]][1]/3;
  center.pz += vertices[faces[0][j]][2]/3;
 }
 struct point3D direction = {2*center.px,center.py,.5*center.pz,0};
 if (!checkHit(icosahedron, (struct point3D){2,0,0,1}, direction,
               1, (struct point3D){2+direction.px,direction.py,direction.pz,1})) return 23;

 // confirm a ray pointing away from the solid is a miss
 struct ray3D missRay;
 struct point3D missOrigin = {3,3,3,1};
 struct point3D missDirection = {1,0,0,0};
 struct point3D point, normal;
 double lambda = -1;
 initRay(&missRay, &missOrigin, &missDirection);
 cube->intersect(cube, &missRay, &lambda, &point, &normal);
 if (lambda >= 0) return 24;
 icosahedron->intersect(icosahedron, &missRay, &lambda, &point, &normal);
 if (lambda >= 0) return 25;
 cone->intersect(cone, &missRay, &lambda, &point, &normal);
 if (lambda >= 0) return 26;

 free(cone);
 free(cube);
 free(icosahedron);
 puts("cone, cube, and icosahedron intersection checks passed");
 return 0;
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
 int newPrimitiveResult = testNewPrimitives();
 if (newPrimitiveResult != 0) return newPrimitiveResult;

 puts("plane, sphere, and cylinder intersection checks passed");
 return 0;
}
