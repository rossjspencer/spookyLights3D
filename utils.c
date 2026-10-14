// spookyLights3D geometry and rendering utilities
#include "utils.h"
// matrix inversion routines are provided by svdDynamic.c
int SVD(const double *A, int m, int n, double **U, double **w, double **V, double **rv1);
void InvertMatrix(const double *U, const double *w, const double *V, int n, double *I);

double TOL = 1e-6;
double eye4x4[4][4] = {
 {1.0, 0.0, 0.0, 0.0},
 {0.0, 1.0, 0.0, 0.0},
 {0.0, 0.0, 1.0, 0.0},
 {0.0, 0.0, 0.0, 1.0}
};

struct pointLS *newPLS(struct point3D *p0, double r, double g, double b)
{

 struct pointLS *ls=(struct pointLS *)calloc(1,sizeof(struct pointLS));
 if (!ls) fprintf(stderr,"Out of memory allocating light source!\n");
 else
 {
  memcpy(&ls->p0,p0,sizeof(struct point3D));

  ls->col.R=r;
  ls->col.G=g;
  ls->col.B=b;
 }
 return(ls);
}

void rayTransform(struct ray3D *ray_orig, struct ray3D *ray_transformed, struct object3D *obj)
{

 ray_transformed->d = ray_orig->d;
 ray_transformed->p0 = ray_orig->p0;
 matVecMult(obj->Tinv, &ray_transformed->p0);
 matVecMult(obj->Tinv, &ray_transformed->d);
}

void normalTransform(struct point3D *n_orig, struct point3D *n_transformed, struct object3D *obj)
{
 // transform the normal using the inverse transpose of the object transformation

 n_transformed->px = n_orig->px;
 n_transformed->py = n_orig->py;
 n_transformed->pz = n_orig->pz;
 n_transformed->pw = 0;

 double inverseTranspose[4][4];

 for (int i = 0; i < 4; i++)
 {
  for (int j = 0; j < 4; j++)
  {
   inverseTranspose[i][j] = obj->Tinv[j][i];
  }
 }

 matVecMult(inverseTranspose, n_transformed);
 normalize(n_transformed);
}

void insertObject(struct object3D *o, struct object3D **list)
{
 if (o==NULL) return;
 // insert this object into the scene
 if (*(list)==NULL)
 {
  *(list)=o;
  (*(list))->next=NULL;
 }
 else
 {
  o->next=(*(list))->next;
  (*(list))->next=o;
 }
}

struct object3D *newPlane(double ra, double rd, double rs, double rg, double r, double g, double b, double shiny)
{
 // create a plane with the supplied phong material and an identity transform

 struct object3D *plane=(struct object3D *)calloc(1,sizeof(struct object3D));

 if (!plane) fprintf(stderr,"Unable to allocate new plane, out of memory!\n");
 else
 {
  plane->alb.ra=ra;
  plane->alb.rd=rd;
  plane->alb.rs=rs;
  plane->alb.rg=rg;
  plane->col.R=r;
  plane->col.G=g;
  plane->col.B=b;
  plane->shinyness=shiny;
  plane->intersect=&planeIntersect;
  memcpy(&plane->T[0][0],&eye4x4[0][0],16*sizeof(double));
  memcpy(&plane->Tinv[0][0],&eye4x4[0][0],16*sizeof(double));
  plane->next=NULL;
}
 return(plane);
}

struct object3D *newSphere(double ra, double rd, double rs, double rg, double r, double g, double b, double shiny)
{
 // create a sphere with the supplied phong material and an identity transform

 struct object3D *sphere=(struct object3D *)calloc(1,sizeof(struct object3D));

 if (!sphere) fprintf(stderr,"Unable to allocate new sphere, out of memory!\n");
 else
 {
  sphere->alb.ra=ra;
  sphere->alb.rd=rd;
  sphere->alb.rs=rs;
  sphere->alb.rg=rg;
  sphere->col.R=r;
  sphere->col.G=g;
  sphere->col.B=b;
  sphere->shinyness=shiny;
  sphere->intersect=&sphereIntersect;
  memcpy(&sphere->T[0][0],&eye4x4[0][0],16*sizeof(double));
  memcpy(&sphere->Tinv[0][0],&eye4x4[0][0],16*sizeof(double));
  sphere->next=NULL;
 }
 return sphere;
}

struct object3D *newCyl(double ra, double rd, double rs, double rg, double r, double g, double b, double shiny)
{
 // we consider the canonical cylinder to be a circle of radius 1 on the x, y plane, extended from z = 0 to z = 1
 struct object3D *cylinder=(struct object3D *)calloc(1,sizeof(struct object3D));

 if (!cylinder) fprintf(stderr,"Unable to allocate new cylinder, out of memory!\n");
 else
 {
  cylinder->alb.ra=ra;
  cylinder->alb.rd=rd;
  cylinder->alb.rs=rs;
  cylinder->alb.rg=rg;
  cylinder->col.R=r;
  cylinder->col.G=g;
  cylinder->col.B=b;
  cylinder->shinyness=shiny;
  cylinder->intersect=&cylIntersect;
  memcpy(&cylinder->T[0][0],&eye4x4[0][0],16*sizeof(double));
  memcpy(&cylinder->Tinv[0][0],&eye4x4[0][0],16*sizeof(double));
  cylinder->next=NULL;
 }
 return cylinder;
}

struct object3D *newCone(double ra, double rd, double rs, double rg, double r, double g, double b, double shiny)
{
 // the canonical cone has a circular base of radius 1 at z = 0 and a tip at z = 1
 struct object3D *cone=(struct object3D *)calloc(1,sizeof(struct object3D));

 if (!cone) fprintf(stderr,"Unable to allocate new cone, out of memory!\n");
 else
 {
  cone->alb.ra=ra;
  cone->alb.rd=rd;
  cone->alb.rs=rs;
  cone->alb.rg=rg;
  cone->col.R=r;
  cone->col.G=g;
  cone->col.B=b;
  cone->shinyness=shiny;
  cone->intersect=&coneIntersect;
  memcpy(&cone->T[0][0],&eye4x4[0][0],16*sizeof(double));
  memcpy(&cone->Tinv[0][0],&eye4x4[0][0],16*sizeof(double));
  cone->next=NULL;
 }
 return cone;
}

struct object3D *newCube(double ra, double rd, double rs, double rg, double r, double g, double b, double shiny)
{
 // the canonical cube occupies the region x, y, z in [0, 1]
 struct object3D *cube=(struct object3D *)calloc(1,sizeof(struct object3D));

 if (!cube) fprintf(stderr,"Unable to allocate new cube, out of memory!\n");
 else
 {
  cube->alb.ra=ra;
  cube->alb.rd=rd;
  cube->alb.rs=rs;
  cube->alb.rg=rg;
  cube->col.R=r;
  cube->col.G=g;
  cube->col.B=b;
  cube->shinyness=shiny;
  cube->intersect=&cubeIntersect;
  memcpy(&cube->T[0][0],&eye4x4[0][0],16*sizeof(double));
  memcpy(&cube->Tinv[0][0],&eye4x4[0][0],16*sizeof(double));
  cube->next=NULL;
 }
 return cube;
}

struct object3D *newIcosahedron(double ra, double rd, double rs, double rg, double r, double g, double b, double shiny)
{
 // the canonical icosahedron is centered at the origin with 20 equilateral triangular faces
 // each edge has length 1 and the vertices follow the golden ratio coordinates
 struct object3D *icosahedron=(struct object3D *)calloc(1,sizeof(struct object3D));

 if (!icosahedron) fprintf(stderr,"Unable to allocate new icosahedron, out of memory!\n");
 else
 {
  icosahedron->alb.ra=ra;
  icosahedron->alb.rd=rd;
  icosahedron->alb.rs=rs;
  icosahedron->alb.rg=rg;
  icosahedron->col.R=r;
  icosahedron->col.G=g;
  icosahedron->col.B=b;
  icosahedron->shinyness=shiny;
  icosahedron->intersect=&icosahedronIntersect;
  memcpy(&icosahedron->T[0][0],&eye4x4[0][0],16*sizeof(double));
  memcpy(&icosahedron->Tinv[0][0],&eye4x4[0][0],16*sizeof(double));
  icosahedron->next=NULL;
 }
 return icosahedron;
}

void planeIntersect(struct object3D *plane, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n)
{

 // we're working with the canonical plane, so a point on the plane that works is (0, 0, 0).
 // this simplifies p1 - a to just -a.
 // now, we just have to solve (p + lambda*d) dot n = 0 for lambda, the formula for which is generously provided in the notes.

 // set the plane coordinate normal to the normal for this canonical plane, (0, 0, 1)
 struct point3D planeCoordinateNormal;
 planeCoordinateNormal.px = 0;
 planeCoordinateNormal.py = 0;
 planeCoordinateNormal.pz = 1;
 planeCoordinateNormal.pw = 0;

 // put the ray into canonical plane coordinates
 struct ray3D transformedRay;
 rayTransform(ray, &transformedRay, plane);

 // if this is 0, the line does not intersect the plane (or intersects it infinitely, neither is useful)
 double dDotN = dot(&transformedRay.d, &planeCoordinateNormal);
 if (fabs(dDotN) < TOL)
 {
  *lambda = -1;
  return;
 }

 // use the equation for lambda from the notes
 struct point3D negativeA;
 negativeA.px = -transformedRay.p0.px;
 negativeA.py = -transformedRay.p0.py;
 negativeA.pz = -transformedRay.p0.pz;
 negativeA.pw = 1;

 double lambdaCandidate = dot(&negativeA,&planeCoordinateNormal) / dDotN;

 // if the intersection is behind the ray origin
 if (lambdaCandidate < TOL)
 {
  *lambda = -1;
  return;
 }

 // calculating z is unnecessary because we know it will be 0
 double xCandidate = transformedRay.p0.px + lambdaCandidate * transformedRay.d.px;
 double yCandidate = transformedRay.p0.py + lambdaCandidate * transformedRay.d.py;

 // if the point doesn't actually lie on the plane, just the extended version of it
 if (xCandidate > 1 || yCandidate > 1 || xCandidate < -1 || yCandidate < -1)
 {
  *lambda = -1;
  return;
 }

 *lambda = lambdaCandidate;

 p->px = xCandidate;
 p->py = yCandidate;
 p->pz = 0;
 p->pw = 1;

 // put the intersection point back into world coordinates
 matVecMult(plane->T, p);

 // put the normal into world coordinates
 normalTransform(&planeCoordinateNormal, n, plane);

 // ensure the normal faces against the incoming ray
 if (dot(n, &ray->d) > 0)
 {
  n->px = -n->px;
  n->py = -n->py;
  n->pz = -n->pz;
 }
}

void sphereIntersect(struct object3D *sphere, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n)
{

 // we're working with the canonical sphere, which is a unit sphere centered at the origin.
 // substitute the transformed ray p + lambda*d into x^2 + y^2 + z^2 = 1.
 // this gives us a quadratic in lambda, which we can solve with the quadratic formula.

 // put the ray into canonical sphere coordinates
 struct ray3D transformedRay;
 rayTransform(ray, &transformedRay, sphere);

 // calculate the coefficients of the quadratic
 double quadraticA = dot(&transformedRay.d, &transformedRay.d);
 double quadraticB = 2 * dot(&transformedRay.p0, &transformedRay.d);
 double quadraticC = dot(&transformedRay.p0, &transformedRay.p0) - 1;

 double discriminant = quadraticB * quadraticB - 4 * quadraticA * quadraticC;

 // if the discriminant is negative, the ray does not intersect the sphere
 if (discriminant < 0)
 {
  *lambda = -1;
  return;
 }

 // calculate both possible intersection values
 double lambda1 = (-quadraticB - sqrt(discriminant)) / (2 * quadraticA);
 double lambda2 = (-quadraticB + sqrt(discriminant)) / (2 * quadraticA);

 // choose the first valid intersection in front of the ray origin
 double lambdaCandidate = -1;

 if (lambda1 >= TOL)
 {
  lambdaCandidate = lambda1;
 }

 if (lambda2 >= TOL && (lambdaCandidate < 0 || lambda2 < lambdaCandidate))
 {
  lambdaCandidate = lambda2;
 }

 // if neither intersection is in front of the ray origin
 if (lambdaCandidate < 0)
 {
  *lambda = -1;
  return;
 }

 *lambda = lambdaCandidate;

 // calculate the intersection point in canonical sphere coordinates
 p->px = transformedRay.p0.px + lambdaCandidate * transformedRay.d.px;
 p->py = transformedRay.p0.py + lambdaCandidate * transformedRay.d.py;
 p->pz = transformedRay.p0.pz + lambdaCandidate * transformedRay.d.pz;
 p->pw = 1;

 // for a unit sphere centered at the origin, the normal is just the
 // direction from the origin to the intersection point
 struct point3D sphereCoordinateNormal;
 sphereCoordinateNormal.px = p->px;
 sphereCoordinateNormal.py = p->py;
 sphereCoordinateNormal.pz = p->pz;
 sphereCoordinateNormal.pw = 0;

 // put the intersection point back into world coordinates
 matVecMult(sphere->T, p);

 // put the normal into world coordinates
 normalTransform(&sphereCoordinateNormal, n, sphere);
}

void cylIntersect(struct object3D *cylinder, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n)
{

 // we're working with the canonical cylinder, having a radius of 1, pointing along z with a z range of [0, 1].
 // to find intersection, we first check intersections with the outer wall described by x^2 + y^2 = 1
 // if there is an intersection, we confirm if it lies within the allowable z range.
 // then, we check if the ray intersects the top or bottom by checking if they intersect the planes z = 0 and z = 1
 // if there is an intersection, check if they lie within the circle x^2 + y^2 = 1.

 // put the ray into canonical cylinder coordinates
 struct ray3D transformedRay;
 rayTransform(ray, &transformedRay, cylinder);

 // start by considering the curved face, so we ignore the contribution of z
 struct point3D xyPointComponents;
 xyPointComponents.px = transformedRay.p0.px;
 xyPointComponents.py = transformedRay.p0.py;
 xyPointComponents.pz = 0;
 xyPointComponents.pw = 0;

 struct point3D xyDirectionComponents;
 xyDirectionComponents.px = transformedRay.d.px;
 xyDirectionComponents.py = transformedRay.d.py;
 xyDirectionComponents.pz = 0;
 xyDirectionComponents.pw = 0;

 // calculate the coefficients of the quadratic
 double quadraticA = dot(&xyDirectionComponents, &xyDirectionComponents);
 double quadraticB = 2 * dot(&xyPointComponents, &xyDirectionComponents);
 double quadraticC = dot(&xyPointComponents, &xyPointComponents) - 1;

 double discriminant = quadraticB * quadraticB - 4 * quadraticA * quadraticC;

 // choose the first valid intersection in front of the ray origin
 double lambdaCandidate = -1;

 double lambda1;
 double lambda2;

 // if the discriminant is negative, the ray does not intersect the cylinder wall
 if (discriminant >= 0 && fabs(quadraticA) >= TOL)
 {
  // calculate both possible intersection values
  lambda1 = (-quadraticB - sqrt(discriminant)) / (2 * quadraticA);
  lambda2 = (-quadraticB + sqrt(discriminant)) / (2 * quadraticA);

  if (lambda1 >= TOL)
  {
   // ensure the lambda1 point actually lies within the z range before updating lambdaCandidate
   double z = transformedRay.p0.pz + lambda1 * transformedRay.d.pz;
   if (z <= 1 + TOL && z >= -TOL)
   {
    lambdaCandidate = lambda1;
   }
  }

  if (lambda2 >= TOL && (lambdaCandidate < 0 || lambda2 < lambdaCandidate))
  {
   // ensure the lambda2 point actually lies within the z range before updating lambdaCandidate
   double z = transformedRay.p0.pz + lambda2 * transformedRay.d.pz;
   if (z <= 1 + TOL && z >= -TOL)
   {
    lambdaCandidate = lambda2;
   }
  }
 }

 // now, we consider the circular faces, or intersections with the planes defined by z = 1 and z = 0.
 // we get to ignore x and y now, we just solve for pz + lambda * dz = 1 and the same for 0!
 // only check the circular caps when the ray crosses their planes

 if (fabs(transformedRay.d.pz) > TOL)
 {
  // intersection with the plane z = 0
  lambda1 = -transformedRay.p0.pz / transformedRay.d.pz;

  // intersection with the plane z = 1
  lambda2 = (1 - transformedRay.p0.pz) / transformedRay.d.pz;

  double x;
  double y;

  if (lambda1 >= TOL && (lambdaCandidate < 0 || lambda1 < lambdaCandidate))
  {
   // ensure the lambda1 point actually lies within the circle before updating lambdaCandidate
   x = transformedRay.p0.px + lambda1 * transformedRay.d.px;
   y = transformedRay.p0.py + lambda1 * transformedRay.d.py;
   if (x*x + y*y <= 1 + TOL)
   {
     lambdaCandidate = lambda1;
   }
  }

  if (lambda2 >= TOL && (lambdaCandidate < 0 || lambda2 < lambdaCandidate))
  {
    // ensure the lambda1 point actually lies within the circle before updating lambdaCandidate
    x = transformedRay.p0.px + lambda2 * transformedRay.d.px;
    y = transformedRay.p0.py + lambda2 * transformedRay.d.py;
    if (x*x + y*y <= 1 + TOL)
    {
      lambdaCandidate = lambda2;
    }
  }
 }

 // if we got no valid intersections
 if (lambdaCandidate < 0)
 {
  *lambda = -1;
  return;
 }

 *lambda = lambdaCandidate;

 // calculate the intersection point in canonical cylinder coordinates
 p->px = transformedRay.p0.px + lambdaCandidate * transformedRay.d.px;
 p->py = transformedRay.p0.py + lambdaCandidate * transformedRay.d.py;
 p->pz = transformedRay.p0.pz + lambdaCandidate * transformedRay.d.pz;
 p->pw = 1;

 // for the cylinder, the normal will depend on whether we have hit the curved face or one of the curved faces.
 struct point3D cylinderCoordinateNormal;

 // top edge
 if (fabs(p->pz - 1) < TOL)
 {
  cylinderCoordinateNormal.px = 0;
  cylinderCoordinateNormal.py = 0;
  cylinderCoordinateNormal.pz = 1;
  cylinderCoordinateNormal.pw = 0;
 }
 // bottom edge
 else if (fabs(p->pz) < TOL)
 {
  cylinderCoordinateNormal.px = 0;
  cylinderCoordinateNormal.py = 0;
  cylinderCoordinateNormal.pz = -1;
  cylinderCoordinateNormal.pw = 0;
 }
 // curved face
 else
 {
  // in this case, the normal is the direction of the vector from the origin
  // to the point on the unit circle made by (x, y, 0). no need for gradients or anything
  cylinderCoordinateNormal.px = p->px;
  cylinderCoordinateNormal.py = p->py;
  cylinderCoordinateNormal.pz = 0;
  cylinderCoordinateNormal.pw = 0;
 }

 // put the intersection point back into world coordinates
 matVecMult(cylinder->T, p);

 // put the normal into world coordinates
 normalTransform(&cylinderCoordinateNormal, n, cylinder);
}

void coneIntersect(struct object3D *cone, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n)
{
 // the canonical cone has x*x + y*y = (1-z)*(1-z) for z in [0, 1], with a circular base at z = 0
 // substitute the ray into this equation to find possible intersections with the curved surface
 // put the ray into canonical cone coordinates
 struct ray3D transformedRay;
 rayTransform(ray, &transformedRay, cone);

 // use shorter names for the ray origin and direction components
 double x = transformedRay.p0.px;
 double y = transformedRay.p0.py;
 double z = transformedRay.p0.pz;
 double dx = transformedRay.d.px;
 double dy = transformedRay.d.py;
 double dz = transformedRay.d.pz;

 // collect the coefficients of the quadratic in lambda
 double quadraticA = dx*dx + dy*dy - dz*dz;
 double quadraticB = 2 * (x*dx + y*dy + (1-z)*dz);
 double quadraticC = x*x + y*y - (1-z)*(1-z);
 // keep track of the closest valid intersection and whether it was on the base
 double lambdaCandidate = -1;
 int hitBase = 0;

 // check the curved face, including the case where the quadratic becomes linear
 // a valid curved-face hit must have positive lambda and a z coordinate between 0 and 1
 if (fabs(quadraticA) < TOL)
 {
  // when the quadratic term is zero, solve the remaining linear equation
  if (fabs(quadraticB) >= TOL)
  {
   double candidate = -quadraticC / quadraticB;
   double hitZ = z + candidate*dz;
   if (candidate >= TOL && hitZ >= -TOL && hitZ <= 1 + TOL)
    lambdaCandidate = candidate;
  }
 }
 else
 {
  // otherwise, solve the quadratic if there are real roots
  double discriminant = quadraticB*quadraticB - 4*quadraticA*quadraticC;
  if (discriminant >= 0)
  {
   double lambda1 = (-quadraticB - sqrt(discriminant)) / (2*quadraticA);
   double lambda2 = (-quadraticB + sqrt(discriminant)) / (2*quadraticA);

   // check both roots and keep whichever valid intersection is closer
   double hitZ = z + lambda1*dz;
   if (lambda1 >= TOL && hitZ >= -TOL && hitZ <= 1 + TOL)
    lambdaCandidate = lambda1;

   hitZ = z + lambda2*dz;
   if (lambda2 >= TOL && hitZ >= -TOL && hitZ <= 1 + TOL &&
       (lambdaCandidate < 0 || lambda2 < lambdaCandidate))
    lambdaCandidate = lambda2;
  }
 }

 // check the circular base at z = 0
 // the ray must reach that plane inside the unit circle and before any previous hit
 if (fabs(dz) >= TOL)
 {
  double baseLambda = -z / dz;
  double baseX = x + baseLambda*dx;
  double baseY = y + baseLambda*dy;
  if (baseLambda >= TOL && baseX*baseX + baseY*baseY <= 1 + TOL &&
      (lambdaCandidate < 0 || baseLambda < lambdaCandidate))
  {
   lambdaCandidate = baseLambda;
   hitBase = 1;
  }
 }

 // return -1 if neither the curved face nor the base was intersected
 if (lambdaCandidate < 0)
 {
  *lambda = -1;
  return;
 }

 // calculate the intersection point in canonical cone coordinates
 *lambda = lambdaCandidate;
 p->px = x + lambdaCandidate*dx;
 p->py = y + lambdaCandidate*dy;
 p->pz = z + lambdaCandidate*dz;
 p->pw = 1;

 // the base normal points downward; the curved face normal follows the surface gradient
 struct point3D coneCoordinateNormal;
 if (hitBase)
 {
  coneCoordinateNormal.px = 0;
  coneCoordinateNormal.py = 0;
  coneCoordinateNormal.pz = -1;
 }
 else
 {
  // use the gradient of x*x + y*y - (1-z)*(1-z) = 0
  coneCoordinateNormal.px = p->px;
  coneCoordinateNormal.py = p->py;
  coneCoordinateNormal.pz = 1 - p->pz;

  // the exact tip has no unique normal; use the upward direction there
  if (coneCoordinateNormal.px*coneCoordinateNormal.px +
      coneCoordinateNormal.py*coneCoordinateNormal.py +
      coneCoordinateNormal.pz*coneCoordinateNormal.pz < TOL*TOL)
   coneCoordinateNormal.pz = 1;
 }
 coneCoordinateNormal.pw = 0;

 // put the intersection point and normal back into world coordinates
 matVecMult(cone->T, p);
 normalTransform(&coneCoordinateNormal, n, cone);
}

void cubeIntersect(struct object3D *cube, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n)
{
 // the canonical cube occupies x, y, z in [0, 1]
 // find where the ray enters and exits each pair of parallel faces
 // put the ray into canonical cube coordinates
 struct ray3D transformedRay;
 rayTransform(ray, &transformedRay, cube);

 // store the ray components in arrays so we can check all three axes in the same loop
 double origin[3] = {transformedRay.p0.px, transformedRay.p0.py, transformedRay.p0.pz};
 double direction[3] = {transformedRay.d.px, transformedRay.d.py, transformedRay.d.pz};
 // nearLambda is the latest entry into the cube, and farLambda is the earliest exit
 // also keep the axis and sign needed to calculate the surface normal
 double nearLambda = -INFINITY;
 double farLambda = INFINITY;
 int nearAxis = -1;
 int farAxis = -1;
 double nearSign = 0;
 double farSign = 0;

 // each axis restricts the allowable interval of intersection distances
 for (int axis = 0; axis < 3; axis++)
 {
  // a parallel ray cannot hit these faces unless it already lies within their boundaries
  if (fabs(direction[axis]) < TOL)
  {
   if (origin[axis] < 0 || origin[axis] > 1)
   {
    *lambda = -1;
    return;
   }
   continue;
  }

  // solve for intersections with the planes at coordinate 0 and 1
  double entry = (0 - origin[axis]) / direction[axis];
  double exit = (1 - origin[axis]) / direction[axis];
  double entrySign = -1;
  double exitSign = 1;

  // if the ray is travelling in the negative direction, swap the entry and exit faces
  if (entry > exit)
  {
   double temporary = entry;
   entry = exit;
   exit = temporary;
   entrySign = 1;
   exitSign = -1;
  }

  // keep the most restrictive entry and exit distances from the axes checked so far
  if (entry > nearLambda)
  {
   nearLambda = entry;
   nearAxis = axis;
   nearSign = entrySign;
  }
  if (exit < farLambda)
  {
   farLambda = exit;
   farAxis = axis;
   farSign = exitSign;
  }
  // when the entry is after the exit, there is no common interval inside the cube
  if (nearLambda > farLambda)
  {
   *lambda = -1;
   return;
  }
 }

 // prefer the entry face if it is in front of the ray origin
 double lambdaCandidate;
 int normalAxis;
 double normalSign;
 if (nearLambda >= TOL)
 {
  lambdaCandidate = nearLambda;
  normalAxis = nearAxis;
  normalSign = nearSign;
 }
 // otherwise, a ray starting inside the cube can still hit an exit face
 else if (farLambda >= TOL && farAxis >= 0)
 {
  // when the ray starts inside the cube, return the exit face
  lambdaCandidate = farLambda;
  normalAxis = farAxis;
  normalSign = farSign;
 }
 else
 {
  *lambda = -1;
  return;
 }

 // calculate the intersection point in canonical cube coordinates
 *lambda = lambdaCandidate;
 p->px = origin[0] + lambdaCandidate*direction[0];
 p->py = origin[1] + lambdaCandidate*direction[1];
 p->pz = origin[2] + lambdaCandidate*direction[2];
 p->pw = 1;

 // the normal points along the axis of the face we hit, with the outward sign
 struct point3D cubeCoordinateNormal = {0, 0, 0, 0};
 if (normalAxis == 0) cubeCoordinateNormal.px = normalSign;
 if (normalAxis == 1) cubeCoordinateNormal.py = normalSign;
 if (normalAxis == 2) cubeCoordinateNormal.pz = normalSign;

 // put the intersection point and normal back into world coordinates
 matVecMult(cube->T, p);
 normalTransform(&cubeCoordinateNormal, n, cube);
}

void icosahedronIntersect(struct object3D *icosahedron, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n)
{
 // we consider the canonical icosahedron to be centered at the origin, with all edges of length 1
 // the twelve vertices are given using the golden ratio, and the twenty faces are triangles
 double phi = (1 + sqrt(5.0)) / 2;
 double vertices[12][3] = {
  {0, .5, phi/2}, {0, -.5, phi/2}, {0, .5, -phi/2}, {0, -.5, -phi/2},
  {.5, phi/2, 0}, {-.5, phi/2, 0}, {.5, -phi/2, 0}, {-.5, -phi/2, 0},
  {phi/2, 0, .5}, {-phi/2, 0, .5}, {phi/2, 0, -.5}, {-phi/2, 0, -.5}
 };

 // each row contains the indices of one triangle's vertices, ordered for outward normals
 static const int faces[20][3] = {
  {0, 1, 8}, {0, 4, 5}, {0, 5, 9}, {0, 8, 4}, {0, 9, 1},
  {1, 6, 8}, {5, 4, 2}, {7, 1, 9}, {7, 3, 6}, {7, 6, 1},
  {10, 2, 4}, {10, 3, 2}, {10, 4, 8}, {10, 6, 3}, {10, 8, 6},
  {11, 2, 3}, {11, 3, 7}, {11, 5, 2}, {11, 7, 9}, {11, 9, 5}
 };

 // put the ray into canonical icosahedron coordinates
 struct ray3D transformedRay;
 rayTransform(ray, &transformedRay, icosahedron);

 double bestLambda = -1;
 struct point3D bestNormal;

 // check whether the ray intersects each of the twenty triangular faces
 for (int i = 0; i < 20; i++)
 {
  const double *v0 = vertices[faces[i][0]];
  const double *v1 = vertices[faces[i][1]];
  const double *v2 = vertices[faces[i][2]];

  // calculate the two edges from the first vertex and their cross product (the face normal)
  struct point3D e1 = {v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2], 0};
  struct point3D e2 = {v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2], 0};
  struct point3D faceNormal = {
   e1.py*e2.pz - e1.pz*e2.py,
   e1.pz*e2.px - e1.px*e2.pz,
   e1.px*e2.py - e1.py*e2.px,
   0
  };

  // first intersect the ray with the infinite plane containing the triangle
  double dDotN = dot(&transformedRay.d, &faceNormal);
  if (fabs(dDotN) < TOL) continue;

  struct point3D originToFace = {
   v0[0] - transformedRay.p0.px,
   v0[1] - transformedRay.p0.py,
   v0[2] - transformedRay.p0.pz,
   0
  };
  double candidate = dot(&originToFace, &faceNormal) / dDotN;

  // ignore intersections behind the ray origin or farther than a previous hit
  if (candidate < TOL || (bestLambda >= 0 && candidate >= bestLambda)) continue;

  // calculate the intersection point on the plane
  struct point3D hit = {
   transformedRay.p0.px + candidate*transformedRay.d.px,
   transformedRay.p0.py + candidate*transformedRay.d.py,
   transformedRay.p0.pz + candidate*transformedRay.d.pz,
   1
  };

  // check all three edges using the cross-product test from the notes
  // the intersection is inside if (edge x vector_to_hit) dot faceNormal is nonnegative for every edge
  int inside = 1;
  for (int j = 0; j < 3; j++)
  {
   const double *start = vertices[faces[i][j]];
   const double *end = vertices[faces[i][(j+1)%3]];
   struct point3D edge = {end[0]-start[0], end[1]-start[1], end[2]-start[2], 0};
   struct point3D toHit = {hit.px-start[0], hit.py-start[1], hit.pz-start[2], 0};
   struct point3D edgeCross = {
    edge.py*toHit.pz - edge.pz*toHit.py,
    edge.pz*toHit.px - edge.px*toHit.pz,
    edge.px*toHit.py - edge.py*toHit.px,
    0
   };
   if (dot(&edgeCross, &faceNormal) < -TOL)
   {
    inside = 0;
    break;
   }
  }

  // remember the closest face that contains the plane intersection
  if (inside)
  {
   bestLambda = candidate;
   bestNormal = faceNormal;
  }
 }

 // return -1 if none of the triangular faces was intersected
 if (bestLambda < 0)
 {
  *lambda = -1;
  return;
 }

 // calculate the canonical intersection point, then transform it and the normal into world coordinates
 *lambda = bestLambda;
 p->px = transformedRay.p0.px + bestLambda*transformedRay.d.px;
 p->py = transformedRay.p0.py + bestLambda*transformedRay.d.py;
 p->pz = transformedRay.p0.pz + bestLambda*transformedRay.d.pz;
 p->pw = 1;
 matVecMult(icosahedron->T, p);
 normalTransform(&bestNormal, n, icosahedron);

}

void insertPLS(struct pointLS *l, struct pointLS **list)
{
 if (l==NULL) return;
 // insert this point light into the scene
 if (*(list)==NULL)
 {
  *(list)=l;
  (*(list))->next=NULL;
 }
 else
 {
  l->next=(*(list))->next;
  (*(list))->next=l;
 }

}

void invert(double *T, double *Tinv)
{

 double *U, *s, *V, *rv1;
 int singFlag, i;

 U=NULL;
 s=NULL;
 V=NULL;
 rv1=NULL;
 singFlag=0;

 SVD(T,4,4,&U,&s,&V,&rv1);
 if (U==NULL||s==NULL||V==NULL)
 {
  fprintf(stderr,"Error: Matrix not invertible for this object, returning identity\n");
  memcpy(Tinv,eye4x4,16*sizeof(double));
  return;
 }

 for (i=0;i<4;i++) if (*(s+i)<1e-9) singFlag=1;
 if (singFlag)
 {
  fprintf(stderr,"Error: Transformation matrix is singular, returning identity\n");
  memcpy(Tinv,eye4x4,16*sizeof(double));
  return;
 }

 InvertMatrix(U,s,V,4,Tinv);

 free(U);
 free(s);
 free(V);
}

void RotateXMat(double T[4][4], double theta)
{

 double R[4][4];
 memset(&R[0][0],0,16*sizeof(double));

 R[0][0]=1.0;
 R[1][1]=cos(theta);
 R[1][2]=-sin(theta);
 R[2][1]=sin(theta);
 R[2][2]=cos(theta);
 R[3][3]=1.0;

 matMult(R,T);
}

void RotateX(struct object3D *o, double theta)
{

 double R[4][4];
 memset(&R[0][0],0,16*sizeof(double));

 R[0][0]=1.0;
 R[1][1]=cos(theta);
 R[1][2]=-sin(theta);
 R[2][1]=sin(theta);
 R[2][2]=cos(theta);
 R[3][3]=1.0;

 matMult(R,o->T);
}

void RotateYMat(double T[4][4], double theta)
{

 double R[4][4];
 memset(&R[0][0],0,16*sizeof(double));

 R[0][0]=cos(theta);
 R[0][2]=sin(theta);
 R[1][1]=1.0;
 R[2][0]=-sin(theta);
 R[2][2]=cos(theta);
 R[3][3]=1.0;

 matMult(R,T);
}

void RotateY(struct object3D *o, double theta)
{

 double R[4][4];
 memset(&R[0][0],0,16*sizeof(double));

 R[0][0]=cos(theta);
 R[0][2]=sin(theta);
 R[1][1]=1.0;
 R[2][0]=-sin(theta);
 R[2][2]=cos(theta);
 R[3][3]=1.0;

 matMult(R,o->T);
}

void RotateZMat(double T[4][4], double theta)
{

 double R[4][4];
 memset(&R[0][0],0,16*sizeof(double));

 R[0][0]=cos(theta);
 R[0][1]=-sin(theta);
 R[1][0]=sin(theta);
 R[1][1]=cos(theta);
 R[2][2]=1.0;
 R[3][3]=1.0;

 matMult(R,T);
}

void RotateZ(struct object3D *o, double theta)
{

 double R[4][4];
 memset(&R[0][0],0,16*sizeof(double));

 R[0][0]=cos(theta);
 R[0][1]=-sin(theta);
 R[1][0]=sin(theta);
 R[1][1]=cos(theta);
 R[2][2]=1.0;
 R[3][3]=1.0;

 matMult(R,o->T);
}

void TranslateMat(double T[4][4], double tx, double ty, double tz)
{

 double tr[4][4];
 memset(&tr[0][0],0,16*sizeof(double));

 tr[0][0]=1.0;
 tr[1][1]=1.0;
 tr[2][2]=1.0;
 tr[0][3]=tx;
 tr[1][3]=ty;
 tr[2][3]=tz;
 tr[3][3]=1.0;

 matMult(tr,T);
}

void Translate(struct object3D *o, double tx, double ty, double tz)
{

 double tr[4][4];
 memset(&tr[0][0],0,16*sizeof(double));

 tr[0][0]=1.0;
 tr[1][1]=1.0;
 tr[2][2]=1.0;
 tr[0][3]=tx;
 tr[1][3]=ty;
 tr[2][3]=tz;
 tr[3][3]=1.0;

 matMult(tr,o->T);
}

void ScaleMat(double T[4][4], double sx, double sy, double sz)
{

 double S[4][4];
 memset(&S[0][0],0,16*sizeof(double));

 S[0][0]=sx;
 S[1][1]=sy;
 S[2][2]=sz;
 S[3][3]=1.0;

 matMult(S,T);
}

void Scale(struct object3D *o, double sx, double sy, double sz)
{

 double S[4][4];
 memset(&S[0][0],0,16*sizeof(double));

 S[0][0]=sx;
 S[1][1]=sy;
 S[2][2]=sz;
 S[3][3]=1.0;

 matMult(S,o->T);
}

void printmatrix(double mat[4][4])
{
 fprintf(stderr,"Matrix contains:\n");
 fprintf(stderr,"%f %f %f %f\n",mat[0][0],mat[0][1],mat[0][2],mat[0][3]);
 fprintf(stderr,"%f %f %f %f\n",mat[1][0],mat[1][1],mat[1][2],mat[1][3]);
 fprintf(stderr,"%f %f %f %f\n",mat[2][0],mat[2][1],mat[2][2],mat[2][3]);
 fprintf(stderr,"%f %f %f %f\n",mat[3][0],mat[3][1],mat[3][2],mat[3][3]);
}

struct view *setupView(struct point3D *e, struct point3D *g, struct point3D *up, double f, double wl, double wt, double wsize)
{
 
 struct view *c;
 struct point3D *u, *v;

 u=v=NULL;

 c=(struct view *)calloc(1,sizeof(struct view));
 if (c==NULL)
 {
  fprintf(stderr,"Out of memory setting up camera model!\n");
  return(NULL);
 }

 c->e.px=e->px;
 c->e.py=e->py;
 c->e.pz=e->pz;
 c->e.pw=1;

 c->w.px=-g->px;
 c->w.py=-g->py;
 c->w.pz=-g->pz;
 c->w.pw=1;
 normalize(&c->w);

 u=cross(&c->w, up);
 normalize(u);
 c->u.px=u->px;
 c->u.py=u->py;
 c->u.pz=u->pz;
 c->u.pw=1;

 v=cross(&c->u, &c->w);
 normalize(v);
 c->v.px=v->px;
 c->v.py=v->py;
 c->v.pz=v->pz;
 c->v.pw=1;

 c->f=f;
 c->wl=wl;
 c->wt=wt;
 c->wsize=wsize;

 c->C2W[0][0]=c->u.px;
 c->C2W[1][0]=c->u.py;
 c->C2W[2][0]=c->u.pz;
 c->C2W[3][0]=0;

 c->C2W[0][1]=c->v.px;
 c->C2W[1][1]=c->v.py;
 c->C2W[2][1]=c->v.pz;
 c->C2W[3][1]=0;

 c->C2W[0][2]=c->w.px;
 c->C2W[1][2]=c->w.py;
 c->C2W[2][2]=c->w.pz;
 c->C2W[3][2]=0;

 c->C2W[0][3]=c->e.px;
 c->C2W[1][3]=c->e.py;
 c->C2W[2][3]=c->e.pz;
 c->C2W[3][3]=1;

 c->W2C[0][0]=c->u.px;
 c->W2C[1][0]=c->v.px;
 c->W2C[2][0]=c->w.px;
 c->W2C[3][0]=0;

 c->W2C[0][1]=c->u.py;
 c->W2C[1][1]=c->v.py;
 c->W2C[2][1]=c->w.py;
 c->W2C[3][1]=0;

 c->W2C[0][2]=c->u.pz;
 c->W2C[1][2]=c->v.pz;
 c->W2C[2][2]=c->w.pz;
 c->W2C[3][2]=0;

 c->W2C[0][3]=-dot(&c->u,&c->e);
 c->W2C[1][3]=-dot(&c->v,&c->e);
 c->W2C[2][3]=-dot(&c->w,&c->e);
 c->W2C[3][3]=1;

 free(u);
 free(v);
 return(c);
}

struct image *newImage(int size_x, int size_y)
{

 struct image *im;

 im=(struct image *)calloc(1,sizeof(struct image));
 if (im!=NULL)
 {
  im->rgbdata=NULL;
  im->sx=size_x;
  im->sy=size_y;
  im->rgbdata=(void *)calloc(size_x*size_y*3,sizeof(unsigned char));
  if (im->rgbdata!=NULL) return(im);
 }
 fprintf(stderr,"Unable to allocate memory for new image\n");
 return(NULL);
}

void imageOutput(struct image *im, const char *filename)
{

 FILE *f;

 if (im!=NULL)
  if (im->rgbdata!=NULL)
  {
   f=fopen(filename,"wb+");
   if (f==NULL)
   {
    fprintf(stderr,"Unable to open file %s for output! No image written\n",filename);
    return;
   }
   fprintf(f,"P6\n");
   fprintf(f,"# Output from spookyLights3D\n");
   fprintf(f,"%d %d\n",im->sx,im->sy);
   fprintf(f,"255\n");
   fwrite((unsigned char *)im->rgbdata,im->sx*im->sy*3*sizeof(unsigned char),1,f);
   fclose(f);
   return;
  }
 fprintf(stderr,"imageOutput(): Specified image is empty. Nothing output\n");
}

void deleteImage(struct image *im)
{

 if (im!=NULL)
 {
  if (im->rgbdata!=NULL) free(im->rgbdata);
  free(im);
 }

}
void cleanup(struct object3D *o_list, struct pointLS *l_list)
{
 // release the scene objects and point lights
 while (o_list != NULL)
 {
  struct object3D *next = o_list->next;
  free(o_list);
  o_list = next;
 }
 while (l_list != NULL)
 {
  struct pointLS *next = l_list->next;
  free(l_list);
  l_list = next;
 }
}
