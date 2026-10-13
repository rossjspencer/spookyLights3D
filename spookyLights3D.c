// spookyLights3D ray tracer
// derived in part from educational code by f. j. estrada; see attribution.md
#include "utils.h"

// scene objects, lights, and maximum reflection depth
struct object3D *object_list;
struct pointLS *light_list;
int MAX_DEPTH;



void rtShade(struct object3D *obj, struct point3D *p, struct point3D *n, struct ray3D *ray, int depth, struct colourRGB *col)
{
 // shade an intersection with ambient, diffuse, specular, and reflected light
 struct colourRGB tmp_col;
 double R,G,B;		

 tmp_col.R=0;
 tmp_col.G=0;
 tmp_col.B=0;

 double ra = obj->alb.ra;
 double rd = obj->alb.rd;
 double rs = obj->alb.rs;
 double rg = obj->alb.rg;

 R = obj->col.R;
 G = obj->col.G;
 B = obj->col.B;

 // prepare the c direction since it does not change between light sources
 struct point3D c;
 c.px = -ray->d.px;
 c.py = -ray->d.py;
 c.pz = -ray->d.pz;
 c.pw = 0;
 normalize(&c);

 // add the ambient phong contribution
 tmp_col.R += ra * R;
 tmp_col.G += ra * G;
 tmp_col.B += ra * B;

 if (rg > 0 && depth < MAX_DEPTH)
 {
  // compute the perfect reflection direction of the incoming ray
  struct point3D reflectedDirection;
  double dDotN = dot(&ray->d, n);
  reflectedDirection.px = ray->d.px - 2 * dDotN * n->px;
  reflectedDirection.py = ray->d.py - 2 * dDotN * n->py;
  reflectedDirection.pz = ray->d.pz - 2 * dDotN * n->pz;
  reflectedDirection.pw = 0;
  normalize(&reflectedDirection);

  // create a reflected ray starting at p
  struct ray3D reflectedRay;
  initRay(&reflectedRay, p, &reflectedDirection);

  struct colourRGB reflectedColour;

  // recursively trace that reflected ray
  rayTrace(&reflectedRay, depth + 1, &reflectedColour, obj);

  // if it hits something, add rg times that returned colour
  if (reflectedColour.R >= 0)
  {
   tmp_col.R += rg * reflectedColour.R;
   tmp_col.G += rg * reflectedColour.G;
   tmp_col.B += rg * reflectedColour.B;
  }
 }

 struct pointLS *current = light_list; // we assume the global pointer is the head of the list

 // for every light, apply its specular and diffuse contributions
 while (current != NULL)
 {
  struct ray3D shadowRay;

  struct point3D s;
  s.px = current->p0.px - p->px;
  s.py = current->p0.py - p->py;
  s.pz = current->p0.pz - p->pz;
  s.pw = 0;

  // don't normalize s yet so that lambda = 1 represents the full ray between p and the light source
  shadowRay.d = s;
  shadowRay.p0 = *p;

  double shadowLambda = -1;
  struct object3D *shadowObj = NULL;
  struct point3D shadowPoint, shadowNormal;

  findFirstHit(&shadowRay, &shadowLambda, obj, &shadowObj, &shadowPoint, &shadowNormal);

  // if the shadow ray collides with an object before reaching the light source
  if (shadowLambda > 0 && shadowLambda < 1)
  {
   // do not add diffuse or specular components for this source
   current = current->next;
   continue;
  }

  // normalize s for the remaining calculations
  normalize(&s);

  // prepare and calculate m
  struct point3D m;
  double sDotN = dot(&s,n);
  m.px = n->px * 2 * sDotN - s.px;
  m.py = n->py * 2 * sDotN - s.py;
  m.pz = n->pz * 2 * sDotN - s.pz;
  m.pw = 0;
  normalize(&m);

  // we assume these colour values correspond to the Id and Is intensities for the Phong model
  double lightR = current->col.R;
  double lightG = current->col.G;
  double lightB = current->col.B;

  // add diffuse component
  tmp_col.R += R * rd * lightR * fmax(0, dot(n, &s));
  tmp_col.G += G * rd * lightG * fmax(0, dot(n, &s));
  tmp_col.B += B * rd * lightB * fmax(0, dot(n, &s));

  // add specular component
  tmp_col.R += rs * lightR * pow(fmax(0, dot(&c, &m)), obj->shinyness);
  tmp_col.G += rs * lightG * pow(fmax(0, dot(&c, &m)), obj->shinyness);
  tmp_col.B += rs * lightB * pow(fmax(0, dot(&c, &m)), obj->shinyness);

  current = current->next;
 }

 // store the final colour for this ray
 col->R = tmp_col.R;
 col->G = tmp_col.G;
 col->B = tmp_col.B;
}

void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p, struct point3D *n)
{
 // find the nearest valid object intersection
 struct object3D *current = object_list;
 struct point3D tempPoint;
 struct point3D tempNormal;
 double bestLambda = -1;
 double currentLambda;

 // check intersection with every object in the scene
 while (current != NULL)
 {
  currentLambda = -1;

  // ignore self intersection
  // skipping the source object avoids numerical self-intersections
  if (current == Os)
  {
   current = current->next;
   continue;
  }

  // check intersection with the given object
  current->intersect(current, ray, &currentLambda, &tempPoint, &tempNormal);

  // update the object pointer, normal, and intersection only if the found lambda is an improvement
  // also update pointer if bestLambda is still -1 and currentLambda found something
  if (currentLambda > 0 && (currentLambda < bestLambda || bestLambda < 0))
  {
   *obj = current;
   *p = tempPoint;
   *n = tempNormal;
   bestLambda = currentLambda;
  }

  current = current->next;
 }

 *lambda = bestLambda;
}

void rayTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os)
{
 // trace a ray and shade its nearest intersection
 double lambda;	
 struct object3D *obj;
 struct point3D p;
 struct point3D n;
 struct colourRGB I;

 if (depth>MAX_DEPTH)	// Max recursion depth reached. Return invalid colour.
 {
  col->R=-1;
  col->G=-1;
  col->B=-1;
  return;
 }


 // find the closest object intersected by this ray
 obj = NULL;
 lambda = -1;
 findFirstHit(ray, &lambda, Os, &obj, &p, &n);

 // if the ray did not hit anything, return an invalid colour
 if (obj == NULL || lambda < 0)
 {
  col->R=-1;
  col->G=-1;
  col->B=-1;
  return;
 }

 // shade the closest intersection point
 rtShade(obj, &p, &n, ray, depth, &I);

 // return the colour found by the shading function
 col->R=I.R;
 col->G=I.G;
 col->B=I.B;
}

int main(int argc, char *argv[])
{
 // set up the camera and render the image
 struct image *im;
 struct view *cam;
 int sx;	
 const char *output_name;
 struct point3D e;	
 struct point3D g;
 struct point3D up;
 double du, dv;		
 struct point3D pc,d;	
 struct ray3D ray;	
 struct colourRGB col;	
 struct colourRGB background;   // Background colour
 int i,j;			// Counters for pixel coordinates
 unsigned char *rgbIm;

 if (argc != 4)
 {
  fprintf(stderr, "usage: %s <image_size> <reflection_depth> <output.ppm>\n", argv[0]);
  return 1;
 }
 sx = atoi(argv[1]);
 MAX_DEPTH = atoi(argv[2]);
 output_name = argv[3];
 if (sx <= 0 || MAX_DEPTH < 0)
 {
  fprintf(stderr, "image size must be positive and reflection depth must be nonnegative\n");
  return 1;
 }
 fprintf(stderr, "rendering %dx%d image with reflection depth %d\n", sx, sx, MAX_DEPTH);

 object_list=NULL;
 light_list=NULL;

 // allocate the output image
 im=newImage(sx, sx);
 if (!im)
 {
  fprintf(stderr,"unable to allocate output image\n");
  return 1;
 }
 else rgbIm=(unsigned char *)im->rgbdata;
 buildScene(); // add scene objects and lights

 // position the camera behind the scene
 e.px=0;
 e.py=0;
 e.pz=-1;
 e.pw=1;

 // point the camera toward the origin
 g.px=-e.px;
 g.py=-e.py;
 g.pz=-e.pz;
 g.pw=1;

 // use the y axis as the upward direction
 up.px=0;
 up.py=1;
 up.pz=0;
 up.pw=1;

 // use a 4-unit camera viewing window
 cam=setupView(&e, &g, &up, -1, -2, 2, 4);

 if (cam==NULL)
 {
  fprintf(stderr,"unable to set up the camera\n");
  cleanup(object_list,light_list);
  deleteImage(im);
  return 1;
 }

 // set a black background
 background.R=0;
 background.G=0;
 background.B=0;

 // calculate the distance between adjacent pixel centers
 du=cam->wsize/sx;
 dv=-cam->wsize/sx;


 fprintf(stderr,"rendering rows...\n");
 for (j=0;j<sx;j++)		// For each of the pixels in the image
 {
  if (j % (sx / 10 + 1) == 0) fprintf(stderr,"  row %d/%d\n",j,sx);
  for (i=0;i<sx;i++)
  {
   // find the center of each pixel in camera coordinates
   pc.px=cam->wl+(i+0.5)*du;
   pc.py=cam->wt+(j+0.5)*dv;
   pc.pz=cam->f;
   pc.pw=1;

   // put the pixel center into world coordinates
   matVecMult(cam->C2W,&pc);

   // get the direction from the camera center through the pixel center
   d.px=pc.px-cam->e.px;
   d.py=pc.py-cam->e.py;
   d.pz=pc.pz-cam->e.pz;
   d.pw=0;
   normalize(&d);

   // initialize and trace the primary ray
   initRay(&ray,&cam->e,&d);
   rayTrace(&ray,0,&col,NULL);

   // if the ray hit nothing, use the background colour
   if (col.R<0 || col.G<0 || col.B<0)
   {
    col.R=background.R;
    col.G=background.G;
    col.B=background.B;
   }

   // make sure the final colour stays in [0,1]
   if (col.R<0) col.R=0;
   if (col.R>1) col.R=1;
   if (col.G<0) col.G=0;
   if (col.G>1) col.G=1;
   if (col.B<0) col.B=0;
   if (col.B>1) col.B=1;

   // store the colour in the output image
   rgbIm[(j*sx+i)*3]=255*col.R;
   rgbIm[(j*sx+i)*3+1]=255*col.G;
   rgbIm[(j*sx+i)*3+2]=255*col.B;
  } // end for i
 } // end for j

 fprintf(stderr,"done\n");

 // write the image as a binary ppm file
 imageOutput(im,output_name);

 // release the scene and image memory
 cleanup(object_list,light_list);
 deleteImage(im);				
 free(cam);					
 return 0;
}

