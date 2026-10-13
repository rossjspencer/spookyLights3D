#include "utils.h"

void buildScene(void)
{
 struct object3D *o;
 struct pointLS *l;
 struct point3D p;

 // create two transformed spheres and a floor
 o=newSphere(.05,.95,.35,.35,1,.25,.25,6);
 Scale(o,1.5,.75,.75);
 RotateZ(o,PI/4);
 Translate(o,2.0,2.5,1.5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newSphere(.05,.95,.95,.75,.75,.95,.55,6);
 Scale(o,.95,1.65,.65);
 RotateZ(o,-PI/1.5);
 Translate(o,-2.2,1.75,1.35);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newPlane(.05,.75,.05,.05,.55,.8,.75,2);
 Scale(o,11,11,11);
 RotateZ(o,PI/4);
 RotateX(o,PI/2);
 Translate(o,0,-4,5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 // position a white point light above the objects
 p.px=0;
 p.py=25.5;
 p.pz=-3.5;
 p.pw=1;
 l=newPLS(&p,.95,.95,.95);
 insertPLS(l,&light_list);
}
