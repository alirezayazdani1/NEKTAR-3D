//                                                                        //
//   Author:    T.Warburton                                               //
//   Design:    T.Warburton && S.Sherwin                                  //
//   Date  :    12/4/96                                                   //
//                                                                        //
//   Copyright notice:  This code shall not be replicated or used without //
//                      the permission of the author.                     //
//                                                                        //
/**************************************************************************/

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <polylib.h>
#include "veclib.h"
#include "hotel.h"

using namespace polylib;

#define TANTOL 1e-10
/* new atan2 function to stop Nan on atan(0,0)*/
static double atan2_proof (double x, double y)
{
  if (fabs(x) + fabs(y) > TANTOL) return (atan2(x,y));
  else return (0.);
}
#define atan2 atan2_proof



typedef struct point    {  /* A 2-D point  */
  double  x,y;             /* coordinate   */
} Point;


double find_spiral_theta(Curve *curve, double x0, double y0, double z0);

void genCylinder(Curve *curve, double *x, double *y, double *z, int q){
  register int i;
  double   *f,theta,phi,cp,sp,ct,st;
  
  f = dvector(0,q-1);

  /* translate co-ordinates so that given point is origin */
  dsadd(q, -curve->info.cyl.xc, x, 1, x, 1);
  dsadd(q, -curve->info.cyl.yc, y, 1, y, 1);
  dsadd(q, -curve->info.cyl.zc, z, 1, z, 1);
    
  /* rotate co-ordinates so that cylinder axis is aligned with z axis */
  phi = atan2(curve->info.cyl.ay, curve->info.cyl.ax);
  cp  = cos(phi); sp = sin(phi);
  drot(1,&curve->info.cyl.ax,1,&curve->info.cyl.ay,1,cp,sp);
  theta = atan2(curve->info.cyl.ax,curve->info.cyl.az);
  ct    = cos(theta); st = sin(theta);

  drot(q,x,1,y,1,cp,sp);
  drot(q,z,1,x,1,ct,st);

  /* move co-ordinate out to surface value */
  for(i = 0; i < q; ++i){
    theta = atan2(y[i],x[i]);
    x[i]  = curve->info.cyl.radius*cos(theta);
    y[i]  = curve->info.cyl.radius*sin(theta);
  }
  
  /* rotate back */
  drot(q,z,1,x,1,ct,-st);
  drot(q,x,1,y,1,cp,-sp);
  
  /* translate co-ordinates back to original position */
  dsadd(q, curve->info.cyl.xc, x, 1, x, 1);
  dsadd(q, curve->info.cyl.yc, y, 1, y, 1);
  dsadd(q, curve->info.cyl.zc, z, 1, z, 1);
    
  free(f);
}

void genCone(Curve *curve, double *x, double *y, double *z, int q){
  register int i;
  double   *f,theta,phi,cp,sp,ct,st;
  
  f = dvector(0,q-1);

  /* translate co-ordinates so that apex is at origin */
  dsadd(q, -curve->info.cone.xc, x, 1, x, 1);
  dsadd(q, -curve->info.cone.yc, y, 1, y, 1);
  dsadd(q, -curve->info.cone.zc, z, 1, z, 1);
    
  /* rotate co-ordinates so that cone axis is aligned with z axis */
  phi   = atan2(curve->info.cone.ay, curve->info.cone.ax);
  cp = cos(phi); sp = sin(phi);
  drot(1,&curve->info.cone.ax,1,&curve->info.cone.ay,1,cp,sp);
  theta = atan2(curve->info.cone.ax,curve->info.cone.az);
  ct = cos(theta); st = sin(theta);
  
  drot(q,x,1,y,1,cp,sp);
  drot(q,z,1,x,1,ct,st);

  /* move co-ordinate out to surface value */
  for(i = 0; i < q; ++i){
    theta = atan2(y[i],x[i]);
    x[i]  = fabs(z[i])*curve->info.cone.alpha*cos(theta);
    y[i]  = fabs(z[i])*curve->info.cone.alpha*sin(theta);
  }
  
  /* rotate back */
  drot(q,z,1,x,1,ct,-st);
  drot(q,x,1,y,1,cp,-sp);
  
  /* translate co-ordinates back to original position */
  dsadd(q, curve->info.cone.xc, x, 1, x, 1);
  dsadd(q, curve->info.cone.yc, y, 1, y, 1);
  dsadd(q, curve->info.cone.zc, z, 1, z, 1);
    
  free(f);
}

void genSphere(Curve *curve, double *x, double *y, double *z, int q){
  register int i;
  double   *f,theta,phi;
  
  f = dvector(0,q-1);

  /* translate co-ordinates so that apex is at origin */
  dsadd(q, -curve->info.sph.xc, x, 1, x, 1);
  dsadd(q, -curve->info.sph.yc, y, 1, y, 1);
  dsadd(q, -curve->info.sph.zc, z, 1, z, 1);
    
  /* move co-ordinate out to surface value */
  for(i = 0; i < q; ++i){
    phi   = atan(z[i]/sqrt(x[i]*x[i]+y[i]*y[i]));
    theta = atan2(y[i],x[i]);
    z[i]  = curve->info.sph.radius*sin(phi);
    x[i]  = y[i]  = curve->info.sph.radius*cos(phi);
    x[i]  *= cos(theta);
    y[i]  *= sin(theta);
  }
  
  /* translate co-ordinates back to original position */
  dsadd(q, curve->info.sph.xc, x, 1, x, 1);
  dsadd(q, curve->info.sph.yc, y, 1, y, 1);
  dsadd(q, curve->info.sph.zc, z, 1, z, 1);
    
  free(f);
}

void genSheet(Curve *curve, double *x, double *y, double *z, int q){
  register int i;
  double   *f,theta,phi,cp,sp,ct,st,rtmp;
  
  f = dvector(0,q-1);

  /* translate co-ordinates so that given point is origin */
  dsadd(q, -curve->info.she.xc, x, 1, x, 1);
  dsadd(q, -curve->info.she.yc, y, 1, y, 1);
  dsadd(q, -curve->info.she.zc, z, 1, z, 1);
    
  /* rotate co-ordinates so that cylinder axis is aligned with z axis */
  phi   = atan2(curve->info.she.ay, curve->info.she.ax);
  cp = cos(phi); sp = sin(phi);
  drot(1,&curve->info.she.ax,1,&curve->info.she.ay,1,cp,sp);
  theta = atan2(curve->info.she.ax,curve->info.she.az);
  ct = cos(theta); st = sin(theta);

  drot(q,x,1,y,1,cp,sp);
  drot(q,z,1,x,1,ct,st);

  /* move co-ordinate out to surface value */
  for(i = 0; i < q; ++i){
    theta = z[i]*curve->info.she.twist+curve->info.she.zerotwistz;
    rtmp = sqrt(x[i]*x[i]+y[i]*y[i]);
    x[i]  = rtmp*cos(theta);
    y[i]  = rtmp*sin(theta);
  }
  /* rotate back */
  drot(q,z,1,x,1,ct,-st);
  drot(q,x,1,y,1,cp,-sp);
  
  /* translate co-ordinates back to original position */
  dsadd(q, curve->info.she.xc, x, 1, x, 1);
  dsadd(q, curve->info.she.yc, y, 1, y, 1);
  dsadd(q, curve->info.she.zc, z, 1, z, 1);
    
  free(f);
}

#ifdef OLDHELIX
void genSpiral(Curve *curve, double *x, double *y, double *z, int q){
  register int i;
  double   *f,theta,phi,cp,sp,ct,st,rtmp;
  
  f = dvector(0,q-1);

  /* translate co-ordinates so that given point is origin */
  dsadd(q, -curve->info.spi.xc, x, 1, x, 1);
  dsadd(q, -curve->info.spi.yc, y, 1, y, 1);
  dsadd(q, -curve->info.spi.zc, z, 1, z, 1);
    
  /* rotate co-ordinates so that cylinder axis is aligned with z axis */
  phi   = atan2(curve->info.spi.ay, curve->info.spi.ax);
  cp = cos(phi); sp = sin(phi);
  drot(1,&curve->info.spi.ax,1,&curve->info.spi.ay,1,cp,sp);
  theta = atan2(curve->info.spi.ax,curve->info.spi.az);
  ct = cos(theta); st = sin(theta);

  drot(q,x,1,y,1,cp,sp);
  drot(q,z,1,x,1,ct,st);

  /* move co-ordinate out to surface value */
  for(i = 0; i < q; ++i){
    theta = z[i]*curve->info.spi.twist+curve->info.spi.zerotwistz;
    rtmp = curve->info.spi.axialradius;
    x[i] -= rtmp*cos(theta);
    y[i] -= rtmp*sin(theta);
    phi   = atan2(y[i],x[i]);
    x[i]  = curve->info.spi.piperadius*cos(phi);
    y[i]  = curve->info.spi.piperadius*sin(phi);
    x[i] += rtmp*cos(theta);
    y[i] += rtmp*sin(theta);
  }
  /* rotate back */
  drot(q,z,1,x,1,ct,-st);
  drot(q,x,1,y,1,cp,-sp);
  
  /* translate co-ordinates back to original position */
  dsadd(q, curve->info.spi.xc, x, 1, x, 1);
  dsadd(q, curve->info.spi.yc, y, 1, y, 1);
  dsadd(q, curve->info.spi.zc, z, 1, z, 1);
    
  free(f);
}
#else


void genSpiral(Curve *curve, double *x, double *y, double *z, int q){
  register int i;
  double   *f,theta,theta0,phi,cp,sp,ct,st;
  double   axialrad,pitch,piperad;
  double   cx,cy,cz;
  
  f = dvector(0,q-1);

  /* translate co-ordinates so that given point is origin */
  dsadd(q, -curve->info.spi.xc, x, 1, x, 1);
  dsadd(q, -curve->info.spi.yc, y, 1, y, 1);
  dsadd(q, -curve->info.spi.zc, z, 1, z, 1);
    
  /* rotate co-ordinates so that cylinder axis is aligned with z axis */
  phi   = atan2(curve->info.spi.ay, curve->info.spi.ax);
  cp    = cos(phi);   sp = sin(phi);
  drot(1,&curve->info.spi.ax,1,&curve->info.spi.ay,1,cp,sp);
  theta = atan2(curve->info.spi.ax,curve->info.spi.az);
  ct    = cos(theta); st = sin(theta);

  drot(q,x,1,y,1,cp,sp);
  drot(q,z,1,x,1,ct,st);

  phi = 0.5*M_PI -
    atan(curve->info.spi.pitch/(2*M_PI*curve->info.spi.axialradius));

  piperad  = curve->info.spi.piperadius;
  axialrad = curve->info.spi.axialradius;
  pitch    = curve->info.spi.pitch;
  
  /* move co-ordinate out to surface value */
  for(i = 0; i < q; ++i){
    
    /* intially find theta0 assuming in first winding */
    theta0  = find_spiral_theta(curve,x[i],y[i],z[i]);
    cz      = 0.0;

    /* if theta0 gives a helix point which is outside pipe radius */
    /* assume that we are in the next winding and so re-calculate */
    while(fabs(z[i] - cz - pitch*theta0/(2*M_PI)) > piperad){
      cz     += pitch;
      theta0  = find_spiral_theta(curve,x[i],y[i],z[i]-cz);
    }

    cx  = axialrad*cos(theta0);
    cy  = axialrad*sin(theta0);
    cz += pitch*theta0/(2*M_PI);

    /* move section centre to origin */
    x[i] -= cx;
    y[i] -= cy;
    z[i] -= cz;
    
    /* rotate theta about z axis */
    drot(1,x+i,1,y+i,1,cos(theta0),sin(theta0));
    
    /* rotate -phi about x axis */
    drot(1,y+i,1,z+i,1,cos(-phi),sin(-phi));

    theta = atan2(y[i],x[i]);
    x[i]  = curve->info.spi.piperadius*cos(theta);
    y[i]  = curve->info.spi.piperadius*sin(theta);

    /* rotate phi about x axis */
    drot(1,y+i,1,z+i,1,cos(phi),sin(phi));

    /* rotate theta about z axis */
    drot(1,x+i,1,y+i,1,cos(-theta0),sin(-theta0));

    /* move section centre back */
    x[i] += cx;
    y[i] += cy;
    z[i] += cz;
  }
  
  /* rotate back */
  drot(q,z,1,x,1,ct,-st);
  drot(q,x,1,y,1,cp,-sp);

  /* translate co-ordinates back to original position */
  dsadd(q, curve->info.spi.xc, x, 1, x, 1);
  dsadd(q, curve->info.spi.yc, y, 1, y, 1);
  dsadd(q, curve->info.spi.zc, z, 1, z, 1);

  free(f);
}

#define TOLTHETA 1e-10
#define ITERTHETA 100
/* find the angle corresponding to the intersection of the helix and the plane
   containing the points x0,y0,z0 */

double find_spiral_theta(Curve *curve, double x0, double y0, double z0){
  register int     count = 0;
  double   dt;
  double   theta, theta0;
  double   pitch  = curve->info.spi.pitch;

  double   rad    = curve->info.spi.axialradius;
  double   A,B;


  theta = theta0 = atan2(y0,x0);
  /* correct for negative value of theta0 */
  theta0 += (fabs(z0 - pitch*theta0/(2*M_PI)) < rad)? 0:2*M_PI;

  theta = theta0;

  A = pitch*pitch/(4*M_PI*M_PI*rad);
  B = -x0*sin(theta0) + y0*cos(theta0) + z0*pitch/(2*M_PI*rad);

  dt = (-A*theta + B)/(rad+A);
 
  while((fabs(dt) > TOLTHETA)&&(count++ < ITERTHETA)){
    theta += dt;
    dt = (rad*sin(theta0 - theta) - A*theta + B)/(rad*cos(theta0-theta)+A);
  }
  if(count == ITERTHETA)
    fprintf(stderr,"Iterations failed to converge in spiral_theta\n");
  
  return theta;
}
#endif

void genTaurus(Curve *curve, double *x, double *y, double *z, int q){
  register int i;
  double   *f,theta,phi,cp,sp,ct,st,rtmp,xtmp,ytmp;
  
  f = dvector(0,q-1);

  /* translate co-ordinates so that given point is origin */
  dsadd(q, -curve->info.tau.xc, x, 1, x, 1);
  dsadd(q, -curve->info.tau.yc, y, 1, y, 1);
  dsadd(q, -curve->info.tau.zc, z, 1, z, 1);
    
  /* rotate co-ordinates so that tauinder axis is aligned with z axis */
  phi   = atan2(curve->info.tau.ay, curve->info.tau.ax);
  cp = cos(phi); sp = sin(phi);
  drot(1,&curve->info.tau.ax,1,&curve->info.tau.ay,1,cp,sp);
  theta = atan2(curve->info.tau.ax,curve->info.tau.az);
  ct = cos(theta); st = sin(theta);

  drot(q,x,1,y,1,cp,sp);
  drot(q,z,1,x,1,ct,st);

  /* move co-ordinate out to surface value */
  for(i = 0; i < q; ++i){
    theta = atan2(y[i],x[i]);
    rtmp = curve->info.tau.axialradius;
    xtmp = x[i];
    ytmp = y[i];
    x[i] = xtmp*cos(theta)+ytmp*sin(theta)-rtmp;
    y[i] = -xtmp*sin(theta)+ytmp*cos(theta);
    phi = atan2(z[i],x[i]);
    x[i] = curve->info.tau.piperadius*cos(phi)+rtmp;
    z[i] = curve->info.tau.piperadius*sin(phi);
    xtmp = x[i];
    ytmp = y[i];
    x[i]  = xtmp*cos(theta)-ytmp*sin(theta);
    y[i]  = xtmp*sin(theta)+ytmp*cos(theta);
  }
  
  /* rotate back */
  drot(q,z,1,x,1,ct,-st);
  drot(q,x,1,y,1,cp,-sp);
  
  /* translate co-ordinates back to original position */
  dsadd(q, curve->info.tau.xc, x, 1, x, 1);
  dsadd(q, curve->info.tau.yc, y, 1, y, 1);
  dsadd(q, curve->info.tau.zc, z, 1, z, 1);
    
  free(f);
}


void trans_coords(Coord *o, Coord *t, Coord *n, Coord *b, 
		  Coord *x, Coord *newx, int np, int dir){
  int i;


  // to (t,n,b)
  if(dir == 1){
    dcopy(np, x->x, 1, newx->x, 1);
    dcopy(np, x->y, 1, newx->y, 1);
    dcopy(np, x->z, 1, newx->z, 1);

    dsadd(np, -o->x[0], newx->x, 1, newx->x, 1);
    dsadd(np, -o->y[0], newx->y, 1, newx->y, 1);
    dsadd(np, -o->z[0], newx->z, 1, newx->z, 1);

    for(i=0;i<np;++i){
      newx->x[i] = t->x[0]*x->x[i] + t->y[0]*x->y[i] + t->z[0]*x->z[i];
      newx->y[i] = n->x[0]*x->x[i] + n->y[0]*x->y[i] + n->z[0]*x->z[i];
      newx->z[i] = b->x[0]*x->x[i] + b->y[0]*x->y[i] + b->z[0]*x->z[i];
    }

  }
  else{
    for(i=0;i<np;++i){
      newx->x[i] = t->x[0]*x->x[i] + n->x[0]*x->y[i] + b->x[0]*x->z[i];
      newx->y[i] = t->y[0]*x->x[i] + n->y[0]*x->y[i] + b->y[0]*x->z[i];
      newx->z[i] = t->z[0]*x->x[i] + n->z[0]*x->y[i] + b->z[0]*x->z[i];
    }

    /* translate co-ordinates back to original position */

    dsadd(np, o->x[0], newx->x, 1, newx->x, 1);
    dsadd(np, o->y[0], newx->y, 1, newx->y, 1);
    dsadd(np, o->z[0], newx->z, 1, newx->z, 1);
  }
}


#define c1       ( 0.29690)
#define c2       (-0.12600)
#define c3       (-0.35160)
#define c4       ( 0.28430)
#define c5       (-0.10360)

/* naca profile -- usage: naca t x  returns points on naca 00 aerofoil of 
   thickness t at position x */
 
static double naca(double L, double x, double t){
  x = x/L;
  if(L==0.)
    return 0.;
  //  return 5.*t*L*(c1*sqrt(x)+ x*(c2 + x*(c3 + x*(c4 + c5*x))));
  return 5.*t*L*(c1*sqrt(x)+ x*(c2 + x*(c3 + x*(c4 + c5*x))));
}

#undef c1
#undef c2
#undef c3
#undef c4
#undef c5

void genNaca3d(Curve *curve, double *x, double *y, double *z, int q){
  register int i;
  Coord X,newX;
  double sg;

  X.x = x;
  X.y = y;
  X.z = z;
  
  newX.x = dvector(0, q-1);
  newX.y = dvector(0, q-1);
  newX.z = dvector(0, q-1);

  trans_coords(curve->info.nac3d.origin,
	       curve->info.nac3d.axis,
	       curve->info.nac3d.lead,
	       curve->info.nac3d.locz,
	       &X,&newX,q,1); // transform to aligned coord.s
	       
  /* move co-ordinate out to surface value */
  for(i = 0; i < q; ++i){
    sg = (newX.z[i]<0) ? -1. : 1.;
    if(newX.x[i] < 0. || newX.x[i] > curve->info.nac3d.length)
      fprintf(stderr, "X: %lf\n", newX.x[i]);
    
    newX.z[i]  = sg*naca(curve->info.nac3d.length,
			 newX.x[i],
			 curve->info.nac3d.thickness);
  }
  trans_coords(curve->info.nac3d.origin,
	       curve->info.nac3d.axis,
	       curve->info.nac3d.lead,
	       curve->info.nac3d.locz,
	       &newX,&X,q,-1); // transform to aligned coord.s
  
  free(newX.x);
  free(newX.y);
  free(newX.z);
}



void Quad_Face_JacProj(Bndry *B){

  Element *E = B->elmt;
  int      face = B->face;
  const    int qa = E->qa, qb = E->qb;
  Coord    S; 
  double   **da,**db,**dc,**dt;
  double   **D, *x, *y, *z, *xr, *xs, *yr, *ys, *zr, *zs;

  D = dmatrix(0,9,0,QGmax*QGmax-1);
  xr = D[0];  xs = D[1];
  yr = D[2];  ys = D[3];
  zr = D[4];  zs = D[5];
  
  S.x = D[0]; x = D[6];
  S.y = D[1]; y = D[7];
  S.z = D[2]; z = D[8];
  
  E->GetFaceCoord(face,&S);

  E->InterpToFace1(face,S.x,x);
  E->InterpToFace1(face,S.y,y);
  E->InterpToFace1(face,S.z,z);

  // db->da
  
  /* calculate derivatives */
  E->getD(&da,&dt,&db,&dt,&dc,&dt);

  /* calculate dx/dr */
  dgemm('T','N',qa,qb,qa,1.0,*da,qa,x,qa,0.0,xr,qa);
  
  /* calculate dx/ds */
  dgemm('N','N',qa,qb,qb,1.0,x,qa,*da,qb,0.0,xs,qa);
  
  /* calculate dy/dr */
  dgemm('T','N',qa,qb,qa,1.0,*da,qa,y,qa,0.0,yr,qa);
  
  /* calculate dy/ds */
  dgemm('N','N',qa,qb,qb,1.0,y,qa,*da,qb,0.0,ys,qa);
  
  /* calculate dz/dr */
  dgemm('T','N',qa,qb,qa,1.0,*da,qa,z,qa,0.0,zr,qa);
  
  /* calculate dz/ds */
  dgemm('N','N',qa,qb,qb,1.0,z,qa,*da,qb,0.0,zs,qa);

  /* x = yr*zs - zr*ys*/
  dvmul (qa*qb,yr,1,zs,1,x,1);
  dvvtvm(qa*qb,zr,1,ys,1,x,1,x,1);

  /* y = zr*xs - xr*zs*/
  dvmul (qa*qb,zr,1,xs,1,y,1);
  dvvtvm(qa*qb,xr,1,zs,1,y,1,y,1);

  /* z = xr*ys - xs*yr*/
  dvmul (qa*qb,xr,1,ys,1,z,1);
  dvvtvm(qa*qb,xs,1,yr,1,z,1,z,1);

  /* Surface Jacobea = sqrt(x^2 + y^2 + z^2) */
  
  dvmul (qa*qb,x,1,x,1,B->sjac.p,1);
  dvvtvp(qa*qb,y,1,y,1,B->sjac.p,1,B->sjac.p,1);
  dvvtvp(qa*qb,z,1,z,1,B->sjac.p,1,B->sjac.p,1);
  
  dvsqrt(qa*qb,B->sjac.p,1,B->sjac.p,1);
  
  free_dmatrix(D,0,0);
}
/* calculate the surface jacobian which is defined for as 2D surface in 
   a 3D space as:

   Surface Jac  = sqrt(Nx^2 + Ny^2 + Nz^2) 

   where [Nx,Ny,Nz] is the vector normal to the surface given by the
   cross product of the two tangent vectors in the r and s direction,
   i.e.

   Nx = y_r z_s - z_r y_s
   Ny = z_r x_s - x_r z_s
   Nz = x_r y_s - y_r x_s

   */

/* this function is defined with qa and qc to be compatible with the 
prism and pyramid triangular faces */
void Tri_Face_JacProj(Bndry *B){
  register int i;
  Element *E = B->elmt;
  int      face = B->face;
  const    int qa = E->qa, qb = E->qb, qc = E->qc;
  Coord    S; 
  double   **da,**db,**dc,**dt;
  double   **D, *x, *y, *z, *xr, *xs, *yr, *ys, *zr, *zs;
  Mode     *v = E->getbasis()->vert;

  D = dmatrix(0,9,0,QGmax*QGmax-1);
  xr = D[0];  xs = D[1];
  yr = D[2];  ys = D[3];
  zr = D[4];  zs = D[5];
  
  S.x = D[0]; x = D[6];
  S.y = D[1]; y = D[7];
  S.z = D[2]; z = D[8];
  
  E->GetFaceCoord(face,&S);

  E->InterpToFace1(face,S.x,x);
  E->InterpToFace1(face,S.y,y);
  E->InterpToFace1(face,S.z,z);

  /* calculate derivatives */
  E->getD(&da,&dt,&db,&dt,&dc,&dt);

  /* calculate dx/dr */
  dgemm('T','N',qa,qc,qa,1.0,*da,qa,x,qa,0.0,xr,qa);
  for(i = 0; i < qc; ++i)  dsmul(qa,1/v->c[i],xr+i*qa,1,xr+i*qa,1);
  
  /* calculate dx/ds */
  for(i = 0; i < qc; ++i) dvmul(qa,v[1].a,1,xr+i*qa,1,xs+i*qa,1);
  dgemm('N','N',qa,qc,qc,1.0,x,qa,*dc,qc,1.0,xs,qa);
  
  /* calculate dy/dr */
  dgemm('T','N',qa,qc,qa,1.0,*da,qa,y,qa,0.0,yr,qa);
  for(i = 0; i < qc; ++i)  dsmul(qa,1/v->c[i],yr+i*qa,1,yr+i*qa,1);
  
  /* calculate dy/ds */
  for(i = 0; i < qc; ++i) dvmul(qa,v[1].a,1,yr+i*qa,1,ys+i*qa,1);
  dgemm('N','N',qa,qc,qc,1.0,y,qa,*dc,qc,1.0,ys,qa);
  
  /* calculate dz/dr */
  dgemm('T','N',qa,qc,qa,1.0,*da,qa,z,qa,0.0,zr,qa);
  for(i = 0; i < qc; ++i)  dsmul(qa,1/v->c[i],zr+i*qa,1,zr+i*qa,1);
  
  /* calculate dz/ds */
  for(i = 0; i < qc; ++i) dvmul(qa,v[1].a,1,zr+i*qa,1,zs+i*qa,1);
  dgemm('N','N',qa,qc,qc,1.0,z,qa,*dc,qc,1.0,zs,qa);

  /* x = yr*zs - zr*ys*/
  dvmul (qa*qc,yr,1,zs,1,x,1);
  dvvtvm(qa*qc,zr,1,ys,1,x,1,x,1);

  /* y = zr*xs - xr*zs*/
  dvmul (qa*qc,zr,1,xs,1,y,1);
  dvvtvm(qa*qc,xr,1,zs,1,y,1,y,1);

  /* z = xr*ys - xs*yr*/
  dvmul (qa*qc,xr,1,ys,1,z,1);
  dvvtvm(qa*qc,xs,1,yr,1,z,1,z,1);

  /* Surface Jacobean = sqrt(x^2 + y^2 + z^2) */
  
  dvmul (qa*qc,x,1,x,1,B->sjac.p,1);
  dvvtvp(qa*qc,y,1,y,1,B->sjac.p,1,B->sjac.p,1);
  dvvtvp(qa*qc,z,1,z,1,B->sjac.p,1,B->sjac.p,1);
  
  dvsqrt(qa*qc,B->sjac.p,1,B->sjac.p,1);
  
  free_dmatrix(D,0,0);
}

#if 0
// only needed for explicit codes -- generate normals at Gauss points 'g'x'h'
void gen_face_normals(Element *E){
  int i;
  Bndry *Bc; // temporary bc 
  
  Bc = (Bndry*) calloc(1,sizeof(Bndry));
  Bc->elmt = E;
  for(i=0;i<E->Nfaces;++i){
    Bc->face = i;
    E->Surface_geofac(Bc);
    if(Bc->sjac.p){
      // need to fix for variable order
      E->face[i].sjac.p = dvector(0, E->lmax*E->lmax-1);
      E->InterpToGaussFace(0, Bc->sjac.p, E->lmax, E->lmax, E->face[i].sjac.p);
      free(Bc->sjac.p); 
      Bc->sjac.p = NULL;

      E->face[i].nx.p   = dvector(0, E->lmax*E->lmax-1);
      E->InterpToGaussFace(0,   Bc->nx.p, E->lmax, E->lmax,   E->face[i].nx.p);
      free(Bc->nx.p);
      Bc->nx.p = NULL;  

      E->face[i].ny.p   = dvector(0, E->lmax*E->lmax-1);
      E->InterpToGaussFace(0,   Bc->ny.p, E->lmax, E->lmax,   E->face[i].ny.p);
      free(Bc->ny.p);
      Bc->ny.p = NULL;  

      E->face[i].nz.p   = dvector(0, E->lmax*E->lmax-1);
      E->InterpToGaussFace(0,   Bc->nz.p, E->lmax, E->lmax,   E->face[i].nz.p);
      free(Bc->nz.p);
      Bc->nz.p = NULL;
    }
    else{
      E->face[i].sjac.d = Bc->sjac.d;
      E->face[i].nx.d   = Bc->nx.d;
      E->face[i].ny.d   = Bc->ny.d;
      E->face[i].nz.d   = Bc->nz.d;
    }
  }
}
#endif



void gen_ellipse(Element *E, Curve *cur, double *x, double *y){
  int i;
  double x0   = cur->info.ellipse.xo;
  double y0   = cur->info.ellipse.yo;
  double rmin = cur->info.ellipse.rmin;
  double rmaj = cur->info.ellipse.rmaj;

  Coord  X;
  X.x = dvector(0, QGmax-1);
  X.y = dvector(0, QGmax-1);
  
  double *xa = dvector(0, QGmax-1);
  double *ya = dvector(0, QGmax-1);

  E->straight_edge(&X, cur->face);
  E->InterpToFace1(cur->face, X.x, xa);
  E->InterpToFace1(cur->face, X.y, ya);

  double t0=0., t1=0.;

#if 1
  t0 = atan2((ya[0]-y0)*rmin, (xa[0]-x0)*rmaj);
  t1 = atan2((ya[E->qa-1]-y0)*rmin, (xa[E->qa-1]-x0)*rmaj);

  if(E->id == 162 || E->id == 387 || E->id == 445)
    fprintf(stderr, "id: %d t0: %lf t1: %lf\n", E->id+1,t0, t1);

#endif

  double t;

  double *z, *w;
  getzw(E->qa, &z, &w, 'a');

  if(t0 > 0 && t1 < 0){
    for(i=0;i<E->qa;++i){
      t = 0.5*(1-z[i])*t0 + 0.5*(1+z[i])*(t1+2.*M_PI);
      
      x[i] = x0 + rmaj*cos(t);
      y[i] = y0 + rmin*sin(t);
      if(E->id == 162 || E->id == 387 || E->id == 445)
	fprintf(stderr, "id: %d t: %lf x[%lf],y[%lf]\n", E->id+1, t, x[i] , y[i]);
    }
  }
  else if(t0 < 0 && t1 > 0){
    for(i=0;i<E->qa;++i){
      t = 0.5*(1-z[i])*(t0+2.*M_PI) + 0.5*(1+z[i])*t1;
      
      x[i] = x0 + rmaj*cos(t);
      y[i] = y0 + rmin*sin(t);
      if(E->id == 162 || E->id == 387 || E->id == 445)
	fprintf(stderr, "id: %d t: %lf x[%lf],y[%lf]\n", E->id+1, t, x[i] , y[i]);
      
    }
  }
  else{
    for(i=0;i<E->qa;++i){
      t = 0.5*(1-z[i])*t0 + 0.5*(1+z[i])*t1;
      
      x[i] = x0 + rmaj*cos(t);
      y[i] = y0 + rmin*sin(t);
      if(E->id == 162 || E->id == 387 || E->id == 445)
	fprintf(stderr, "id: %d t: %lf x[%lf],y[%lf]\n", E->id+1,t, x[i] , y[i]);
      
    }
  }
    
  
  free(xa);  free(ya);
  free(X.x); free(X.y);

  return;
}


void gen_sin(Element *E, Curve *cur, double *x, double *y){
  int i;
  double x0 = cur->info.sin.xo;
  double y0 = cur->info.sin.yo;
  double A  = cur->info.sin.amp;
  double lambda = cur->info.sin.wavelength;
  Coord  X;
  
  double *xa = dvector(0, QGmax-1);
  double *ya = dvector(0, QGmax-1);

  X.x = xa; X.y =ya;
  
  E->straight_edge(&X, cur->face);
  E->InterpToFace1(cur->face, X.x, x);
  for(i=0;i<E->qa;++i)
    y[i] = y0+A*sin(2.*M_PI*(x[i]-x0)/lambda);
  
  free(xa); free(ya);
  
  return;
}






#define c1       ( 0.29690)
#define c2       (-0.12600)
#define c3       (-0.35160)
#define c4       ( 0.28430)
#define c5       (-0.10360)

/* naca profile -- usage: naca t x  returns points on naca 00 aerofoil of 
   thickness t at position x 
 
   to compile
   
   cc -o naca naca.c -lm
   */
 
double Tri_naca(double L, double x, double t){
  x = x/L;
  if(L==0.)
    return 0.;
  //  return 5.*t*L*(c1*sqrt(x)+ x*(c2 + x*(c3 + x*(c4 + c5*x))));
  return 5.*t*L*(c1*sqrt(x)+ x*(c2 + x*(c3 + x*(c4 + c5*x))));
}

#undef c1
#undef c2
#undef c3
#undef c4
#undef c5

void Tri_genNaca(Element *E, Curve *curve, double *x, double *y){
  int   i;
  Coord X;
  X.x = dvector(0, QGmax-1);
  X.y = dvector(0, QGmax-1);

  E->GetFaceCoord(curve->face, &X);
  
  dcopy(E->qa, X.x, 1, x, 1);
  for(i=0;i<E->qa;++i){
    y[i] = Tri_naca(curve->info.nac2d.length,
		X.x[i]-curve->info.nac2d.xo,
		curve->info.nac2d.thickness);
    if(X.y[i] < curve->info.nac2d.yo)
      y[i] =  curve->info.nac2d.yo-y[i];
    else
      y[i] =  curve->info.nac2d.yo+y[i];
  }
  
  free(X.x);
  free(X.y);
}


#define distance(p1,p2) (sqrt((p2.x-p1.x)*(p2.x-p1.x)+(p2.y-p1.y)*(p2.y-p1.y)))

#define ITERNACA 1000
#define TOLNACA  1e-5
#define c1       ( 0.29690)
#define c2       (-0.12600)
#define c3       (-0.35160)
#define c4       ( 0.28430)
// #define c5       (-0.10150)
#define c5       (-0.10360)
/*all that follows is to set up a spline fitting routine from a data file*/

typedef struct geomf  {    /* Curve defined in a file */
  int           npts  ;    /* number of points        */
  int           pos   ;    /* last confirmed position */
  char         *name  ;    /* File/curve name         */
  double       *x, *y ;    /* coordinates             */
  double       *sx,*sy;    /* spline coefficients     */
  double       *arclen;    /* arclen along the curve  */
  struct geomf *next  ;    /* link to the next        */
} Geometry;

typedef struct vector {    /* A 2-D vector */
  double     x, y     ;    /* components   */
  double     length   ;    /* length       */
} Vector;

#define _MAX_NC         1024   /* Points describing a curved side   */
static int    closest    (Point p, Geometry *g);
static void   bracket    (double s[], double f[], Geometry *g, Point a,
			  Vector ap);
static Vector setVector  (Point p1, Point p2);

static double searchGeom (Point a, Point p, Geometry *g),
              brent      (double s[], Geometry *g, Point a, Vector ap,
			  double tol);

static Geometry *lookupGeom (char *name),
                *loadGeom   (char *name);

static Geometry *geomlist;

static Point setPoint (double x, double y)
{
  Point p;
  p.x = x; 
  p.y = y; 
  return p;
}


static Vector setVector (Point p1, Point p2)
{
  Vector v;

  v.x      = p2.x - p1.x;
  v.y      = p2.y - p1.y;
  v.length = sqrt (v.x*v.x + v.y*v.y);

  return v;
}

/* Compute the angle between the vector ap and the vector from a to
 * a point s on the curv.  Uses the small-angle approximation */
 
static double getAngle (double s, Geometry *g, Point a, Vector ap)
{
  Point  c;
  Vector ac;

  c  = setPoint (splint(g->npts, s, g->arclen, g->x, g->sx),
                 splint(g->npts, s, g->arclen, g->y, g->sy));
  ac = setVector(a, c);
			
  return 1. - ((ap.x * ac.x + ap.y * ac.y) / (ap.length * ac.length));
}

/* Search for the named Geometry */

static Geometry *lookupGeom (char *name)
{
  Geometry *g = geomlist;

  while (g) {
    if (strcmp(name, g->name) == 0) 
      return g;
    g = g->next;
  }

  return (Geometry *) NULL;
}

/* Load a geometry file */

static Geometry *loadGeom (char *name){
  const int verbose = option("verbose");
  Geometry *g   = (Geometry *) calloc (1, sizeof(Geometry));
  char      buf [BUFSIZ];
  double    tmp[_MAX_NC];
  Point     p1, p2, p3, p4;
  FILE     *fp;  
  register  int i;
  double  xscal = dparam("XSCALE");
  double  yscal = dparam("YSCALE");
  double  xmove = dparam("XMOVE");
  double  ymove = dparam("YMOVE");

  if (verbose > 1)
    printf ("Loading geometry file %s...", name);
  if ((fp = fopen(name, "r")) == (FILE *) NULL) {
    fprintf (stderr, "couldn't find the curved-side file %s", name);
    exit (-1);
  }

  while (fgets (buf, BUFSIZ, fp))    /* Read past the comments */
    if (*buf != '#') break;
  
  /* Allocate space for the coordinates */

  g -> x = (double*) calloc (_MAX_NC, sizeof(double));
  g -> y = (double*) calloc (_MAX_NC, sizeof(double));

  strcpy (g->name = (char *) malloc (strlen(name)+1), name);

  /* Read the coordinates.  The first line is already in *
   * the input buffer from the comment loop above.       */
  
  i = 0;
  while (i <= _MAX_NC && sscanf (buf,"%lf%lf", g->x + i, g->y + i) == 2) {
    i++;
    if (!fgets(buf, BUFSIZ, fp)) break;
  }
  g->npts = i;

  if(xmove)  dsadd(g->npts,xmove,g->x,1,g->x,1);
  if(ymove)  dsadd(g->npts,ymove,g->y,1,g->y,1);
  if(xscal)  dscal(g->npts,xscal,g->x,1);
  if(yscal)  dscal(g->npts,yscal,g->y,1);

  if (i < 2 ) error_msg (geometry file does not have enough points);

  if (i > _MAX_NC) error_msg (geometry file has too many points);

  if (verbose > 1) printf ("%d points", g->npts);

  /* Allocate memory for the other quantities */

  g->sx     = (double*) calloc (g->npts, sizeof(double));
  g->sy     = (double*) calloc (g->npts, sizeof(double));
  g->arclen = (double*) calloc (g->npts, sizeof(double));

  /* Compute spline information for the (x,y)-coordinates.  The vector "tmp"
     is a dummy independent variable for the function x(eta), y(eta).  */

  tmp[0] = 0.; 
  tmp[1] = 1.;
  dramp  (g->npts, tmp, tmp + 1, tmp, 1);
  spline (g->npts, 1.e30, 1.e30, tmp, g->x, g->sx);
  spline (g->npts, 1.e30, 1.e30, tmp, g->y, g->sy);

  /* Compute the arclength of the curve using 4 points per segment */

  for (i = 0; i < (*g).npts-1; i++) {
    p1 = setPoint (g->x[i], g->y[i] );
    p2 = setPoint (splint (g->npts, i+.25, tmp, g->x, g->sx),
		   splint (g->npts, i+.25, tmp, g->y, g->sy));
    p3 = setPoint (splint (g->npts, i+.75, tmp, g->x, g->sx),
		   splint (g->npts, i+.75, tmp, g->y, g->sy));
    p4 = setPoint (g->x[i+1], g->y[i+1]);

    g->arclen [i+1] = g->arclen[i] + distance (p1, p2) + distance (p2, p3) +
                                     distance (p3, p4);
  }

  /* Now that we have the arclength, compute x(s), y(s) */

  spline (g->npts, 1.e30, 1.e30, g->arclen, g->x, g->sx);
  spline (g->npts, 1.e30, 1.e30, g->arclen, g->y, g->sy);

  if (verbose > 1) 
    printf (", arclength  = %f\n", g->arclen[i]);


  /* add to the list of geometries */

  g ->next = geomlist;
  geomlist = g;

  fclose (fp);
  return g;
}

/* 
 * Find the point at which a line passing from the anchor point "a" 
 * through the search point "p" intersects the curve defined by "g".
 * Always searches from the last point found to the end of the curve.
 */

static double searchGeom (Point a, Point p, Geometry *g)
{
  Vector   ap;
  double   tol = dparam("TOLCURV"), s[3], f[3];
  register int ip;

  /* start the search at the closest point */

  ap   = setVector (a, p);
  s[0] = g -> arclen[ip = closest (p, g)];
  s[1] = g -> arclen[ip + 1];

  bracket (s, f, g, a, ap);
  if (fabs(f[1]) > tol) 
    brent (s, g, a, ap, tol); 

  return s[1];
}

int id_min(int n, double *d, int skip);
/* ---------------  Bracketing and Searching routines  --------------- */

static int closest (Point p, Geometry *g)
{
  const
  double  *x = g->x    + g->pos,
          *y = g->y    + g->pos;
  const    int n = g->npts - g->pos;
  double   len[_MAX_NC];
  register int i;
  
  for (i = 0; i < n; i++)
    len[i] = sqrt (pow(p.x - x[i],2.) + pow(p.y - y[i],2.));

  i = id_min (n, len, 1) + g->pos;
  i = min(i, g->npts-2);

  /* If we found the same position and it's not the very first *
   * one, start the search over at the beginning again.  The   *
   * test for i > 0 makes sure we only do the recursion once.  */

  if (i && i == g->pos) { g->pos = 0; i = closest (p, g); }

  return g->pos = i;
}

#define GOLD      1.618034
#define CGOLD     0.3819660
#define GLIMIT    100.
#define TINY      1.e-20
#define ZEPS      1.0e-10
#define ITMAX     100

#define SIGN(a,b)     ((b) > 0. ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d)  (a)=(b);(b)=(c);(c)=(d);
#define SHFT2(a,b,c)   (a)=(b);(b)=(c);

#define fa f[0]
#define fb f[1]
#define fc f[2]
#define xa s[0]
#define xb s[1]
#define xc s[2]

static void bracket (double s[], double f[], Geometry *g, Point a, Vector ap)
{
  double ulim, u, r, q, fu;

  fa = getAngle (xa, g, a, ap);
  fb = getAngle (xb, g, a, ap);

  if (fb > fa) { SHFT (u, xa, xb, u); SHFT (fu, fb, fa, fu); }

  xc = xb + GOLD*(xb - xa);
  fc = getAngle (xc, g, a, ap);

  while (fb > fc) {
    r = (xb - xa) * (fb - fc);
    q = (xb - xc) * (fb - fa);
    u =  xb - ((xb - xc) * q - (xb - xa) * r) /
              (2.*SIGN(max(fabs(q-r),TINY),q-r));
    ulim = xb * GLIMIT * (xc - xb);

    if ((xb - u)*(u - xc) > 0.) {      /* Parabolic u is bewteen b and c */
      fu = getAngle (u, g, a, ap);
      if (fu < fc) {                    /* Got a minimum between b and c */
	SHFT2 (xa,xb, u);
	SHFT2 (fa,fb,fu);
	return;
      } else if (fu > fb) {             /* Got a minimum between a and u */
	xc = u;
	fc = fu;
	return;
      }
      u  = xc + GOLD*(xc - xb);    /* Parabolic fit was no good. Use the */
      fu = getAngle (u, g, a, ap);             /* default magnification. */

    } else if ((xc-u)*(u-ulim) > 0.) {   /* Parabolic fit is bewteen c   */
      fu = getAngle (u, g, a, ap);                         /* and ulim   */
      if (fu < fc) {
	SHFT  (xb, xc, u, xc + GOLD*(xc - xb));
	SHFT  (fb, fc, fu, getAngle(u, g, a, ap));
      }
    } else if ((u-ulim)*(ulim-xc) >= 0.) {  /* Limit parabolic u to the  */
      u   = ulim;                           /* maximum allowed value     */
      fu  = getAngle (u, g, a, ap);
    } else {                                       /* Reject parabolic u */
      u   = xc + GOLD * (xc - xb);
      fu  = getAngle (u, g, a, ap);
    }
    SHFT  (xa, xb, xc, u);      /* Eliminate the oldest point & continue */
    SHFT  (fa, fb, fc, fu);
  }
  return;
}

/* Brent's algorithm for parabolic minimization */

static double brent (double s[], Geometry *g, Point ap, Vector app, double tol)
{
  int    iter;
  double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;

  a  = min (xa, xc);               /* a and b must be in decending order */
  b  = max (xa, xc);
  d  = 1.;
  x  = w  = v  = xb;
  fw = fv = fx = getAngle (x, g, ap, app);

  for (iter = 1; iter <= ITMAX; iter++) {    /* ....... Main Loop ...... */
    xm   = 0.5*(a+b);
    tol2 = 2.0*(tol1 = tol*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {             /* Completion test */
      xb = x;
      return fx;
    }
    if (fabs(e) > tol1) {             /* Construct a trial parabolic fit */
      r = (x-w) * (fx-fv);
      q = (x-v) * (fx-fw);
      p = (x-v) * q-(x-w) * r;
      q = (q-r) * 2.;
      if (q > 0.) p = -p;
      q = fabs(q);
      etemp=e;
      e = d;

      /* The following conditions determine the acceptability of the    */
      /* parabolic fit.  Following we take either the golden section    */
      /* step or the parabolic step.                                    */

      if (fabs(p) >= fabs(.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	d = CGOLD * (e = (x >= xm ? a-x : b-x));
      else {
	d = p / q;
	u = x + d;
	if (u-a < tol2 || b-u < tol2)
	  d = SIGN(tol1,xm-x);
      }
    } else
      d = CGOLD * (e = (x >= xm ? a-x : b-x));

    u  = (fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu = getAngle(u,g,ap,app);                     

    /* That was the one function evaluation per step.  Housekeeping... */

    if (fu <= fx) {
      if (u >= x) a = x; else b = x;
      SHFT(v ,w ,x ,u );
      SHFT(fv,fw,fx,fu)
    } else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
	v  = w;
	w  = u;
	fv = fw;
	fw = fu;
      } else if (fu <= fv || v == x || v == w) {
	v  = u;
	fv = fu;
      }
    }
  }                        /* .......... End of the Main Loop .......... */
  
  error_msg(too many iterations in brent());
  xb = x;
  return fx;
}

#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SIGN
#undef fa
#undef fb
#undef fc
#undef xa
#undef xb
#undef xc
