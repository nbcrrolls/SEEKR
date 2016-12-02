/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <stdlib.h> 
#ifndef _NO_MALLOC_H
#include <malloc.h>
#endif
#include <errno.h>
#include "TclCommands.h"
#include "Vector.h"
#include "NamdTypes.h"

#ifdef NAMD_TCL

#define USE_COMPAT_CONST
#include <tcl.h>

#include "Matrix4.C"
#include "TclVec.C"



// Get a 3-D vector from a TCL list
static int get_3D_vector(Tcl_Interp *interp, Tcl_Obj *const list, Vector &result)
{
  int num;
  Tcl_Obj **data;
  
  if (Tcl_ListObjGetElements(interp, list, &num, &data) != TCL_OK)
    return 0;
  if (Tcl_GetDoubleFromObj(interp,data[0],&(result.x)) != TCL_OK)
    return 0;
  if (Tcl_GetDoubleFromObj(interp,data[1],&(result.y)) != TCL_OK)
    return 0;
  if (Tcl_GetDoubleFromObj(interp,data[2],&(result.z)) != TCL_OK)
    return 0;

  return 1;
}


// Append a 3-D vector to the result string
static Tcl_Obj* obj_3D_vector(const Vector &v)
{
  Tcl_Obj* doublev[3];
  
  doublev[0] = Tcl_NewDoubleObj(v.x);
  doublev[1] = Tcl_NewDoubleObj(v.y);
  doublev[2] = Tcl_NewDoubleObj(v.z);
  
  return Tcl_NewListObj(3,doublev);
}


// Function: getbond coor1 coor2
//  Returns: the length of the bond formed by the two atoms (i.e., the distance between them)
int proc_getbond(ClientData, Tcl_Interp *interp, int argc, Tcl_Obj *const argv[])
{
  Vector r1, r2;
  
  if (argc != 3)
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[1],r1))
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[2],r2))
    return TCL_ERROR;
  Tcl_SetObjResult(interp, Tcl_NewDoubleObj((r2-r1).length()));
  return TCL_OK;
}


// Function: getangle coor1 coor2 coor3
//  Returns: the angle formed by the three atoms
int proc_getangle(ClientData, Tcl_Interp *interp, int argc, Tcl_Obj *const argv[])
{
  Vector r1, r2, r3, r12, r32;
  
  if (argc != 4)
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[1],r1))
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[2],r2))
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[3],r3))
    return TCL_ERROR;
  r12 = r1 - r2;
  r32 = r3 - r2;
  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(
		acos((r12*r32)/(r12.length()*r32.length()))*180/PI  ));
  return TCL_OK;
}


// Function: getdihedral coor1 coor2 coor3 coor4
//  Returns: the dihedral formed by the four atoms
int proc_getdihedral(ClientData, Tcl_Interp *interp, int argc, Tcl_Obj *const argv[])
{
  BigReal rA, rB, rC;
  Vector r1, r2, r3, r4, r12, r23, r34, A, B, C;
  
  if (argc != 5)
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[1],r1))
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[2],r2))
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[3],r3))
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[4],r4))
    return TCL_ERROR;
  r12 = r1 - r2;
  r23 = r2 - r3;
  r34 = r3 - r4;
  A = cross(r12,r23);
  B = cross(r23,r34);
  C = cross(r23,A);
  rA = A.length();
  rB = B.length();
  rC = C.length();
  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(
  		-atan2((C*B)/(rC*rB),(A*B)/(rA*rB))*180/PI  ));
  return TCL_OK;
}


// Function: anglegrad coor1 coor2 coor3
//  Returns: a list of gradients for each atom
// The code was basically copied from ComputeAngles.C
int proc_anglegrad(ClientData, Tcl_Interp *interp, int argc, Tcl_Obj *const argv[])
{
  Vector r1, r2, r3, r12, r32;
  
  if (argc != 4)
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[1],r1))
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[2],r2))
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[3],r3))
    return TCL_ERROR;
  
  r12 = r1 - r2;
  BigReal d12 = r12.length();
  r32 = r3 - r2;
  BigReal d32 = r32.length();
  
  BigReal cos_theta = (r12*r32)/(d12*d32);
  
  //  Normalize vector r12 and r32
  BigReal d12inv = 1. / d12;
  BigReal d32inv = 1. / d32;
  
  BigReal sin_theta = sqrt(1.0 - cos_theta*cos_theta);
  BigReal diff = 1/sin_theta;
  BigReal c1 = diff * d12inv;
  BigReal c2 = diff * d32inv;
  
  //  Calculate the actual forces
  Force force1 = c1*(r12*(d12inv*cos_theta) - r32*d32inv);
  Force force2 = force1;
  Force force3 = c2*(r32*(d32inv*cos_theta) - r12*d12inv);
  force2 += force3;  force2 *= -1;
  
  Tcl_Obj* forcev[3];

  forcev[0] = obj_3D_vector(force1);
  forcev[1] = obj_3D_vector(force2);
  forcev[2] = obj_3D_vector(force3);

  Tcl_SetObjResult(interp, Tcl_NewListObj(3, forcev));
  
  return TCL_OK;
}


// Function: dihedralgrad coor1 coor2 coor3 coor4
//  Returns: a list of gradients for each atom
// The code was basically copied from ComputeDihedrals.C
int proc_dihedralgrad(ClientData, Tcl_Interp *interp, int argc, Tcl_Obj *const argv[])
{
  BigReal K1;
  Vector r1, r2, r3, r4, r12, r23, r34;
  
  if (argc != 5)
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[1],r1))
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[2],r2))
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[3],r3))
    return TCL_ERROR;
  if (!get_3D_vector(interp,argv[4],r4))
    return TCL_ERROR;
  
  r12 = r1 - r2;
  r23 = r2 - r3;
  r34 = r3 - r4;
  
  //  Calculate the cross products and distances
  Vector A = cross(r12,r23);
  BigReal rA = A.length();
  Vector B = cross(r23,r34);
  BigReal rB = B.length();
  Vector C = cross(r23,A);
  BigReal rC = C.length();
  
  //  Calculate the sin and cos
  BigReal cos_phi = (A*B)/(rA*rB);
  BigReal sin_phi = (C*B)/(rC*rB);
  
  Force f1,f2,f3;
  
  //  Normalize B
  rB = 1.0/rB;
  B *= rB;
  
  //  We first need to figure out whether the
  //  sin or cos form will be more stable.  For this,
  //  just look at the value of phi
  if (fabs(sin_phi) > 0.1)
  {
    //  use the sin version to avoid 1/cos terms
    
    //  Normalize A
    rA = 1.0/rA;
    A *= rA;
    Vector dcosdA = rA*(cos_phi*A-B);
    Vector dcosdB = rB*(cos_phi*B-A);
    
    K1 = -1/sin_phi;
    
    f1.x = K1*(r23.y*dcosdA.z - r23.z*dcosdA.y);
    f1.y = K1*(r23.z*dcosdA.x - r23.x*dcosdA.z);
    f1.z = K1*(r23.x*dcosdA.y - r23.y*dcosdA.x);
    
    f3.x = K1*(r23.z*dcosdB.y - r23.y*dcosdB.z);
    f3.y = K1*(r23.x*dcosdB.z - r23.z*dcosdB.x);
    f3.z = K1*(r23.y*dcosdB.x - r23.x*dcosdB.y);
    
    f2.x = K1*(r12.z*dcosdA.y - r12.y*dcosdA.z
             + r34.y*dcosdB.z - r34.z*dcosdB.y);
    f2.y = K1*(r12.x*dcosdA.z - r12.z*dcosdA.x
             + r34.z*dcosdB.x - r34.x*dcosdB.z);
    f2.z = K1*(r12.y*dcosdA.x - r12.x*dcosdA.y
             + r34.x*dcosdB.y - r34.y*dcosdB.x);
  }
  else
  {
    //  This angle is closer to 0 or 180 than it is to
    //  90, so use the cos version to avoid 1/sin terms
    
    //  Normalize C
    rC = 1.0/rC;
    C *= rC;
    Vector dsindC = rC*(sin_phi*C-B);
    Vector dsindB = rB*(sin_phi*B-C);
    
    K1 = 1/cos_phi;
    
    f1.x = K1*((r23.y*r23.y + r23.z*r23.z)*dsindC.x
              - r23.x*r23.y*dsindC.y
              - r23.x*r23.z*dsindC.z);
    f1.y = K1*((r23.z*r23.z + r23.x*r23.x)*dsindC.y
              - r23.y*r23.z*dsindC.z
              - r23.y*r23.x*dsindC.x);
    f1.z = K1*((r23.x*r23.x + r23.y*r23.y)*dsindC.z
              - r23.z*r23.x*dsindC.x
              - r23.z*r23.y*dsindC.y);
    
    f3 = cross(K1,dsindB,r23);
    
    f2.x = K1*(-(r23.y*r12.y + r23.z*r12.z)*dsindC.x
           +(2.0*r23.x*r12.y - r12.x*r23.y)*dsindC.y
           +(2.0*r23.x*r12.z - r12.x*r23.z)*dsindC.z
           +dsindB.z*r34.y - dsindB.y*r34.z);
    f2.y = K1*(-(r23.z*r12.z + r23.x*r12.x)*dsindC.y
           +(2.0*r23.y*r12.z - r12.y*r23.z)*dsindC.z
           +(2.0*r23.y*r12.x - r12.y*r23.x)*dsindC.x
           +dsindB.x*r34.z - dsindB.z*r34.x);
    f2.z = K1*(-(r23.x*r12.x + r23.y*r12.y)*dsindC.z
           +(2.0*r23.z*r12.x - r12.z*r23.x)*dsindC.x
           +(2.0*r23.z*r12.y - r12.z*r23.y)*dsindC.y
           +dsindB.y*r34.x - dsindB.x*r34.y);
  }
  
  Tcl_Obj* forcev[4];

  forcev[0] = obj_3D_vector(f1);
  forcev[1] = obj_3D_vector(f2-f1);
  forcev[2] = obj_3D_vector(f3-f2);
  forcev[3] = obj_3D_vector(-f3);

  Tcl_SetObjResult(interp, Tcl_NewListObj(4, forcev));
  
  return TCL_OK;
}

int tcl_vector_math_init(Tcl_Interp *interp) {

  // first import from TclVec.C stolen from VMD
  Vec_Init(interp);

  Tcl_CreateObjCommand(interp, "getbond", proc_getbond,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "getangle", proc_getangle,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "getdihedral", proc_getdihedral,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "anglegrad", proc_anglegrad,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "dihedralgrad", proc_dihedralgrad,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);

  return TCL_OK;
}


#endif
