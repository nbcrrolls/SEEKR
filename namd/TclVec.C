/***************************************************************************
 *cr                                                                       
 *cr            (C) Copyright 1995-2008 The Board of Trustees of the           
 *cr                        University of Illinois                       
 *cr                         All Rights Reserved                        
 *cr                                                                   
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: TclVec.C,v $
 *      $Author: jim $        $Locker:  $             $State: Exp $
 *      $Revision: 1.2 $      $Date: 2008/09/17 16:19:54 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   A C-based implementation of some performance-critical Tcl callable
 *   routines in VMD.  The C implementation outperforms a raw Tcl version
 *   by a factor of three or so.  The performance advantage helps 
 *   significantly when doing analysis in VMD.
 ***************************************************************************/

#include <tcl.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include "eigeng.c"
#include "gauss.C"
// #include "TclCommands.h"
// #include "Matrix4.h"
// #include "utilities.h"

#define VMD_PI      3.14159265358979323846
#define VMD_TWOPI   (2.0 * VMD_PI)
#define DEGTORAD(a)     (a*VMD_PI/180.0)
#define RADTODEG(a)     (a*180.0/VMD_PI)

/***************** override some of the vector routines for speed ******/
/* These should be the exact C equivalent to the corresponding Tcl    */
/* vector commands */

// Function:  vecadd v1 v2 {v3 ...}
//  Returns: the sum of vectors; must all be the same length
//  The increase in speed from Tcl to C++ is 4561 / 255 == 18 fold
static int obj_vecadd(ClientData, Tcl_Interp *interp, int argc, 
		       Tcl_Obj * const objv[]) {
  if (argc < 3) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"vec1 vec2 ?vec3? ?vec4? ...");
    return TCL_ERROR;
  }
  int num;
  Tcl_Obj **data;
  if (Tcl_ListObjGetElements(interp, objv[1], &num, &data) != TCL_OK) {
    return TCL_ERROR;
  }
  double *sum = new double[num];
  int i;
  for (i=0; i<num; i++) {
    if (Tcl_GetDoubleFromObj(interp, data[i], sum+i) != TCL_OK) {
      delete [] sum;
      return TCL_ERROR;
    }
  }
  // do the sums on the rest
  int num2;
  for (int term=2; term < argc; term++) {
    if (Tcl_ListObjGetElements(interp, objv[term], &num2, &data) != TCL_OK) {
      delete [] sum;
      return TCL_ERROR;
    }
    if (num != num2) {
      Tcl_SetResult(interp, (char *) "vecadd: two vectors don't have the same size", TCL_STATIC);
      delete [] sum;
      return TCL_ERROR;
    }
    for (i=0; i<num; i++) {
      double df;
      if (Tcl_GetDoubleFromObj(interp, data[i], &df) != TCL_OK) {
	delete [] sum;
	return TCL_ERROR;
      }
      sum[i] += df;
    }
  }

  
  // and return the result
  Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
  for (i=0; i<num; i++) {
    Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewDoubleObj(sum[i]));
  }
  Tcl_SetObjResult(interp, tcl_result);
  delete [] sum;
  return TCL_OK;
}

// Function:  vecdot v1 v2
// added by Lane Votapka, Amaro lab UCSD 2014
//  Returns: the dot product of vectors; must all be the same length
//  The increase in speed from Tcl to C++ is unknown
static int obj_vecdot(ClientData, Tcl_Interp *interp, int argc, Tcl_Obj *const objv[]) {

  if (argc != 3) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"v1 v2");
    return TCL_ERROR;
  }

  int num1 = 0;
  Tcl_Obj **data1;
  if (Tcl_ListObjGetElements(interp, objv[1], &num1, &data1) != TCL_OK) {
    return TCL_ERROR;
  }
  
  int num2 = 0;
  Tcl_Obj **data2;
  if (Tcl_ListObjGetElements(interp, objv[2], &num2, &data2) != TCL_OK) {
    return TCL_ERROR;
  }
  
  if (num1 != num2) {
    Tcl_SetResult(interp, (char *) "vecdot: two vectors don't have the same size", TCL_STATIC);
    //delete [] sum;
    return TCL_ERROR;
  }
  
  double dot = 0.0;
  for (int i=0; i<num1; i++) {
    double tmp1; 
    double tmp2;
    if ((Tcl_GetDoubleFromObj(interp, data1[i], &tmp1) != TCL_OK) || (Tcl_GetDoubleFromObj(interp, data2[i], &tmp2) != TCL_OK)) {
      Tcl_SetResult(interp, (char *) "vecdot: non-numeric in vector", TCL_STATIC);
      return TCL_ERROR;
    } else {
      dot += tmp1*tmp2;
    }
  }
  
  Tcl_Obj *tcl_result = Tcl_GetObjResult(interp);
  Tcl_SetDoubleObj(tcl_result, dot);
  return TCL_OK; 
}

// Function:  veccross v1 v2
// added by Lane Votapka, Amaro lab UCSD 2014
//  Returns: the cross product of vectors; must all be the same length of 3
//  The increase in speed from Tcl to C++ is unknown
static int obj_veccross(ClientData, Tcl_Interp *interp, int argc, Tcl_Obj *const objv[]) {
  int i;
  if (argc != 3) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"v1 v2");
    return TCL_ERROR;
  }

  int num1 = 0;
  Tcl_Obj **data1;
  if (Tcl_ListObjGetElements(interp, objv[1], &num1, &data1) != TCL_OK) {
    return TCL_ERROR;
  }
  
  int num2 = 0;
  Tcl_Obj **data2;
  if (Tcl_ListObjGetElements(interp, objv[2], &num2, &data2) != TCL_OK) {
    return TCL_ERROR;
  }
  
  if ((num1 != 3 ) || ( num2 != 3)) {
    Tcl_SetResult(interp, (char *) "veccross: one of the two vectors don't have a size of 3", TCL_STATIC);
    //delete [] sum;
    return TCL_ERROR;
  }
  
  double *a = new double[num1];
  double *b = new double[num1];
  double cross[3];
  for (i=0; i<num1; i++) {
    
    if ((Tcl_GetDoubleFromObj(interp, data1[i], a+i) != TCL_OK) || (Tcl_GetDoubleFromObj(interp, data2[i], b+i) != TCL_OK)) {
      Tcl_SetResult(interp, (char *) "veccross: non-numeric in vector", TCL_STATIC);
      return TCL_ERROR;
    }
  }
  
  cross[0] = a[1]*b[2] - a[2]*b[1];
  cross[1] = a[2]*b[0] - a[0]*b[2];
  cross[2] = a[0]*b[1] - a[1]*b[0];
  
  // and return the result
  Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
  for (i=0; i<num1; i++) {
    Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewDoubleObj(cross[i]));
  }
  Tcl_SetObjResult(interp, tcl_result);
  delete [] a;
  delete [] b;
  return TCL_OK;
}


// Function:  vecsub  v1 v2
//  Returns:   v1 - v2


static int obj_vecsub(ClientData, Tcl_Interp *interp, int argc, Tcl_Obj *const objv[])
{
  if (argc != 3) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"?x? ?y?");
    return TCL_ERROR;
  }
  int num1=0, num2=0;
  Tcl_Obj **data1, **data2;
  if (Tcl_ListObjGetElements(interp, objv[1], &num1, &data1) != TCL_OK)
    return TCL_ERROR;
  if (Tcl_ListObjGetElements(interp, objv[2], &num2, &data2) != TCL_OK)
    return TCL_ERROR;

  if (num1 != num2) {
    Tcl_SetResult(interp, (char *)"vecsub: two vectors don't have the same size", TCL_STATIC);
    return TCL_ERROR;
  }

  Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
  for (int i=0; i<num1; i++) {
    double d1=0, d2=0;
    if (Tcl_GetDoubleFromObj(interp, data1[i], &d1) != TCL_OK) {
      Tcl_SetResult(interp, (char *)"vecsub: non-numeric in first argument", TCL_STATIC);
      return TCL_ERROR; 
    }
    if (Tcl_GetDoubleFromObj(interp, data2[i], &d2) != TCL_OK) {
      Tcl_SetResult(interp, (char *)"vecsub: non-numeric in second argument", TCL_STATIC);
      return TCL_ERROR; 
    }
    Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewDoubleObj(d1-d2));
  }
  Tcl_SetObjResult(interp, tcl_result);
  return TCL_OK;
}


// Function: vecscale
//  Returns: scalar * vector or vector * scalar
// speedup is 1228/225 = 5.5 fold
static int obj_vecscale(ClientData, Tcl_Interp *interp, int argc, 
		       Tcl_Obj * const objv[]) {
  if (argc != 3) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"?c? ?v?");
    return TCL_ERROR;
  }
    
  int num1 = 0, num2 = 0;
  Tcl_Obj **data1, **data2;
  if (Tcl_ListObjGetElements(interp, objv[1], &num1, &data1) != TCL_OK) {
    return TCL_ERROR;
  }
  if (Tcl_ListObjGetElements(interp, objv[2], &num2, &data2) != TCL_OK) {
    return TCL_ERROR;
  }
  if (num1 == 0 || num2 == 0) {
    Tcl_SetResult(interp, (char *) "vecscale: parameters must have data", TCL_STATIC);
    return TCL_ERROR;
  } else if (num1 != 1 && num2 != 1) {
    Tcl_SetResult(interp, (char *) "vecscale: one parameter must be a scalar value", TCL_STATIC);
    return TCL_ERROR;
  }
  
  int num = 0;
  Tcl_Obj *scalarobj, **vector;
  if (num1 == 1) {
    scalarobj = data1[0];
    vector = data2;
    num = num2;
  } else {
    scalarobj = data2[0];
    vector = data1;
    num = num1;
  }
 
  double scalar = 0.0;
  if (Tcl_GetDoubleFromObj(interp, scalarobj, &scalar) != TCL_OK)
    return TCL_ERROR;

  Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
  for (int i=0; i<num; i++) {
    double val = 0.0;
    if (Tcl_GetDoubleFromObj(interp, vector[i], &val) != TCL_OK) {
      Tcl_SetResult(interp, (char *) "vecscale: non-numeric in vector", TCL_STATIC);
      return TCL_ERROR;
    }
    val *= scalar;
    Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewDoubleObj(val));
  }
  Tcl_SetObjResult(interp, tcl_result);
  return TCL_OK;
}

/// Given a string with a matrix in it, return the matrix
// returns TCL_OK if good
// If bad, returns TCL_ERROR and sets the Tcl result to the error message
// The name of the function should be passed in 'fctn' so the error message
// can be constructed correctly



#if 0
int tcl_get_vector(const char *s, double *val, Tcl_Interp *interp)
{
  int num;
  const char **pos;
  if (Tcl_SplitList(interp, s, &num, &pos) != TCL_OK) {
    Tcl_SetResult(interp, (char *) "need three data elements for a vector", 
                  TCL_STATIC);
    return TCL_ERROR;
  }
  if (num != 3) {
    Tcl_SetResult(interp, (char *) "need three numbers for a vector", TCL_STATIC);
    return TCL_ERROR;
  }
  double a[3];
  if (Tcl_GetDouble(interp, pos[0], a+0) != TCL_OK ||
      Tcl_GetDouble(interp, pos[1], a+1) != TCL_OK ||
      Tcl_GetDouble(interp, pos[2], a+2) != TCL_OK) {
    ckfree((char *) pos); // free of tcl data
    return TCL_ERROR;
  }
  val[0] = (double) a[0];
  val[1] = (double) a[1];
  val[2] = (double) a[2];
  ckfree((char *) pos); // free of tcl data
  return TCL_OK;
}
#endif

int tcl_get_matrix(const char *fctn, Tcl_Interp *interp, 
			  Tcl_Obj *s, double *mat)
{ 
  int num_rows;
  Tcl_Obj **data_rows;
  if (Tcl_ListObjGetElements(interp, s, &num_rows, &data_rows) != TCL_OK) {
    char tmpstring[1024];
    sprintf(tmpstring, "%s: badly formed matrix", fctn);
    Tcl_SetResult(interp, tmpstring, TCL_VOLATILE);
    return TCL_ERROR;
  }
  if (num_rows != 4) {
    char tmpstring[1024];
    sprintf(tmpstring, "%s: need a 4x4 matrix", fctn);
    Tcl_SetResult(interp, tmpstring, TCL_VOLATILE);
    return TCL_ERROR;
  }
  int num_row[4];
  Tcl_Obj **data_row[4];
  if (Tcl_ListObjGetElements(interp, data_rows[0], num_row+0, data_row+0) != TCL_OK ||
      num_row[0] != 4 ||
      Tcl_ListObjGetElements(interp, data_rows[1], num_row+1, data_row+1) != TCL_OK ||
      num_row[1] != 4 ||
      Tcl_ListObjGetElements(interp, data_rows[2], num_row+2, data_row+2) != TCL_OK ||
      num_row[2] != 4 ||
      Tcl_ListObjGetElements(interp, data_rows[3], num_row+3, data_row+3) != TCL_OK ||
      num_row[3] != 4) {
    Tcl_AppendResult(interp, fctn, ": poorly formed matrix", NULL);
    return TCL_ERROR;
  }
  // now get the numbers
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      double tmp = 0.0;
      if (Tcl_GetDoubleFromObj(interp, data_row[i][j], &tmp) != TCL_OK) {
        char tmpstring[1024];
	sprintf(tmpstring, "%s: non-numeric in matrix", fctn);
        Tcl_SetResult(interp, tmpstring, TCL_VOLATILE);
        return TCL_ERROR;
      } else {
	mat[4*j+i] = (double)tmp;  // Matrix4 is transpose of Tcl's matrix
      }
    }
  }
  return TCL_OK;
}

// append the matrix into the Tcl result
void tcl_append_matrix(Tcl_Interp *interp, const double *mat) {
  Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
  for (int i=0; i<4; i++) {
    Tcl_Obj *m = Tcl_NewListObj(0, NULL);
    for (int j=0; j<4; j++) 
      Tcl_ListObjAppendElement(interp, m, Tcl_NewDoubleObj(mat[4*j+i]));
    Tcl_ListObjAppendElement(interp, tcl_result, m);
  }
  Tcl_SetObjResult(interp, tcl_result);
}

// append the 3x3 matrix into the Tcl result
void tcl_append_3matrix(Tcl_Interp *interp, const double *mat) {
  Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
  for (int i=0; i<3; i++) {
    Tcl_Obj *m = Tcl_NewListObj(0, NULL);
    for (int j=0; j<3; j++) 
      Tcl_ListObjAppendElement(interp, m, Tcl_NewDoubleObj(mat[3*j+i]));
    Tcl_ListObjAppendElement(interp, tcl_result, m);
  }
  Tcl_SetObjResult(interp, tcl_result);
}

int tcl_get_tallx3_matrix(const char *fctn, Tcl_Interp *interp, 
			  Tcl_Obj *s, double *mat, int *ptr_num_rows)
{ 
  int i, j;  
  //int num_rows;
  Tcl_Obj **data_rows;
  if (Tcl_ListObjGetElements(interp, s, ptr_num_rows, &data_rows) != TCL_OK) {
    char tmpstring[1024];
    sprintf(tmpstring, "%s: badly formed matrix", fctn);
    Tcl_SetResult(interp, tmpstring, TCL_VOLATILE);
    return TCL_ERROR;
  }
  int num_rows = ptr_num_rows[0];
  
  /*if (num_rows != 4) {
    char tmpstring[1024];
    sprintf(tmpstring, "%s: need a 4x4 matrix", fctn);
    Tcl_SetResult(interp, tmpstring, TCL_VOLATILE);
    return TCL_ERROR;
  }*/ // there can be as many rows as we want
  int num_row;
  //num_row = (int *)malloc(sizeof(int)*num_rows); // allocate the array to contain the row size for each
  Tcl_Obj ***data_row;
  data_row = (Tcl_Obj***)malloc(sizeof(Tcl_Obj**)*num_rows);
  //*mat = (double*)malloc(sizeof(double)*num_rows*3); // this matrix will contain the data itself
  for (i=0; i<num_rows; i++) {
    if (Tcl_ListObjGetElements(interp, data_rows[i], &num_row, data_row+i) != TCL_OK ||
      num_row != 3) {
      Tcl_AppendResult(interp, fctn, ": poorly formed matrix", NULL);
      return TCL_ERROR;
    } else {
      for (int j=0; j<3; j++) {
        double tmp;
        if (Tcl_GetDoubleFromObj(interp, data_row[i][j], &tmp) != TCL_OK) {
          char tmpstring[1024];
	  sprintf(tmpstring, "%f: non-numeric in matrix", tmp);
          Tcl_SetResult(interp, tmpstring, TCL_VOLATILE);
          return TCL_ERROR;
        } else {
          mat[3*i+j] = tmp;  // Matrix4 is transpose of Tcl's matrix
        }
      }
    }
  }
  free(data_row);
  return TCL_OK;
}

// Function: make_quat_ls_matrix coords1 coords2
// added by Lane Votapka, Amaro lab UCSD 2014
//  Returns: the 4x4 matrix for calculating quaternion least-squares
// speedup is unknown
static int obj_make_quat_ls_matrix(ClientData, Tcl_Interp *interp, int argc, 
		   Tcl_Obj * const objv[]) {
  // make there there are at least two values
  int i, j, k, el;
  if (argc < 2) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"m1");
    return TCL_ERROR;
  }
  // Get the first matrix
  int n = 0;
  if (Tcl_ListObjLength(interp, objv[1], &n) != TCL_OK)
    return TCL_ERROR;
  
  double *coords1 = new double[n*3]; // M is our skew-symmetric matrix
  double *coords2 = new double[n*3];
  int num_coords = 0;
  
  if (tcl_get_tallx3_matrix("make_quat_ls_matrix: ", interp, objv[1], coords1, &num_coords) != TCL_OK) {
    return TCL_ERROR;
  }
  if (tcl_get_tallx3_matrix("make_quat_ls_matrix: ", interp, objv[2], coords2, &num_coords) != TCL_OK) {
    return TCL_ERROR;
  }
  
  double M[16] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double Pi_trans[16];
  double Qi[16];
  double pi1, pi2, pi3, qi1, qi2, qi3;
  double tmp;
  char tmpstring[1024];
  
  for (i=0; i<n; i++) {
    pi1 = coords1[i*3]; pi2 = coords1[i*3 + 1]; pi3 = coords1[i*3 + 2];
    qi1 = coords2[i*3]; qi2 = coords2[i*3 + 1]; qi3 = coords2[i*3 + 2];
    Pi_trans[0] = 0.0;  Pi_trans[1] = pi1; Pi_trans[2] = pi2; Pi_trans[3] = pi3; 
    Pi_trans[4] = -pi1; Pi_trans[5] = 0.0; Pi_trans[6] =-pi3; Pi_trans[7] = pi2;
    Pi_trans[8] = -pi2; Pi_trans[9] = pi3; Pi_trans[10]= 0.0; Pi_trans[11]=-pi1;
    Pi_trans[12]= -pi3; Pi_trans[13]=-pi2;Pi_trans[14]= pi1; Pi_trans[15]= 0.0;
    //{0.0, pi1, pi2, pi3, 
    //-pi1, 0.0, -pi3, pi2, 
    //-pi2, pi3, 0.0, -pi1, 
    //-pi3, -pi2, pi1, 0.0 };
    Qi[0] = 0.0; Qi[1] =-qi1; Qi[2] =-qi2; Qi[3] =-qi3; 
    Qi[4] = qi1; Qi[5] = 0.0; Qi[6] =-qi3; Qi[7] = qi2;
    Qi[8] = qi2; Qi[9] = qi3; Qi[10]= 0.0; Qi[11]=-qi1;
    Qi[12]= qi3; Qi[13]=-qi2; Qi[14]= qi1; Qi[15]= 0.0;
    //{0.0, -qi1, -qi2, -qi3, 
    // qi1, 0.0, -qi3, qi2, 
    // qi2, qi3, 0.0, -qi1, 
    // qi3, -qi2, qi1, 0.0 };
    
    // matrix multiplcation and addition
    for (j=0; j<4; j++) { // first matrix's row
      for (k=0; k<4; k++) { // second matrix's column
        tmp = 0.0;
        for (el=0; el<4; el++) { // first matrix's column, second matrix's row
          tmp += Pi_trans[j*4 + el] * Qi[el*4 + k];
        }
        //sprintf(tmpstring, "%f: result in matrix[%d][%d]", tmp, k, j);
        //printf("%s\n", tmpstring);
        if (tmp != tmp) // if tmp is Nan 
	  tmp = 0.0;
        M[j*4 + k] += tmp;
      }   
    }
  }
  
  tcl_append_matrix(interp, M);
  //free(coords);
  return TCL_OK;
}



// speed up the matrix * vector routines -- DIFFERENT ERROR MESSAGES
// THAN THE TCL VERSION
// speedup is nearly 25 fold
static int obj_vectrans(ClientData, Tcl_Interp *interp, int argc, 
		  Tcl_Obj * const objv[])
{
  if (argc != 3) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"?matrix? ?vector?");
    return TCL_ERROR;
  }

  // get the matrix data
  double mat[16];
  if (tcl_get_matrix(
    Tcl_GetStringFromObj(objv[0],NULL), interp, objv[1], mat) != TCL_OK) {
    return TCL_ERROR;
  }
  
  // for the vector
  Tcl_Obj **vec;
  int vec_size;
  if (Tcl_ListObjGetElements(interp, objv[2], &vec_size, &vec) != TCL_OK)
    return TCL_ERROR;

  if (vec_size != 3 && vec_size != 4) {
    Tcl_SetResult(interp, (char *) "vectrans: vector must be of size 3 or 4",
                  TCL_STATIC);
    return TCL_ERROR;
  }
  double opoint[4];
  opoint[3] = 0;
  for (int i=0; i<vec_size; i++) {
    double tmp;
    if (Tcl_GetDoubleFromObj(interp, vec[i], &tmp) != TCL_OK) {
      Tcl_SetResult(interp, (char *) "vectrans: non-numeric in vector", TCL_STATIC);
      return TCL_ERROR;
    }
    opoint[i] = (double)tmp;
  }
  // vector data is in vec_data
  double npoint[4];
 
  npoint[0]=opoint[0]*mat[0]+opoint[1]*mat[4]+opoint[2]*mat[8]+opoint[3]*mat[12]
;
  npoint[1]=opoint[0]*mat[1]+opoint[1]*mat[5]+opoint[2]*mat[9]+opoint[3]*mat[13]
;
  npoint[2]=opoint[0]*mat[2]+opoint[1]*mat[6]+opoint[2]*mat[10]+opoint[3]*mat[14
];
  npoint[3]=opoint[0]*mat[3]+opoint[1]*mat[7]+opoint[2]*mat[11]+opoint[3]*mat[15
];
  // return it

  {
  Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
  for (int i=0; i<vec_size; i++) 
    Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewDoubleObj(npoint[i]));
  Tcl_SetObjResult(interp, tcl_result);
  }
  return TCL_OK;
}


// Function: transmult m1 m2 ... mn
//  Returns: the product of the matricies
// speedup is 136347 / 1316 = factor of 104
static int obj_transmult(ClientData, Tcl_Interp *interp, int argc, 
		   Tcl_Obj * const objv[]) {
  // make there there are at least two values
  if (argc < 3) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"mx my ?m1? ?m2? ...");
    return TCL_ERROR;
  }
  // Get the first matrix
  double mult[16];
  if (tcl_get_matrix("transmult: ", interp, objv[1], mult) != TCL_OK) {
    return TCL_ERROR;
  }
  int i = 2;
  double pre[16];
  while (i < argc) {
    if (tcl_get_matrix("transmult: ", interp, objv[i], pre) != TCL_OK) {
      return TCL_ERROR;
    }
    // premultiply mult by tmp
    double tmp[4];
    for (int k=0; k<4; k++) {
      tmp[0] = mult[k];
      tmp[1] = mult[4+k];
      tmp[2] = mult[8+k];
      tmp[3] = mult[12+k];
      for (int j=0; j<4; j++) {
        mult[4*j+k] = pre[4*j]*tmp[0] + pre[4*j+1]*tmp[1] +
          pre[4*j+2]*tmp[2] + pre[4*j+3]*tmp[3];
      }
    }
    i++;
  }
  tcl_append_matrix(interp, mult);
  return TCL_OK;
}

// Function: transadd m1 m2 ... mn 
// added by Lane Votapka, Amaro lab UCSD 2014
//  Returns: the sum of the matricies
// speedup is unknown
static int obj_transadd(ClientData, Tcl_Interp *interp, int argc, 
		   Tcl_Obj * const objv[]) {
  // make there there are at least two values
  if (argc < 3) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"mx my ?m1? ?m2? ...");
    return TCL_ERROR;
  }
  // Get the first matrix
  double sum_mat[16];
  if (tcl_get_matrix("transadd: ", interp, objv[1], sum_mat) != TCL_OK) {
    return TCL_ERROR;
  }
  int i = 2;
  double pre[16];
  while (i < argc) { // for each matrix
    if (tcl_get_matrix("transadd: ", interp, objv[i], pre) != TCL_OK) {
      return TCL_ERROR;
    }
    // add pre to sum_mat
    for (int k=0; k<16; k++) {
      sum_mat[k] += pre[k];
    }    
    /* multiply portion
    // premultiply mult by tmp
    double tmp[4];
    for (int k=0; k<4; k++) {
      tmp[0] = mult[k];
      tmp[1] = mult[4+k];
      tmp[2] = mult[8+k];
      tmp[3] = mult[12+k];
      for (int j=0; j<4; j++) {
        mult[4*j+k] = pre[4*j]*tmp[0] + pre[4*j+1]*tmp[1] +
          pre[4*j+2]*tmp[2] + pre[4*j+3]*tmp[3];
      }
    }
    */
    i++;
  }
  tcl_append_matrix(interp, sum_mat);
  return TCL_OK;
}

// Function: trans_principle_eig m1
// added by Lane Votapka, Amaro lab UCSD 2014
//  Returns: the principle eigenvalue of a matrix
// speedup is unknown
static int obj_trans_principal_eig(ClientData, Tcl_Interp *interp, int argc, 
		   Tcl_Obj * const objv[]) {
  int i, j;
  // make there there are at least two values
  if (argc != 2) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"mx");
    return TCL_ERROR;
  }
  // Get the first matrix
  double pre_mat[16];
  double mat[4][4];
  if (tcl_get_matrix("trans_principal_eig: ", interp, objv[1], pre_mat) != TCL_OK) {
    return TCL_ERROR;
  }
  for (i=0; i<4; i++) {
    for (j=0; j<4; j++) {
      mat[i][j] = pre_mat[i*4 + j]; // fill out a proper 4x4 matrix instead of a 16 length array
    }
  }
  //double array_evec[4]; // the eigenvector we are trying to find
  
  double *mem, *evalr, *evali, *evec; // memory, real parts of eigenvalues, imaginary parts of eigenvalues, eigenvectors
  mem = (double*)calloc(4 * 4 + 8, sizeof(double)); // allocate the memory for the eigenvals/vecs
  evec = mem; // the first 16 places reserved for eigenvectors
  evalr = evec + 16; // ... then 4 real components of eigenvals
  evali = evalr + 4; // ... then 4 imag components of eigenvals
  //n_eigeng(mat[0], 4, evalr, evali, evec); // find all evecs/evals
  /* // print all results
  for (i = 0; i < 5; ++i) {
    printf("%le + %le J\n", evalr[i], evali[i]);
  }
  printf("\n");
  for (i = 0; i <= 4; i++) {
    for (j = 0; j <= 4; j++)
      printf("%12.6e ", evec[i*5 + j]);
    printf("\n");
  }
  printf("\n");
  free(mem);*/
  /*for (i=0; i<4; i++) {
    array_evec[i][0] = evec[i*4];
  }*/
  // Find principal eigenvector by searching for the largest eigenvalue
  float largest_eval = -999999.9;
  int largest_eval_index = -1;
  for (i=0;i<4;i++) {
    if ( evalr[i] > largest_eval ) {
      largest_eval = evalr[i];
      largest_eval_index = i;
    }
  }
  if (largest_eval_index == -1) {
    //printf("C-CODE: largest eval: %f \n", largest_eval);
    //printf("C-CODE: largest eval index: %f \n", largest_eval_index);
    return TCL_ERROR; // then we have a problem
  }
  Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
  for (int i=0; i<4; i++) 
    Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewDoubleObj(evec[i*4 + largest_eval_index]));
  Tcl_SetObjResult(interp, tcl_result);
  
  return TCL_OK;
}



static int obj_transvec(ClientData, Tcl_Interp *interp, int argc, 
		   Tcl_Obj * const objv[]) {
  if (argc != 2) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"?vector?");
    return TCL_ERROR;
  }
  
  int num;
  Tcl_Obj **data;
  if (Tcl_ListObjGetElements(interp, objv[1], &num, &data) != TCL_OK) 
    return TCL_ERROR;
  if (num != 3) {
    Tcl_AppendResult(interp, "transvec: vector must have three elements",NULL);
    return TCL_ERROR;
  }
  double x,y,z;
  if (Tcl_GetDoubleFromObj(interp, data[0], &x) != TCL_OK ||
      Tcl_GetDoubleFromObj(interp, data[1], &y) != TCL_OK ||
      Tcl_GetDoubleFromObj(interp, data[2], &z) != TCL_OK) {
    Tcl_SetResult(interp, (char *)"transvec: non-numeric in vector", TCL_STATIC);
    return TCL_ERROR;
  }
  Matrix4 mat;
  mat.transvec((double) x,(double) y,(double) z);
  tcl_append_matrix(interp, mat.mat);
  return TCL_OK;
}

static int obj_transvecinv(ClientData, Tcl_Interp *interp, int argc, 
		   Tcl_Obj * const objv[]) {
  if (argc != 2) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"?vector?");
    return TCL_ERROR;
  }
  
  int num;
  Tcl_Obj **data;
  if (Tcl_ListObjGetElements(interp, objv[1], &num, &data) != TCL_OK) 
    return TCL_ERROR;
  if (num != 3) {
    Tcl_AppendResult(interp, "transvecinv: vector must have three elements",NULL);
    return TCL_ERROR;
  }
  double x,y,z;
  if (Tcl_GetDoubleFromObj(interp, data[0], &x) != TCL_OK ||
      Tcl_GetDoubleFromObj(interp, data[1], &y) != TCL_OK ||
      Tcl_GetDoubleFromObj(interp, data[2], &z) != TCL_OK) {
    Tcl_SetResult(interp, (char *)"transvecinv: non-numeric in vector", TCL_STATIC);
    return TCL_ERROR;
  }
  Matrix4 mat;
  mat.transvecinv((double) x,(double) y,(double) z);
  tcl_append_matrix(interp, mat.mat);
  return TCL_OK;
}

// Returns the transformation matrix needed to rotate by a certain
// angle around a given axis.
// Tcl syntax:
// transabout v amount [deg|rad|pi]
// The increase in speed from Tcl to C++ is 15 fold
static int obj_transabout(ClientData, Tcl_Interp *interp, int argc, 
		   Tcl_Obj * const objv[]) {
  if (argc != 3 && argc != 4) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"axis amount [deg|rad|pi]");
    return TCL_ERROR;
  }
  
  int num;
  Tcl_Obj **data;
  // get the axis
  if (Tcl_ListObjGetElements(interp, objv[1], &num, &data) != TCL_OK) 
    return TCL_ERROR;
  if (num != 3) {
    Tcl_AppendResult(interp, "transabout: vector must have three elements",NULL);
    return TCL_ERROR;
  }
  double x,y,z;
  if (Tcl_GetDoubleFromObj(interp, data[0], &x) != TCL_OK ||
      Tcl_GetDoubleFromObj(interp, data[1], &y) != TCL_OK ||
      Tcl_GetDoubleFromObj(interp, data[2], &z) != TCL_OK) {
    Tcl_SetResult(interp, (char *)"transabout: non-numeric in vector", TCL_STATIC);
    return TCL_ERROR;
  }

  // get the amount
  double amount;
  if (Tcl_GetDoubleFromObj(interp, objv[2], &amount) != TCL_OK) {
    Tcl_SetResult(interp, (char *)"transabout: non-numeric angle", TCL_STATIC);
    return TCL_ERROR;
  }

  // get units
  if (argc == 4) {
    if (!strcmp(Tcl_GetStringFromObj(objv[3], NULL), "deg")) {
      amount = DEGTORAD(amount);
    } else if (!strcmp(Tcl_GetStringFromObj(objv[3], NULL), "rad")) {
      // amount = amount; 
    } else if (!strcmp(Tcl_GetStringFromObj(objv[3], NULL), "pi")) {
      amount = amount*VMD_PI;
    } else {
      Tcl_AppendResult(interp, "transabout: unit must be deg|rad|pi",NULL);
      return TCL_ERROR;
    }
  } else {
    // If no unit was specified assume that we have degrees
    amount = DEGTORAD(amount);
  }

  double axis[3];
  axis[0] = (double) x;
  axis[1] = (double) y;
  axis[2] = (double) z;

  Matrix4 mat;
  mat.rotate_axis(axis, amount);
  tcl_append_matrix(interp, mat.mat);
  return TCL_OK;
}



static int obj_veclength(ClientData, Tcl_Interp *interp, int argc, Tcl_Obj *const objv[]) {

  if (argc != 2) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"?vector?");
    return TCL_ERROR;
  }

  int num;
  Tcl_Obj **data;
  if (Tcl_ListObjGetElements(interp, objv[1], &num, &data) != TCL_OK) 
    return TCL_ERROR;

  double length = 0.;
  for (int i=0; i<num; i++) {
    double tmp;
    if (Tcl_GetDoubleFromObj(interp, data[i], &tmp) != TCL_OK) {
      Tcl_SetResult(interp, (char *) "veclength: non-numeric in vector", TCL_STATIC);
      return TCL_ERROR;
    } else {
      length += tmp*tmp;
    }
  }

  length = sqrt(length);
  Tcl_Obj *tcl_result = Tcl_GetObjResult(interp);
  Tcl_SetDoubleObj(tcl_result, length);
  return TCL_OK; 
}


static double* obj_getdoublearray(Tcl_Interp *interp, Tcl_Obj *const objv[], int *len) {
  int num;

  Tcl_Obj **data;
  if (Tcl_ListObjGetElements(interp, objv[1], &num, &data) != TCL_OK)
    return NULL;
 
  double *list = (double*) malloc(num*sizeof(double));
  if (list == NULL)
    return NULL;

  for (int i=0; i<num; i++) {
    double tmp;
    if (Tcl_GetDoubleFromObj(interp, data[i], &tmp) != TCL_OK) {
      Tcl_SetResult(interp, (char *) "veclength: non-numeric in vector", TCL_STATIC);
      free(list);
      return NULL;
    }
    list[i] = tmp;
  }

  *len = num;

  return list;
}

// Function: getCOM_fast weights coords
// added by Lane Votapka, Amaro lab UCSD 2014
//  Returns: the center of a mass of a set of atoms
// speedup is unknown
static int obj_getCOM(ClientData, Tcl_Interp *interp, int argc, 
		   Tcl_Obj * const objv[]) {
  // make there there are at least two values
  // NOTE: THIS CODE IS BUGGY
  int i, j, k, el;
  if (argc < 3) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"weights coordinates");
    return TCL_ERROR;
  }
  // Get the first matrix
  int n;
  if (Tcl_ListObjLength(interp, objv[2], &n) != TCL_OK)
    return TCL_ERROR;
  double *coords1 = new double[n*3]; // M is our skew-symmetric matrix
  //double *coords2 = new double[n*3];
  int num_coords;
  if (tcl_get_tallx3_matrix("getCOM: ", interp, objv[2], coords1, &num_coords) != TCL_OK) {
    return TCL_ERROR;
  }
  int w_num;
  double *weights = obj_getdoublearray(interp, objv, &w_num);
  if (weights == NULL) 
    return TCL_ERROR;
  double totalweight = 0.0;
  for (i=0; i<n; i++) {
    totalweight += weights[i];
  }
  double inv_totalweights = 1.0 / totalweight;
  double COM[3] = {0.0,0.0,0.0}; // this is necessary, otherwise memory would contain random (potentially nonzero) values...
  //for (int i=0; i<3; i++)
  //  COM[i] = 0.0; 
  for (i=0; i<n; i++) {
    for (j=0; j<3; j++) {
      COM[j] += coords1[i*3 + j] * weights[i] * inv_totalweights;
    }
  }
  
  Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
  for (int i=0; i<3; i++) {
    //printf("COM[i]: %f\n", COM[i]);
    Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewDoubleObj(COM[i]));
  }
  Tcl_SetObjResult(interp, tcl_result);
  
  free(weights);
  free(coords1);
  return TCL_OK;
}

// Function: get_angular_accel torque weights coords
// added by Lane Votapka, Amaro lab UCSD 2014
//  Returns: the angular acceleration by calculating and inverting the inertia matrix and multiplying times the torque
// speedup is unknown
static int obj_get_angular_accel(ClientData, Tcl_Interp *interp, int argc, 
		   Tcl_Obj * const objv[]) {
  // make there there are at least two values
  int i, j, k, el;
  if (argc < 4) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"torque weights coordinates");
    return TCL_ERROR;
  }
  // Get the first matrix
  int n;
  if (Tcl_ListObjLength(interp, objv[4], &n) != TCL_OK)
    return TCL_ERROR;
  double *coords1 = new double[n*3]; // M is our skew-symmetric matrix
  int num_coords;
  // get the coordinates
  if (tcl_get_tallx3_matrix("get_angular_accel: ", interp, objv[4], coords1, &num_coords) != TCL_OK) {
    return TCL_ERROR;
  }
  if (n != num_coords) 
    return TCL_ERROR;
  
  // get the weights
  int w_num;
  double *weights = obj_getdoublearray(interp, &objv[1], &w_num);
  if (weights == NULL) 
    return TCL_ERROR;
  // get the torque
  int t_num;
  double *torque = obj_getdoublearray(interp, objv, &t_num);
  if (torque == NULL) 
    return TCL_ERROR;
  if (t_num != 3) 
    return TCL_ERROR;
  // get the center of mass
  int com_num;
  double *com = obj_getdoublearray(interp, &objv[2], &com_num);
  if (torque == NULL) 
    return TCL_ERROR;
  if (com_num != 3) 
    return TCL_ERROR;
  //double Imat[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
  // Construct the inertia matrix
  double **Imat;
  Imat = new double *[3];
  for (i=0; i<3; i++) {
    Imat[i] = new double[3];
    for (j=0; j<3; j++)
      Imat[i][j] = 0.0;
  }
  double x, y, z, m;
  for (i=0; i<n; i++) {
    x = coords1[i*3] - com[0]; y = coords1[i*3 + 1] - com[1]; z = coords1[i*3 + 2] - com[2]; m = weights[i];
    Imat[0][0] += m * (y*y + z*z); //Ixx
    Imat[0][1] -= m * x * y; // Ixy
    Imat[0][2] -= m * x * z; // Ixz
    //I[3] = I[1]
    Imat[1][1] += m * (x*x + z*z); //Iyy
    Imat[1][2] -= m * y * z; // Iyz
    //I[6] = I[2]
    //I[7] = I[5]
    Imat[2][2] += m * (x*x + y*y); //Iyy
  }
  Imat[1][0] = Imat[0][1];
  Imat[2][0] = Imat[0][2];
  Imat[2][1] = Imat[1][2];
  double ang_accel[3] = {0.0,0.0,0.0};
  // solve the inertia matrix for the torque to find the angular acceleration
  solve_matrix(3, Imat, ang_accel, torque);
  // return all results
  Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
  for (int i=0; i<3; i++) 
    Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewDoubleObj(ang_accel[i]));
  Tcl_SetObjResult(interp, tcl_result);
  free(weights);
  free(torque);
  free(com);
  return TCL_OK;
}

static int obj_vecsum(ClientData, Tcl_Interp *interp, int argc, Tcl_Obj *const objv[]) {
  if (argc != 2) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"?vector?");
    return TCL_ERROR;
  }

  int num;
  double *list = obj_getdoublearray(interp, objv, &num);
  if (list == NULL) 
    return TCL_ERROR;

  double sum = 0.;
  for (int i=0; i<num; i++) {
    sum += list[i];
  }
  free(list);

  Tcl_Obj *tcl_result = Tcl_GetObjResult(interp);
  Tcl_SetDoubleObj(tcl_result, sum);
  return TCL_OK;
}


static int obj_vecmean(ClientData, Tcl_Interp *interp, int argc, Tcl_Obj *const objv[]) {
  if (argc != 2) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"?vector?");
    return TCL_ERROR;
  }

  int num;
  double *list = obj_getdoublearray(interp, objv, &num);
  if (list == NULL) 
    return TCL_ERROR;

  double sum = 0.;
  for (int i=0; i<num; i++) {
    sum += list[i];
  }
  sum /= (double) num;
  free(list);

  Tcl_Obj *tcl_result = Tcl_GetObjResult(interp);
  Tcl_SetDoubleObj(tcl_result, sum);
  return TCL_OK;
}


static int obj_vecstddev(ClientData, Tcl_Interp *interp, int argc, Tcl_Obj *const objv[]) {
  if (argc != 2) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"?vector?");
    return TCL_ERROR;
  }

  int i, num;
  double* list = obj_getdoublearray(interp, objv, &num);
  if (list == NULL) 
    return TCL_ERROR;

  double mean = 0.;
  for (i=0; i<num; i++) {
    mean += list[i];
  }
  mean /= (double) num;

  double stddev = 0.;
  for (i=0; i<num; i++) {
    double tmp = list[i] - mean;
    stddev += tmp * tmp; 
  }
  stddev /= (double) num;
  stddev = sqrt(stddev);
  free(list);

  Tcl_Obj *tcl_result = Tcl_GetObjResult(interp);
  Tcl_SetDoubleObj(tcl_result, stddev);
  return TCL_OK;
}


int Vec_Init(Tcl_Interp *interp) {
  Tcl_CreateObjCommand(interp, "vecadd", obj_vecadd,
                    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "vecdot", obj_vecdot,
                    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "veccross", obj_veccross,
                    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "vecsub", obj_vecsub,
                    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "vecscale", obj_vecscale,
                    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "transmult", obj_transmult,
                    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "transadd", obj_transadd,
                    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "trans_principal_eig", obj_trans_principal_eig,
                    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "vectrans", obj_vectrans,
                    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "veclength", obj_veclength,
                    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "getCOM_fast", obj_getCOM,
                    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "vecmean", obj_vecmean,
                    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "vecstddev", obj_vecstddev,
                    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "vecsum", obj_vecsum,
                    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "transvec", obj_transvec,
                    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "transvecinv", obj_transvecinv,
                    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "transabout", obj_transabout,
                    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "make_quat_ls_matrix", obj_make_quat_ls_matrix,
                    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);  // not yet ready     
  Tcl_CreateObjCommand(interp, "get_angular_accel", obj_get_angular_accel,
                    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);  // not yet ready       
                    
  return TCL_OK;
}
 
