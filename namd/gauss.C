#include<stdlib.h>
#include<stdio.h>
#include<math.h>

int solve_matrix(int n, double **a, double *x, double *b)
{
  double *temp, btemp,app,sum,mult;
  int i,j,k,p;
  temp = new double[n];
  
  for(i=0;i<n;i++){
    app = a[i][i];
    //initialization of p
    p = i;
    //find largest no of the columns and row no. of largest no.
    for(k = i+1; k < n; k++)
      if(fabs(app) < fabs(a[k][i])){
        app = a[k][i] ;
        p = k;
      }
      //swaping the elements of diagonal row and row containing largest no
      for(j = 0; j < n; j++) {
        temp[j] = a[p][j];
        a[p][j] = a[i][j];
        a[i][j] = temp[j];
      }
      btemp = b[p];
      b[p] = b[i];
      b[i] = btemp;
            
      //calculating triangular matrix
      for(j=i+1;j<n;j++){
        mult = a[j][i]/a[i][i];
        for(k=0;k<n;k++)
          a[j][k] -= mult*a[i][k];
          b[j] -= mult*b[i];
        }
      }
      //for calculating value of x via backward substitution method
      for(i=n-1;i>=0;i--){
        sum = 0;
        for(j=i+1;j<n;j++)
          sum += a[i][j]*temp[j];
        temp[i] = (b[i]-sum)/a[i][i];
      }
      for(i=0;i<n;i++) { // fill out the values of x
        x[i] = temp[i];
      }
      
    return 0;
}

/*int main(int argc, char *argv[]) {
  int i;
  double **a;
  a = new double *[3];
  for (i=0; i<3; i++) {
    a[i] = new double[3];
  }
  a[0][0] = 2.0; a[0][1] = 1.0; a[0][2] = 3.0; 
  a[1][0] = 2.0; a[1][1] = 6.0; a[1][2] = 8.0; 
  a[2][0] = 6.0; a[2][1] = 8.0; a[2][2] = 18.0; 
  double *b = new double[3];
  b[0] = 1.0; b[1] = 3.0; b[2] = 5.0;
  double *x = new double[3];
  solve_matrix(3,a,x,b);
  return 0;
}*/
