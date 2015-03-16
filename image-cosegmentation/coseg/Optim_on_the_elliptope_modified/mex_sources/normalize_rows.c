#include "mex.h"
#include "math.h"

void wxyz(double *x, double *z, mwSize m, mwSize n)
{
 mwSize i,j,count1,count2=0;
 
    for (i=0; i<m; i++) {
         double d = 0;
         for (j=0; j<n; j++) {
            count1=m*j+i;
            d += *(x+count1) * *(x+count1);            
        }
         d = 1/sqrt(d);
         for (j=0; j<n; j++) {
             count2=m*j+i;
            *(z+count2)= d * *(x+count2);
        }
    }
}



/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *x,*z;
  mwSize mrows,ncols;
  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
  
  /*  create a pointer to the input matrix y */
  x = mxGetPr(prhs[0]);
  
  /*  get the dimensions of the matrix input y */
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  
  /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix(mrows,ncols, mxREAL);
  
  /*  create a C pointer to a copy of the output matrix */
  z = mxGetPr(plhs[0]);
  
  /*  call the C subroutine */
  wxyz(x,z,mrows,ncols);
  
}
