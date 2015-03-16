#include "mex.h"
#include "math.h"
#include "mexOliUtil.h"


void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{ 
  enum{ Ai , X};
  /*enum{ };*/
  /* return A+XX^T */
  oliCheckArgNumber(nrhs,2,nlhs,0);
  int m,n;
  float* pA = (float*)(oliCheckArg(prhs,Ai,&m,&n));
  float* pX = (float*)(oliCheckArg(prhs,X,m,1));
  
  for(int k=0;k<n;k++)
      for(int j=0;j<m;j++)
        pA[j+m*k]-=pX[k]*pX[j];

  
}
