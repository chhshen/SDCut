
#include "mex.h"
#include "math.h"
#include "mexOliUtil.h"


void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{ 
  enum{ Ai , Bi};
  /*enum{ };*/
  oliCheckArgNumber(nrhs,2,nlhs,0);
  int m,n;
  float* pA = (float*)(oliCheckArg(prhs,Ai,&m,&n));
  float* pB = (float*)(oliCheckArg(prhs,Bi,m,n));
  
  for(int i=0;i<m*n;i++)
  {
    pA[i]+=pB[i];
  }

  
}


