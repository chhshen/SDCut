// same as icd_ch2.c but with memory accesses instead of using RAM
// => slower but can handle huge matrices

#include "mex.h"
#include <math.h>

#include "icd_chi2_light.h"

 void mexFunction(int nlhs,
                  mxArray *plhs[],
                  int nrhs,
                  const mxArray *prhs[])
// args in :
//    - file names
//    - tolerance
//    - kernel parameter
//    - max size


n_pic = mxGetN(prhs[0]); /* number of pictures */
temp=mxGetPr(prhs[1]);
alpha=*temp;            /* kernel parameter */
temp=mxGetPr(prhs[2]);
tol=*temp;              /* approximation parameter */
temp=mxGetPr(prhs[3]);
nmax=*temp;             /* maximal rank */



