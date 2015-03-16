#include "mex.h"
#include <math.h>

void mexFunction(int nlhs,//Number of expected mxArrays (Left Hand Side)
                  mxArray *plhs[],//Array of pointers to expected outputs
                  int nrhs,//Number of inputs (Right Hand Side)
                  const mxArray *prhs[])/*Array of pointers to input data. 
		  The input data is read-only and should not be altered 
		  by your mexFunction .*/
{

/* 2 inputs :
		- a vector x
		- chi2 parameter (exp(-lambda*chi2(x))
*/

/* dimensions of x : */
int d = mxGetM(prhs[0]);
int n = mxGetN(prhs[0]);

/* get x : */
double *x;
x = mxGetPr(prhs[0]);



/* get chi2 parameters : */
double lambda;
double *tmp;
tmp= mxGetPr(prhs[1]);
lambda=*tmp;

/* 1 output : 
		-  the kernel K  (x,x)
*/
plhs[0]=mxCreateDoubleMatrix(n,n,mxREAL);
/* 	matlab stocke en colonne donc si on fait :
 	mxCreateDoubleMatrix(m,n,0) => en C on cree un Array de taille m*n
	qui represente la mise a plat de la matrice suivant les colonne i.e
	comme il y a n colonne de taille m le C-Array doit etre considere 
	comme n block de taille m
*/
double *K;
K=mxGetPr(plhs[0]);

int i,j,dd;
double k;
for (i=0;i<n;i++)
{
	K[i+n*i]=1;
	K[n*i+i]=1;
	for(j=i+1;j<n;j++)
	{
		k=0.;
		for(dd=0;dd<d;dd++)
		{
			k+=(x[dd+d*i]-x[dd+d*j])*(x[dd+d*i]-x[dd+d*j])/
				((x[dd+d*i]+x[dd+d*j]!=0)?x[dd+d*i]+x[dd+d*j]:1);
		}
		K[j+n*i]=exp(-lambda*k);
		K[i+n*j]=K[j+n*i];
	}
}
	
}
		  
		  
