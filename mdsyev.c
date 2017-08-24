/******************************************************************                       
* Matlab mex function to solve the eigenvalue problem
*      A*x = (lambda)*x
* with A a real simmetric matrix using the function dsyev of lapack.
*
* For a correct use of the function first compile using mex compiler
* from Matlab with the command
*
* >> mex -v mdsyev.c -lmwlapack
*******************************************************************/

#include "mex.h"
#include "lapack.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *A, *W, *WORK; /* pointers to input & output matrices*/
	size_t n, lda, lwork, info;      /* matrix dimensions */
	
	char *JOBZ = "N", *UPLO = "U";


	A = mxGetPr(prhs[0]); /* first input matrix */

	n = mxGetN(prhs[0]); /* number of rows or columns */
	lda = n;
	lwork = 3*n-1;

	plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
	W = mxGetPr(plhs[0]);

	plhs[1] = mxCreateDoubleMatrix(lwork, 1, mxREAL);
	WORK = mxGetPr(plhs[1]);

	/* Calling lapack function */
	dsyev(JOBZ, UPLO, &n, A, &lda, W, WORK, &lwork, &info);

}
