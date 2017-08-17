
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

	/* Llamo a la funcion dsyev */
	dsyev(JOBZ, UPLO, &n, A, &lda, W, WORK, &lwork, &info);

}
