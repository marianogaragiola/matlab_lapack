/******************************************************************
* Matlab mex function to solve the generalized eigenvalue problem
*      A*x = (lambda)*B*x
* with A a general complex matrix and B a definite positive complex 
* matrix using the function zggev of lapack.
* 
* For a correct use of the function first compile using mex compiler 
* from Matlab with the command
* 
* >> mex -v mzggev.c -lmwlapack
*******************************************************************/
#include <stdlib.h>
#include "mex.h"
#include "lapack.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *Ar, *Ai, *Br, *Bi, *A, *B; /* pointers to input & output matrices*/
	double *alpha, *beta;
	double *alphar, *alphai, *betar, *betai;
	double *vr, *vl;
	double *vrr, *vri;
	double *work, *rwork;
	size_t n, lda, ldb;
	size_t i, j;
	size_t ldvl, ldvr;
	size_t lwork, info;
	char *jobvl = "N", *jobvr = "V";

	/* Real and imaginary part of the first input */
	Ar = mxGetPr(prhs[0]);
	Ai = mxGetPi(prhs[0]);
	/* Real and imaginary part of the second input */
	Br = mxGetPr(prhs[1]);
	Bi = mxGetPi(prhs[1]);
	/* dimensions of the matrices */
	n = mxGetN(prhs[0]);
	lda = n; ldb = n;

	/* Pointer of the first matrix */
	A = (double *) mxMalloc(2*n*n*sizeof(double));
	/* Pointer of the second matrix */
	B = (double *) mxMalloc(2*n*n*sizeof(double));

	/* The matrices are passed to vector array */
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			A[2*(i*n+j)]   = Ar[i*n+j]; /* The real part in the even elements */
			A[2*(i*n+j)+1] = Ai[i*n+j]; /* The imaginary part in the odd elements */
			B[2*(i*n+j)]   = Br[i*n+j]; /* The real part in the even elements */
			B[2*(i*n+j)+1] = Bi[i*n+j]; /* The imaginary part in the odd elements */
		}
	}

	/* Pointer to alpha */
	alpha = (double *) mxMalloc(2*n*sizeof(double));

	/* Pointer to beta */
	beta = (double *) mxMalloc(2*n*sizeof(double));

	plhs[0] = mxCreateDoubleMatrix(n, 1, mxCOMPLEX);
	plhs[1] = mxCreateDoubleMatrix(n, 1, mxCOMPLEX);

	alphar = mxGetPr(plhs[0]);
	alphai = mxGetPi(plhs[0]);
	betar  = mxGetPr(plhs[1]);
	betai  = mxGetPi(plhs[1]);

	ldvl = 1; /* leading dimension for the left eigenvectors */
	ldvr = n; /* leading dimension for the right eigenvectors */

	vl = (double *) mxMalloc(2*ldvl*n*sizeof(double));
	vr = (double *) mxMalloc(2*ldvr*n*sizeof(double));

	plhs[2] = mxCreateDoubleMatrix(ldvr, n, mxCOMPLEX);
	vrr = mxGetPr(plhs[2]);
	vri = mxGetPi(plhs[2]);

	lwork = 2*n;
	work = (double *) mxMalloc(2*lwork*sizeof(double));

	rwork = (double *) mxMalloc(8*n*sizeof(double));

	/* Calling the lapack function */
	zggev(jobvl, jobvr, &n, A, &lda, B, &ldb, alpha,
	beta, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);

	for(i = 0; i < n; i++){
		alphar[i] = alpha[2*i];
		alphai[i] = alpha[2*i+1];
		betar[i]  = beta[2*i];
		betai[i]  = beta[2*i+i];
	}

	for(i = 0; i < ldvr; i++){
		for(j = 0; j < n; j++){
			vrr[i*n+j] = vr[2*(i*n+j)];
			vri[i*n+j] = vr[2*(i*n+j)+1];
		}
	}

	return;
}
