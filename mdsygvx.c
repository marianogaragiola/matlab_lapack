/*********************************************************************
* Matlab mex function to solve the generalized eigenvalue problem 
*      A*x = (lambda)*B*x
* with A and B are real simetric matrix and B is also a positive 
* definite.
* 
* For a correct use of the function first compile using mex compiler 
* from matlab with the command
* 
* mex -v mdsygvx.m -lmwlapack
* **********************************************************************\
#include <stdlib.h>
#include "mex.h"
#include "lapack.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *A, *B, *Z, *w, *work; /* pointers to input & output matrices*/
	double *vl, *vu, abstol;
	size_t il, iu;
	size_t itype, n, m;
	size_t lda, ldb, ldz, lwork, info;
	size_t *iwork, *ifail;

	char *jobz = "V", *range = "A", *uplo = "U";

	itype = 1; /* solve the problem A*x = (lambda)*B*x */

	A = mxGetPr(prhs[0]); /* first input matrix */
	B = mxGetPr(prhs[1]);

	n = mxGetN(prhs[0]); /* number of rows or columns */
	lda = n;
	ldb = n;
	ldz = n;
	m = n;
	lwork = 8*n;

	vl = mxGetPr(prhs[2]); /* lower bound of the interval */
	vu = mxGetPr(prhs[3]); /* upper bound of the interval */

	il = 1;
	iu = n;

	abstol = 1e-12;

	plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
	w = mxGetPr(plhs[0]);

	plhs[1] = mxCreateDoubleMatrix(ldz, m, mxREAL);
	Z = mxGetPr(plhs[1]);

	plhs[2] = mxCreateDoubleMatrix(lwork, 1, mxREAL);
	work = mxGetPr(plhs[2]);

	iwork = (size_t *) mxMalloc(5*n*sizeof(size_t));

	ifail = (size_t *) mxMalloc(n*sizeof(size_t));

	/* Calling lapack function */
	dsygvx(&itype, jobz, range, uplo, &n, A, &lda, B, &ldb,
	vl, vu, &il, &iu, &abstol, &m, w, Z, &ldz, work, &lwork,
  iwork, ifail, &info);

}
