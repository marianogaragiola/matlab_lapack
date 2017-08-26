# matlab-lapack

This project exposes a set of _C_ functions to solve the eigenvalue and eigenvector problem with _Lapack_ and _Matlab_.

_C_ files are _mex_ functions which should be compiled in _Matlab_ using _mex_ compiler.

Example of `call dsyev`:

First compile _C_ function with

    mex -v mdsyev.c -lmwlapack

and the run `test_dsyev.m`
