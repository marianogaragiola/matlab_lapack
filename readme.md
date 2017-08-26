Set of C functions to solve the eigenvalue and eigenvector problem with Lapack and Matlab.

The C files are mex functions which must to be compiled in Matlab using mex compiler.

Example of call dsyev:

First compile C function with

mex -v mdsyev.c -lmwlapack

and the run test_dsyev.m
