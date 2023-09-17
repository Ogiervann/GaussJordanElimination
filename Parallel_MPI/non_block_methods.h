#ifndef NON_BLOCK_METHODS_H
#define NON_BLOCK_METHODS_H
#include "matrinit.h"
int matrixInverse(int n, double *a, double *res, int* mpos, double epsA);
void returnToRegularMonke(int n, double* res, int* mpos);

double norml1(double* a, int n);
double normMatr(double* a, int n);
void matrixMult(int n, int s, int m, double *a, double* b, double* res);
void matrSubstr_eq(int n, int m, double* a, double* b);
void vectSubstr_eq(int n, double* a, double* b);
void switchLines(int i, int j, int n, double* a, double* b);
void switchColumns(int i, int j, int n, double* a, int* mpos);

void printMatrix_old(double* a, int l, int n, int r);


#endif
