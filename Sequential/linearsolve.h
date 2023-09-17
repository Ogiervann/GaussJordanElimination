#ifndef LINEARSOLVE_H
#define LINEARSOLVE_H
#include <cmath>
#include "lin_operations.h"
#include "matrinit.h"


void returnToRegularMonke(int n, double* res, int* mpos);



int matrixInverse(int n, double *a, double *res, int* mpos, double epsA);

int findMaxBlock(int i0, int n, int m, int k, int l, double* a, double* block1,
  double* block2, int* mpos, double epsA);

void returnToMonke(int n, int m,  double* x, double* b, int *pos);

int linearSolve(int n, int m, double* a, double* b, double* x, int* pos,
  double* block1, double* block2, double *block3, int* mpos, double epsA);

double findr1(int n, int m, double* a, double* b, double* x, double* block1,
  double* block2, double* block3);

double findr2(int n, double* x);
#endif
