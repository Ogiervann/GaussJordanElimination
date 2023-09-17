#ifndef LINEARSOLVE_H
#define LINEARSOLVE_H
#include <cmath>
#include "lin_operations.h"
#include "matrinit.h"
#include "basic_funcs.h"
#include "used_classes.h"


void returnToRegularMonke(int n, double* res, int* mpos);

int matrixInverse(int n, double *a, double *res, int* mpos, double epsA);



void findMaxBlock_byThread(int i0, int n, int m, int k, int l, double* a, double* block1,
  double* block2, int* mpos, maxBlockResults *mbres, double epsA, int t, int p);

int findMaxBlock(int i0, int n, int m, int k, int l, double* a, double* block1,
  double* block2, int* mpos, maxBlockResults *mbres, double epsA, int t, int p);

void returnToMonke(int n, int m,  double* x, double* b, int *pos, double* block, int t, int p);

int linearSolve(int n, int m, double* a, double* b, double* x, int* pos,
   double* block1, double* block2, double *block3, int* mpos, maxBlockResults *mbres, double epsA, int t, int p);


double findr1(int n, int m, double* a, double* b, double* x, double* block1,
  double* block2, double* block3, int t, int p);

double findr2(int n, double* x);

#endif
