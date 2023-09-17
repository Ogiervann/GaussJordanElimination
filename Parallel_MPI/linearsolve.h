#ifndef LINEARSOLVE_H
#define LINEARSOLVE_H
#include <cmath>
#include <cstring>
#include "lin_operations.h"
#include "matrinit.h"
#include "mpi_methods.h"
#include "basic_funcs.h"
#include "used_classes.h"
#include "non_block_methods.h"



Pos findMaxBlock(int q, double* a, double* block1, double* block2, int* mpos,
   double epsA, int n, int m, int k, int p, MPI_Comm comm, MPI_Op op_max_Pos,
   MPI_Datatype MY_MPI_POS);

void returnToMonke(double* x, double* b, int *pos, double* block, int n, int m,
  int k, int p, MPI_Comm comm);

int linearSolve(double* a, double* b, double* x, double* buff, int* pos,
  double* block1, double* block2, double *block3, int* mpos,
  double epsA, int n, int m, int k, int p, MPI_Comm comm);


double findr1(double* a, double* b, double* x,
  int n, int m, int k, int p, MPI_Comm comm);
double findr2(double* x, int n, int m, int k, int p, MPI_Comm comm);


#endif
