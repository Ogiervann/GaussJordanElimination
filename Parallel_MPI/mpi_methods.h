#ifndef MPI_METHODS_H
#define MPI_METHODS_H
#include "mpi.h"
#include "basic_funcs.h"
#include <cstring>
#include <cmath>

void matrix_mult_vector(double* a, double* b, double*c, int n, int m, int k,
  int p, MPI_Comm comm);

double mpi_normMatr(double* a, int n, int m, int k, int p, MPI_Comm comm);
double mpi_norml1(double* b, int n, int m, int k, int p, MPI_Comm comm);


#endif
