#ifndef MATRINIT
#define MATRINIT
#include "mpi.h"
#include <cstring>
#include <cstdio>

void setId(double* a, int n);
void setPos(int* pos, int n);
void setPos_mpi(int* pos, int n, int m, int p, int k);

int readArray(FILE *fp, double* a, int n);
int printArray(double* a, int cols, int rows, int printed_rows, int max_print);



double initElement(int n,int s, int i, int j);

int readMatrix(double* a, int n, int m, int p, int k, const char* name,
  double* buff/* размера блочной строки */, MPI_Comm comm);

void initMatrix(double* a, int n, int m, int p, int k, int s);

void initVector(double* a, double* b, int n, int m, int t, int p);

void printMatrix(double *a, int n, int m, int p, int k,
  double* buf/*n*m buf*/, int max_print, MPI_Comm comm);


void printVector(double *b, int n, int m, int p, int k,
  double* buff/* m buff */, int max_print, MPI_Comm comm);

#endif
