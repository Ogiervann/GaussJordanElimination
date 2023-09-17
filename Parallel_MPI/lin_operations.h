#ifndef MATOPERATIONS
#define MATOPERATIONS
#include "mpi.h"


void get_block(int i, int j, int n, int m, int p, int k, double *a, double *block);
void put_block(int i, int j, int n, int m, int p, int k, double *a, double *block);

void get_vect_block(int i, int  n, int m, int p, int k, double *b, double *block);
void put_vect_block(int i, int  n, int m, int p, int k, double *b, double *block);


void get_buff_block(int j, int rows, int n, int m,
  double *buff, double *block);
void put_buff_block(int j, int rows, int n, int m,
  double *buff, double *block);

void switchBlockLines(int i, int j, double *a, double* b, double* block1,
  double* block2, int n, int m, int k, int p, MPI_Comm comm);

void switchBlockColumns(int i, int j, double *a, int *pos, double* block1,
  double* block2, int n, int m, int k, int p, MPI_Comm comm);

#endif
