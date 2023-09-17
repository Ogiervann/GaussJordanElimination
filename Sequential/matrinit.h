#ifndef MATRINIT
#define MATRINIT
#include <cstdio>

void setZero(double *a, int n);
void setId(double* a, int n);
void setPos(int* pos, int n);


double initElement(int n,int s, int i, int j);

int readMatrix(double* a, int n, char* filename);

int initMatrix(double* a, int n, int s, char* filename);

void initVector(double* a, double* b, int n);

void printMatrix(double* a, int l, int n, int r);
#endif
