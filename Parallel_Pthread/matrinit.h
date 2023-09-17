#ifndef MATRINIT
#define MATRINIT
#include <cstdio>

void setZero(double *a, int n);
void setId(double* a, int n);
void setPos(int* pos, int n);


double initElement(int n,int s, int i, int j);

int readMatrix(double* a, int n, char* filename);

int initMatrix(double* a, int n, int m, int s, char* filename, int t, int p);

void initVector(double* a, double* b, int n, int m, int t, int p);

void printMatrix(double* a, int l, int n, int r);
#endif
