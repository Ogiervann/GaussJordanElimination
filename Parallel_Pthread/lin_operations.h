#ifndef MATOPERATIONS
#define MATOPERATIONS

void get_block(int i, int j, int n, int m, int k, int l, double *a, double *block);

void put_block(int i, int j, int n, int m, int k, int l, double *a, double *block);

void get_vect_block(int i, int m, int k, int l, double *b, double *block);
void put_vect_block(int i, int m, int k, int l, double *b, double *block);

double norml1(double* a, int n);

double normMatr(double* a, int n);


void matrixMult(int n, int s, int m, double *a, double* b, double* res);


void matrSubstr_eq(int n, int m, double* a, double* b);


void vectSubstr_eq(int n, double* a, double* b);


void switchLines(int i, int j, int n, double* a, double* b);

void switchColumns(int i, int j, int n, double* a, int* mpos);


void switchBlockLines(int i, int j, int n, int m, int k, int l, double *a,
  double* b, double* block1, double* block2, int t, int p);

void switchBlockColumns(int i, int j, int n, int m, int k, int l, double *a,
  int *pos, double* block1, double* block2, int t, int p);
#endif
