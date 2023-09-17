#include <iostream>
#include <cstdio>
#include <sys/time.h>
#include <cstring>
#include "linearsolve.h"

#define EPS 1e-16

double get_full_time(){
  struct timeval buf;
  gettimeofday(&buf, NULL);
  return buf.tv_sec + buf.tv_usec / 1.e6;
}


int main(int argc, char* argv[]){
  int n, m, r, s, res, task = 17;
  char* filename = 0;
  double r1 = -1, r2 = -1, t1 = 0, t2 = 0, normA;
  double start, end;

  double *a;
  double *b;
  double *x;
  int *pos;
  int *mpos;
  double *block1;
  double *block2;
  double *block3;


  if(!((argc == 5 || argc == 6) && sscanf(argv[1], "%d", &n) == 1 && sscanf(argv[2], "%d", &m) == 1 && sscanf(argv[3], "%d", &r) == 1 && sscanf(argv[4], "%d", &s) == 1) || (s == 0 && argc == 5)){
    printf("Usage %s n m r s (file)\n", argv[0]);
    return 1;
  }

  if(s == 0){
    filename = argv[5];
  }

  start = get_full_time();

  a = new double[n*n];
  if(a == 0){
    printf("Not enough memory\n");
    return 1;
  }
  memset(a, 0, n*n*sizeof(double));
  res = initMatrix(a, n, s, filename);

  if(res == 1){
    delete[] a;
    printf("File not found or data in the file is incorrect.\n");
    return 1;
  }

  printMatrix(a, n, n, r);
  printf("\n");


  b = new double[n];
  x = new double[n];
  pos = new int[n/m];
  mpos = new int[m];
  block1 = new double[m*m];
  block2 = new double[m*m];
  block3 = new double[m*m];

  if(b == 0 || x == 0 || pos == 0 || block1 == 0 || block2 == 0 || block3 == 0 || mpos == 0){
    printf("Not enough memory\n");
    if(a != 0)delete[] a;
    if(b != 0)delete[] b;
    if(x != 0)delete[] x;
    if(block1 != 0)delete[] block1;
    if(block2 != 0)delete[] block2;
    if(block3 != 0)delete[] block3;
    if(pos != 0)delete[] pos;
    if(mpos != 0)delete[] mpos;
    return 1;
  }
  memset(b, 0, n*sizeof(double));
  memset(x, 0, n*sizeof(double));
  setPos(pos, n/m);
  memset(block1, 0, m*m*sizeof(double));
  memset(block2, 0, m*m*sizeof(double));
  memset(block3, 0, m*m*sizeof(double));
  initVector(a, b, n);

  /*
  printMatrix(a, n, n, r);
  printf("\n");
  printMatrix(b, n, 1, r);
  printf("\n");

  printf("\n");
  matrixInverse(n, a, block1, block2);
  printMatrix(block1, n, n, r);
  printf("\n");

  initMatrix(a, n, s, filename);
*/
  normA = normMatr(a, n);
  res = linearSolve(n, m, a, b, x, pos, block1, block2, block3, mpos, normA * EPS);
  end = get_full_time();

  t1 = end - start;
  if(res == 0){
    printMatrix(x, 1, n, r);

    setZero(b, n);

    initMatrix(a, n, s, filename);
    initVector(a, b, n);
/*
    printMatrix(a, n, n, r);
    printf("\n");
    printMatrix(b, n, 1, r);
    printf("\n");
*/

    start = get_full_time();

    r1 = findr1(n, m, a, b, x, block1, block2, block3);
    r2 = findr2(n, x);
    end = get_full_time();

    //std::cout << start << " " << end << std::endl;

    t2 = end - start;

  }

  printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n",
argv[0], task, r1, r2, t1, t2, s, n, m);

  delete[] a;
  delete[] b;
  delete[] x;
  delete[] block1;
  delete[] block3;
  delete[] block2;
  delete[] pos;
  delete[] mpos;
  return 0;
}
