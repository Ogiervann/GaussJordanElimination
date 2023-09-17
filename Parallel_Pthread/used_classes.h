#ifndef USED_CLASSES_H
#define USED_CLASSES_H
#include <pthread.h>

class Results{
public:
  int status = -100;
  double alg_CPU_time = 0;
  double alg_FULL_time = 0;
  double dis_CPU_time = 0;
  double dis_FULL_time = 0;
  double r1 = -1;
  double r2 = -1;
};


class maxBlockResults{
public:
  int pos = -1;
  double norm = -1;
};

class Args{
public:
  Results *res = nullptr;
  maxBlockResults *mbres = nullptr;
  int t = 0; int p = 0;
  int n = 0; int m = 0;
  int s = 0; int r = 0;
  char *filename = nullptr;
  double *a = nullptr;
  double *b = nullptr;
  double *x = nullptr;
  int *pos = nullptr;
  int *mpos = nullptr;
  double *block1 = nullptr;
  double *block2 = nullptr;
  double *block3 = nullptr;
  pthread_t tid = -1;
};

#endif
