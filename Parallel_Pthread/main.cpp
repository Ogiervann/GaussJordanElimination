#include <pthread.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include "my_thread.h"




int main(int argc, char* argv[]){
  int n, m, p, r, s, task = 17;
  char* filename = 0;
  //unsigned int start, end;

  double *a;
  double *b;
  double *x;
  int *pos;
  Args *args; Results *res;
  maxBlockResults *mbres;


  if(!((argc == 6 || argc == 7) && sscanf(argv[1], "%d", &n) == 1 && sscanf(argv[2], "%d", &m) == 1 && sscanf(argv[3], "%d", &p) == 1 && sscanf(argv[4], "%d", &r) == 1 && sscanf(argv[5], "%d", &s) == 1) || (s == 0 && argc == 6)){
    printf("Usage %s n m p r s (file)\n", argv[0]);
    return 1;
  }

  if(s == 0){
    filename = argv[6];
  }


  a = new double[n*n];
  b = new double[n];
  x = new double[n];
  pos = new int[n/m];
  args = new Args[p]; res = new Results[p];
  mbres = new maxBlockResults[p];

  if(a == NULL || b == NULL || x == NULL || pos == NULL || args == NULL || res == NULL || mbres == NULL){
    printf("Not enough memory\n");
    if(a != NULL) delete[] a;
    if(b != NULL) delete[] b;
    if(x != NULL) delete[] x;
    if(pos != NULL) delete[] pos;
    if(args != NULL) delete[] args;
    if(res != NULL) delete[] res;
    if(mbres != NULL) delete[] mbres;
    return 1;
  }

  setZero(a, n*n);
  setZero(b, n);
  setZero(x, n);

  setPos(pos, n/m);

  for(int k = 0; k < p; k++){
    args[k].res = res;
    args[k].mbres = mbres;
    args[k].filename = filename;
    args[k].a = a;
    args[k].b = b;
    args[k].x = x;

    args[k].r = r;
    args[k].s = s;
    args[k].n = n;
    args[k].m = m;
    args[k].t = k;
    args[k].p = p;
    args[k].pos = pos;
    args[k].mpos = new int[m];
    args[k].block1 = new double[m*m];
    args[k].block2 = new double[m*m];
    args[k].block3 = new double[m*m];
    if(args[k].block1 == NULL || args[k].block2 == NULL || args[k].block3 == NULL || args[k].mpos == NULL){
      printf("Not enough memory\n");
      for(int i = 0; i < k; i++){
        delete[] args[i].block1;
        delete[] args[i].block2;
        delete[] args[i].block3;
        delete[] args[i].mpos;
      }
      if(args[k].block1 != NULL) delete[] args[k].block1;
      if(args[k].block2 != NULL) delete[] args[k].block2;
      if(args[k].block3 != NULL) delete[] args[k].block3;
      if(args[k].mpos != NULL) delete[] args[k].mpos;
      delete[] a;
      delete[] b;
      delete[] x;
      delete[] pos;
      delete[] res;
      delete[] mbres;
      delete[] args;
      return 1;
    }
  }


  //start = get_full_time();
  for(int k = 1; k < p; k++){
    if(pthread_create(&args[k].tid, 0, thread_func, args+k)){
      printf("Cannot create thread %d\n", k);
      std::abort();
    }//if
  }//for
  thread_func(args+0);
  for (int k = 1; k < p; k++) {
    pthread_join(args[k].tid, 0);
  }
  //end = get_full_time();
  //double global_time = end-start;

  if(res[0].status == -1){
    printf("Error with file\n");
    for(int i = 0; i < p; i++){
      delete[] args[i].block1;
      delete[] args[i].block2;
      delete[] args[i].block3;
      delete[] args[i].mpos;
    }
    delete[] a;
    delete[] b;
    delete[] x;
    delete[] pos;
    delete[] res;
    delete[] mbres;
    delete[] args;
    return -1;
  }


  if(res[0].status == 0){
    printf("RES:\n");
    printMatrix(x, 1, n, r);
  }

  printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n",
    argv[0], task, res[0].r1, res[0].r2, res[0].alg_FULL_time, res[0].dis_FULL_time, s, n, m, p);

  for(int i = 0; i < p; i++){
    delete[] args[i].block1;
    delete[] args[i].block2;
    delete[] args[i].block3;
    delete[] args[i].mpos;
  }
  delete[] a;
  delete[] b;
  delete[] x;
  delete[] pos;
  delete[] res;
  delete[] mbres;
  delete[] args;
  return 0;
}
