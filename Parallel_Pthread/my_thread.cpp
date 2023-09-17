#include "my_thread.h"
#include <cmath>

void* thread_func(void* ptr){
  double start_CPU, start_FULL, end_CPU, end_FULL;
  Args* args = (Args*) ptr;
  double* a = args->a;
  double *b = args->b;
  double *x = args->x;
  double *block1 = args->block1;
  double *block2 = args->block2;
  double *block3 = args->block3;
  int *pos = args->pos;
  int *mpos = args->mpos;
  int n = args->n;
  int m = args->m;
  int p = args->p;
  int t = args->t;
  int s = args->s;
  int r = args->r;
  char *filename = args->filename;
  Results *res = args->res;
  maxBlockResults *mbres = args->mbres;
  int status = 0;
  double r1 = -1, r2 = -1;

  cpu_set_t cpu;
  CPU_ZERO(&cpu);
  int n_cpus = get_nprocs();
  int cpu_id = n_cpus - 1 - (t % n_cpus);
  CPU_SET(cpu_id, &cpu);
  pthread_t tid = pthread_self();
  pthread_setaffinity_np(tid, sizeof(cpu), &cpu);
  //printf("here%d\n", t);
  for(int i = t * m; i < n; i += p * m) {
      int h = (i + m <= n? m : n - i);
      memset(a + i * n, 0, h * n * sizeof(double));
      memset(x + i , 0, h * sizeof(double));
      memset(b + i , 0, h * sizeof(double));
  }
  reduce_sum(p);

  start_CPU = get_CPU_time();
  start_FULL = get_full_time();
  status = initMatrix(a, n, m, s, filename, t, p);

  if(status == 1){
    end_CPU = get_CPU_time();
    end_FULL = get_full_time();
    res[t].alg_CPU_time = end_CPU-start_CPU;
    res[t].alg_FULL_time = end_FULL-start_FULL;
    res[t].status = -1;
    return 0;
  }
  if(t == 0){
    printf("A:\n");
    printMatrix(a, n, n, r);
    printf("\n");
  }

  initVector(a, b, n, m, t, p);
  if(t == 0){
    printf("b:\n");

    printMatrix(b, 1, n, r);
    printf("\n");
  }
  double normA = normMatr(a, n);
  status = linearSolve(n, m, a, b, x, pos, block1, block2, block3, mpos, mbres, normA * EPS, t, p);
  end_CPU = get_CPU_time();
  end_FULL = get_full_time();
  res[t].alg_CPU_time = end_CPU-start_CPU;
  res[t].alg_FULL_time = end_FULL-start_FULL;
  res[t].status = status;

  if(status == 0){
    initMatrix(a, n, m, s, filename, t, p);
    initVector(a, b, n, m, t, p);
    //if(t==0)printMatrix(a, n, n, n);
    //if(t == 0)printMatrix(b, 1, n, n);
    start_CPU = get_CPU_time();
    start_FULL = get_full_time();
    r1 = findr1(n, m, a, b, x, block1, block2, block3, t, p);
    r2 = findr2(n, x);
    end_CPU = get_CPU_time();
    end_FULL = get_full_time();
    res[t].dis_CPU_time = end_CPU-start_CPU;
    res[t].dis_FULL_time = end_FULL-start_FULL;
    res[t].r1 = r1;
    res[t].r2 = r2;
  }

  return 0;
}
