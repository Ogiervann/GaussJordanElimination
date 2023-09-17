#include "basic_funcs.h"

int l2g(int /*n*/, int m, int p, int k, int i_loc){
  int i_loc_m = i_loc/m;
  int i_glob_m = i_loc_m * p + k;
  return i_glob_m * m + i_loc%m;
}

int g2l(int /*n*/, int m, int p, int &k, int i_glob){
  int i_glob_m = i_glob/m;
  int i_loc_m = i_glob_m/p;
  k = i_glob_m%p;
  return i_loc_m * m + i_glob%m;
}

int get_rows(int n, int m, int p, int k){
  int b = (n+m-1)/m;
  return (b%p <= k ? b/p : b/p + 1);
}

int get_max_rows(int n, int m, int p){
  int b = (n+m-1)/m;
  return b/p+1;  
}



int get_k(int /*n*/, int m, int p, int i_glob){
  int i_glob_m = i_glob/m; return i_glob_m%p;
}



double get_full_time(){
  struct timeval buf;
  gettimeofday(&buf, NULL);
  return buf.tv_sec + buf.tv_usec / 1.e6;
}

double get_CPU_time(){
  struct rusage buf;
  getrusage(RUSAGE_THREAD, &buf);
  return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec / 1.e6;
}
