#include "basic_funcs.h"
void reduce_sum(int p, double* a, int n){
  static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
  static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
  static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
  static int t_in = 0;
  static int t_out = 0;
  static double*r = nullptr;
  int i;
  if(p <= 1)return;
  pthread_mutex_lock(&m);
  if(r == nullptr){
    r = a;
  }
  else{
    for(i = 0; i < n; i++) r[i]+=a[i];
  }
  t_in++;
  if(t_in >= p){
    t_out = 0;
    pthread_cond_broadcast(&c_in);
  }else{
    while(t_in<p){
      pthread_cond_wait(&c_in, &m);
    }
  }


  if(r != a){
    for(i = 0; i < n; i++) a[i] = r[i];
  }
  t_out++;
  if(t_out >=p){
    t_in = 0;
    r = nullptr;
    pthread_cond_broadcast(&c_out);
  }else{
    while(t_out < p){
      pthread_cond_wait(&c_out, &m);
    }
  }
  pthread_mutex_unlock(&m);
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
