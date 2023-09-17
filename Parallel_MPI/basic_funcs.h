#ifndef BASIC_FUNCS_H
#define BASIC_FUNCS_H
#include "mpi.h"
#include <cstdio>
#include <cstdlib>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <sys/resource.h>



int l2g(int n, int m, int p, int k, int i_loc);
int g2l(int n, int m, int p, int &k, int i_glob);

int get_k(int /*n*/, int m, int p, int i_glob);
int get_rows(int n, int m, int p, int k);
int get_max_rows(int n, int m, int p);

double get_full_time();
double get_CPU_time();

#endif
