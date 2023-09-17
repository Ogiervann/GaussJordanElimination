#ifndef BASIC_FUNCS_H
#define BASIC_FUNCS_H
#include <pthread.h>
#include <cstdio>
#include <cstdlib>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "used_classes.h"


void reduce_sum(int p, double* a = nullptr, int n = 0);

double get_full_time();
double get_CPU_time();

#endif
