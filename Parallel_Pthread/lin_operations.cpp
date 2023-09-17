#include "lin_operations.h"
#include "basic_funcs.h"
#include <cmath>


void get_block(int i, int j, int n, int m, int k, int l, double *a, double *block){
  int s, t, n1 = m, n2 = m;
  if(i == k){
    n1 = l;
  }
  if( j == k){
    n2 = l;
  }
  for(s = 0; s < n1; s++){
    for(t = 0; t < n2; t++){
      block[s*n2 + t] = a[i*n*m + s*n + j*m+t];
    }
  }
}

void put_block(int i, int j, int n, int m, int k, int l, double *a, double *block){
  int s, t, n1 = m, n2 = m;
  if(i == k){
    n1 = l;
  }
  if( j == k){
    n2 = l;
  }
  for(s = 0; s < n1; s++){
    for(t = 0; t < n2; t++){
      a[i*n*m + s*n + j*m+t] = block[s*n2 + t];
    }
  }
}

void get_vect_block(int i, int m, int k, int l, double *b, double *block){
  int s, n1 = m;
  if(i == k){
    n1 = l;
  }
  for(s = 0; s < n1; s++){
    block[s] = b[i*m + s];
  }
}

void put_vect_block(int i, int m, int k, int l, double *b, double *block){
  int s, n1 = m;
  if(i == k){
    n1 = l;
  }
  for(s = 0; s < n1; s++){
    b[i*m + s] = block[s];
  }
}


double norml1(double* a, int n){
  int i, k, l;
  double res = 0;

  k = n / 3;
  l = n - k*3;

  for(i = 0; i < k; i++){
    res += fabs(a[i*3]) + fabs(a[i*3+1]) + fabs(a[i*3+2]);
  }
  for(i = 0; i < l; i++){
    res += fabs(a[k*3+i]);
  }
  return res;

}



double normMatr(double* a, int n){
  double norm = 0, tmp0, tmp1, tmp2, t00, t01, t02, t10, t11, t12, t20, t21, t22;
  int i, j, k, l;
  k = n/3;
  l = n - k*3;
  for(i = 0; i < k; i++){
    tmp0 = 0;
    tmp1 = 0;
    tmp2 = 0;
    t00 = 0;
    t01 = 0;
    t02 = 0;
    t10 = 0;
    t11 = 0;
    t12 = 0;
    t20 = 0;
    t21 = 0;
    t22 = 0;
    for(j = 0; j < k; j++){
      t00 += fabs(a[3*i*n+3*j]);
      t01 += fabs(a[3*i*n+3*j+1]);
      t02 += fabs(a[3*i*n+3*j+2]);
      t10 += fabs(a[(3*i+1)*n+3*j]);
      t11 += fabs(a[(3*i+1)*n+3*j+1]);
      t12 += fabs(a[(3*i+1)*n+3*j+2]);
      t20 += fabs(a[(3*i+2)*n+3*j]);
      t21 += fabs(a[(3*i+2)*n+3*j+1]);
      t22 += fabs(a[(3*i+2)*n+3*j+2]);
    }
    for(j = 0; j < l; j++){
      tmp0 += fabs(a[3*i*n+3*k+j]);
      tmp1 += fabs(a[(3*i+1)*n+3*k+j]);
      tmp2 += fabs(a[(3*i+2)*n+3*k+j]);
    }
    tmp0 = tmp0 + t00 + t01 + t02;
    tmp1 = tmp1 + t10 + t11 + t12;
    tmp2 = tmp2 + t20 + t21 + t22;
    tmp0 = (norm>tmp0?norm:tmp0);
    tmp1 = (tmp1>tmp2?tmp1:tmp2);
    norm = (tmp0 > tmp1?tmp0:tmp1);
  }
  tmp0 = 0;
  tmp1 = 0;
  tmp2 = 0;
  t00 = 0;
  t01 = 0;
  t02 = 0;
  for(i = 0; i < l; i++){
    for(j = 0; j < k; j++){
      t00 += fabs(a[(i+3*k)*n+3*j]);
      t01 += fabs(a[(i+3*k)*n+3*j+1]);
      t02 += fabs(a[(i+3*k)*n+3*j+2]);
    }
    for(j = 0; j < l; j++){
      //printf("Here10 %lf\n", a[0]);
      tmp0 += fabs(a[(i+3*k)*n+3*k+j]);
    }
    tmp0 = tmp0 + t00 + t01 + t02;
    norm = (norm>tmp0?norm:tmp0);
    //printf("%lf %lf", tmp0, norm);

  }
  return norm;
}


void matrixMult(int n, int s, int m, double *a, double* b, double* res){
  int k1, l1, k2, l2, k3, l3, i, j, t;
  double t00, t01, t02, t10, t11, t12, t20, t21, t22;
  k1 = n/3;
  k2 = s/3;
  k3 = m/3;
  l1 = n - k1*3;
  l2 = s - k2*3;
  l3 = m - k3*3;

  for(i = 0; i < k1; i++){
    for(j = 0; j < k3; j++){
      t00 = 0;
      t01 = 0;
      t02 = 0;
      t10 = 0;
      t11 = 0;
      t12 = 0;
      t20 = 0;
      t21 = 0;
      t22 = 0;
      for(t = 0; t < k2; t++){
        t00 += a[3*i*s+3*t]*b[3*t*m+3*j] + a[3*i*s+3*t+1]*b[(3*t+1)*m+3*j] + a[3*i*s+3*t+2]*b[(3*t+2)*m+3*j];
        t01 += a[3*i*s+3*t]*b[3*t*m+3*j+1] + a[3*i*s+3*t+1]*b[(3*t+1)*m+3*j+1] + a[3*i*s+3*t+2]*b[(3*t+2)*m+3*j+1];
        t02 += a[3*i*s+3*t]*b[3*t*m+3*j+2] + a[3*i*s+3*t+1]*b[(3*t+1)*m+3*j+2] + a[3*i*s+3*t+2]*b[(3*t+2)*m+3*j+2];
        t10 += a[(3*i+1)*s+3*t]*b[3*t*m+3*j] + a[(3*i+1)*s+3*t+1]*b[(3*t+1)*m+3*j] + a[(3*i+1)*s+3*t+2]*b[(3*t+2)*m+3*j];
        t11 += a[(3*i+1)*s+3*t]*b[3*t*m+3*j+1] + a[(3*i+1)*s+3*t+1]*b[(3*t+1)*m+3*j+1] + a[(3*i+1)*s+3*t+2]*b[(3*t+2)*m+3*j+1];
        t12 += a[(3*i+1)*s+3*t]*b[3*t*m+3*j+2] + a[(3*i+1)*s+3*t+1]*b[(3*t+1)*m+3*j+2] + a[(3*i+1)*s+3*t+2]*b[(3*t+2)*m+3*j+2];
        t20 += a[(3*i+2)*s+3*t]*b[3*t*m+3*j] + a[(3*i+2)*s+3*t+1]*b[(3*t+1)*m+3*j] + a[(3*i+2)*s+3*t+2]*b[(3*t+2)*m+3*j];
        t21 += a[(3*i+2)*s+3*t]*b[3*t*m+3*j+1] + a[(3*i+2)*s+3*t+1]*b[(3*t+1)*m+3*j+1] + a[(3*i+2)*s+3*t+2]*b[(3*t+2)*m+3*j+1];
        t22 += a[(3*i+2)*s+3*t]*b[3*t*m+3*j+2] + a[(3*i+2)*s+3*t+1]*b[(3*t+1)*m+3*j+2] + a[(3*i+2)*s+3*t+2]*b[(3*t+2)*m+3*j+2];
      }
      for(t = 0; t < l2; t++){
        t00 += a[3*i*s+t+3*k2]*b[(t+3*k2)*m+3*j];
        t01 += a[3*i*s+(t+3*k2)]*b[(t+3*k2)*m+3*j+1];
        t02 += a[3*i*s+(t+3*k2)]*b[(t+3*k2)*m+3*j+2];
        t10 += a[(3*i+1)*s+(t+3*k2)]*b[(t+3*k2)*m+3*j];
        t11 += a[(3*i+1)*s+(t+3*k2)]*b[(t+3*k2)*m+3*j+1];
        t12 += a[(3*i+1)*s+(t+3*k2)]*b[(t+3*k2)*m+3*j+2];
        t20 += a[(3*i+2)*s+(t+3*k2)]*b[(t+3*k2)*m+3*j];
        t21 += a[(3*i+2)*s+(t+3*k2)]*b[(t+3*k2)*m+3*j+1];
        t22 += a[(3*i+2)*s+(t+3*k2)]*b[(t+3*k2)*m+3*j+2];
      }
      res[3*i*m+3*j] = t00;
      res[3*i*m+3*j+1] = t01;
      res[3*i*m+3*j+2] = t02;
      res[(3*i+1)*m+3*j] = t10;
      res[(3*i+1)*m+3*j+1] = t11;
      res[(3*i+1)*m+3*j+2] = t12;
      res[(3*i+2)*m+3*j] = t20;
      res[(3*i+2)*m+3*j+1] = t21;
      res[(3*i+2)*m+3*j+2] = t22;
    }

    for(j = 0; j < l3; j++){
      t00 = 0;
      t10 = 0;
      t20 = 0;
      for(t = 0; t < k2; t++){

        t00 += a[3*i*s+3*t]*b[3*t*m+j+3*k3] + a[3*i*s+3*t+1]*b[(3*t+1)*m+j+3*k3] + a[3*i*s+3*t+2]*b[(3*t+2)*m+j+3*k3];
        t10 += a[(3*i+1)*s+3*t]*b[3*t*m+j+3*k3] + a[(3*i+1)*s+3*t+1]*b[(3*t+1)*m+j+3*k3] + a[(3*i+1)*s+3*t+2]*b[(3*t+2)*m+j+3*k3];
        t20 += a[(3*i+2)*s+3*t]*b[3*t*m+j+3*k3] + a[(3*i+2)*s+3*t+1]*b[(3*t+1)*m+j+3*k3] + a[(3*i+2)*s+3*t+2]*b[(3*t+2)*m+j+3*k3];

      }
      for(t = 0; t < l2; t++){
        t00 += a[3*i*s+t+3*k2]*b[(t+3*k2)*m+j+3*k3];
        t10 += a[(3*i+1)*s+(t+3*k2)]*b[(t+3*k2)*m+j+3*k3];
        t20 += a[(3*i+2)*s+(t+3*k2)]*b[(t+3*k2)*m+j+3*k3];
      }
      res[3*i*m+j+3*k3] = t00;
      res[(3*i+1)*m+j+3*k3] = t10;
      res[(3*i+2)*m+j+3*k3] = t20;
    }
  }

  for(i = 0; i < l1; i++){
    for(j = 0; j < k3; j++){
      t00 = 0;
      t01 = 0;
      t02 = 0;
      for(t = 0; t < k2; t++){
        t00 += a[(i+3*k1)*s+3*t]*b[3*t*m+3*j] + a[(i+3*k1)*s+3*t+1]*b[(3*t+1)*m+3*j] + a[(i+3*k1)*s+3*t+2]*b[(3*t+2)*m+3*j];
        t01 += a[(i+3*k1)*s+3*t]*b[3*t*m+3*j+1] + a[(i+3*k1)*s+3*t+1]*b[(3*t+1)*m+3*j+1] + a[(i+3*k1)*s+3*t+2]*b[(3*t+2)*m+3*j+1];
        t02 += a[(i+3*k1)*s+3*t]*b[3*t*m+3*j+2] + a[(i+3*k1)*s+3*t+1]*b[(3*t+1)*m+3*j+2] + a[(i+3*k1)*s+3*t+2]*b[(3*t+2)*m+3*j+2];

      }
      for(t = 0; t < l2; t++){
        t00 += a[(i+3*k1)*s+t+3*k2]*b[(t+3*k2)*m+3*j];
        t01 += a[(i+3*k1)*s+(t+3*k2)]*b[(t+3*k2)*m+3*j+1];
        t02 += a[(i+3*k1)*s+(t+3*k2)]*b[(t+3*k2)*m+3*j+2];
      }
      res[(i+3*k1)*m+3*j] = t00;
      res[(i+3*k1)*m+3*j+1] = t01;
      res[(i+3*k1)*m+3*j+2] = t02;
    }
    for(j = 0; j < l3; j++){
      t00 = 0;
      for(t = 0; t < k2; t++){
        t00 += a[(i+3*k1)*s+3*t]*b[3*t*m+j+3*k3] + a[(i+3*k1)*s+3*t+1]*b[(3*t+1)*m+j+3*k3] + a[(i+3*k1)*s+3*t+2]*b[(3*t+2)*m+j+3*k3];
      }
      for(t = 0; t < l2; t++){
        t00 += a[(i+3*k1)*s+t+3*k2]*b[(t+3*k2)*m+j+3*k3];
      }
      res[(i+3*k1)*m+j+3*k3] = t00;
    }
  }
}


void matrSubstr_eq(int n, int m, double* a, double* b){
  int i, j, k1, l1, k2, l2;
  k1 = n/3;
  l1 = n - k1*3;
  k2 = m/3;
  l2 = m - k2*3;
  for(i = 0; i < k1; i++){
    for(j = 0; j < k2; j++){
      a[3*i*m+3*j] -= b[3*i*m+3*j];
      a[3*i*m+3*j+1] -= b[3*i*m+3*j+1];
      a[3*i*m+3*j+2] -= b[3*i*m+3*j+2];
      a[(3*i+1)*m+3*j] -= b[(3*i+1)*m+3*j];
      a[(3*i+1)*m+3*j+1] -= b[(3*i+1)*m+3*j+1];
      a[(3*i+1)*m+3*j+2] -= b[(3*i+1)*m+3*j+2];
      a[(3*i+2)*m+3*j] -= b[(3*i+2)*m+3*j];
      a[(3*i+2)*m+3*j+1] -= b[(3*i+2)*m+3*j+1];
      a[(3*i+2)*m+3*j+2] -= b[(3*i+2)*m+3*j+2];
    }
    for(j = 0; j < l2; j++){
      a[3*i*m+3*k2+j] -= b[3*i*m+3*k2+j];
      a[(3*i+1)*m+3*k2+j] -= b[(3*i+1)*m+3*k2+j];
      a[(3*i+2)*m+3*k2+j] -= b[(3*i+2)*m+3*k2+j];
    }
  }
  for(i = 0; i < l1; i++){
    for(j = 0; j < k2; j++){
      a[(i+3*k1)*m+3*j] -= b[(i+3*k1)*m+3*j];
      a[(i+3*k1)*m+3*j+1] -= b[(i+3*k1)*m+3*j+1];
      a[(i+3*k1)*m+3*j+2] -= b[(i+3*k1)*m+3*j+2];
    }
    for(j = 0; j < l2; j++){
      a[(i+3*k1)*m+3*k2+j] -= b[(i+3*k1)*m+3*k2+j];
    }
  }
}

void vectSubstr_eq(int n, double* a, double* b){
  int i, k, l;
  k = n/3;
  l = n - k*3;

  for(i = 0; i < k; i++){
    a[3*i] -= b[3*i];
    a[3*i+1] -= b[3*i+1];
    a[3*i+2] -= b[3*i+2];
  }
  for(i = 0; i < l; i++){
    a[3*k+i] -= b[3*k+i];
  }
}


void switchLines(int i, int j, int n, double* a, double* b){
  double k;
  int t;
  if(i == j){
    return;
  }
  for(t = 0; t < n; t++){
    k = a[i*n+t];
    a[i*n+t] = a[j*n+t];
    a[j*n+t] = k;
    k = b[i*n+t];
    b[i*n+t] = b[j*n+t];
    b[j*n+t] = k;
  }
}

void switchColumns(int i, int j, int n, double* a, int* mpos){
  int k;
  double s;
  int t;
  if(i == j){
    return;
  }
  for(t = 0; t < n; t++){
    s = a[t*n+i];
    a[t*n+i] = a[t*n+j];
    a[t*n+j] = s;
  }
  k = mpos[i];
  mpos[i] = mpos[j];
  mpos[j] = k;
}


void switchBlockLines(int i, int j, int n, int m, int k, int l, double *a,
  double* b, double* block1, double* block2, int t, int p){
  int s;
  if(i == j){
    return;
  }
  for(s = t; s < k+1; s += p){
    get_block(i, s, n, m, k, l, a, block1);
    get_block(j, s, n, m, k, l, a, block2);
    put_block(i, s, n, m, k, l, a, block2);
    put_block(j, s, n, m, k, l, a, block1);
  }
  if(t == 0){
    get_vect_block(i, m, k, l, b, block1);
    get_vect_block(j, m, k, l, b, block2);
    put_vect_block(i, m, k, l, b, block2);
    put_vect_block(j, m, k, l, b, block1);
  }
  reduce_sum(p);
}

void switchBlockColumns(int i, int j, int n, int m, int k, int l, double *a,
  int *pos, double* block1, double* block2, int t, int p){
  int s;
  if(i == j){
    return;
  }
  for(s = t; s < k+1; s += p){
    get_block(s, i, n, m, k, l, a, block1);
    get_block(s, j, n, m, k, l, a, block2);
    put_block(s, i, n, m, k, l, a, block2);
    put_block(s, j, n, m, k, l, a, block1);
  }
  if(t == 0){
    s = pos[i];
    pos[i] = pos[j];
    pos[j] = s;
  }
  reduce_sum(p);
}
