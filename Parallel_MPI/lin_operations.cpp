#include "lin_operations.h"
#include "basic_funcs.h"
#include <cmath>


void get_block(int i, int j, int n, int m, int p, int k, double *a, double *block){
  int s, t;
  int i_glob = l2g(n, m, p, k, i*m)/m;
  int h1 = (i_glob*m+m<n?m:n-i_glob*m);
  int h2 = (j*m+m<n?m:n-j*m);

  for(s = 0; s < h1; s++){
    for(t = 0; t < h2; t++){
      block[s*h2 + t] = a[i*n*m + s*n + j*m+t];
    }
  }
}

void put_block(int i, int j, int n, int m, int p, int k, double *a, double *block){
  int s, t;
  int i_glob = l2g(n, m, p, k, i*m)/m;
  int h1 = (i_glob*m+m<n?m:n-i_glob*m);
  int h2 = (j*m+m<n?m:n-j*m);

  for(s = 0; s < h1; s++){
    for(t = 0; t < h2; t++){
      a[i*n*m + s*n + j*m+t] = block[s*h2 + t];
    }
  }
}

void get_buff_block(int j, int rows, int n, int m,
  double *buff, double *block){
  int s, t, h = (j*m+m < n? m:n-j*m);
  for(s = 0; s < rows; s++){
    for(t = 0; t < h; t++){
      block[s*h + t] = buff[s*n + j*m+t];
    }
  }
}

void put_buff_block(int j, int rows, int n, int m,
  double *buff, double *block){
  int s, t, h = (j*m+m < n? m:n-j*m);
  for(s = 0; s < rows; s++){
    for(t = 0; t < h; t++){
      buff[s*n + j*m+t] = block[s*h + t];
    }
  }
}



void get_vect_block(int i, int  n, int m, int p, int k, double *b, double *block){
  int s;
  int i_glob = l2g(n, m, p, k, i*m)/m;
  int h = (i_glob*m+m<n?m:n-m*i_glob);
  for(s = 0; s < h; s++){
    block[s] = b[i*m + s];
  }
}

void put_vect_block(int i, int  n, int m, int p, int k, double *b, double *block){
  int s;
  int i_glob = l2g(n, m, p, k, i*m)/m;
  int h = (i_glob*m+m<n?m:n-m*i_glob);
  for(s = 0; s < h; s++){
    b[i*m + s] = block[s];
  }
}


void switchBlockLines(int i, int j, double *a, double* b, double* block1,
  double* block2, int n, int m, int k, int p, MPI_Comm comm){
  int s;
  if(i == j){
    return;
  }

  int max_b = (n + m - 1)/m;

  int k_i, k_j;

  int i_loc = g2l(n, m, p, k_i, i*m)/m;
  int j_loc = g2l(n, m, p, k_j, j*m)/m;

  if(k_i == k_j && k_i == k){
    for(s = 0; s < max_b; s++){
      get_block(i_loc, s, n, m, p, k, a, block1);
      get_block(j_loc, s, n, m, p, k, a, block2);
      put_block(i_loc, s, n, m, p, k, a, block2);
      put_block(j_loc, s, n, m, p, k, a, block1);
    }
    get_vect_block(i_loc, n, m, p, k, b, block1);
    get_vect_block(j_loc, n, m, p, k, b, block2);
    put_vect_block(i_loc, n, m, p, k, b, block2);
    put_vect_block(j_loc, n, m, p, k, b, block1);

  }
  else if(k == k_i){
    MPI_Status st;
    MPI_Sendrecv_replace(a+i_loc*n*m, n*m, MPI_DOUBLE, k_j, 0/*tag*/,
      k_j, 0/*tag*/, comm, &st);
    MPI_Sendrecv_replace(b+i_loc*m, m, MPI_DOUBLE, k_j, 0/*tag*/,
      k_j, 0/*tag*/, comm, &st);
  }
  else if(k == k_j){
    MPI_Status st;
    MPI_Sendrecv_replace(a+j_loc*n*m, n*m, MPI_DOUBLE, k_i, 0/*tag*/,
      k_i, 0/*tag*/, comm, &st);
    MPI_Sendrecv_replace(b+j_loc*m, m, MPI_DOUBLE, k_i, 0/*tag*/,
      k_i, 0/*tag*/, comm, &st);
  }
}

void switchBlockColumns(int i, int j, double *a, int *pos, double* block1,
  double* block2, int n, int m, int k, int p, MPI_Comm comm){
  int s;
  if(i == j){
    return;
  }

  int rows = get_rows(n, m, p, k);

  int k_i, k_j;

  int i_loc = g2l(n, m, p, k_i, i*m)/m;
  int j_loc = g2l(n, m, p, k_j, j*m)/m;
  //printf("i: %d j: %d i_loc: %d j_loc: %d k_i: %d k_j: %d| %d\n", i, j, i_loc, j_loc, k_i, k_j, k);
  for(s = 0; s < rows; s++){
    get_block(s, i, n, m, p, k, a, block1);
    get_block(s, j, n, m, p, k, a, block2);
    put_block(s, i, n, m, p, k, a, block2);
    put_block(s, j, n, m, p, k, a, block1);
  }

  if(k_i == k_j && k_i == k){
    s = pos[i_loc];
    pos[i_loc] = pos[j_loc];
    pos[j_loc] = s;
  }
  else if(k == k_i){
    MPI_Status st;
    //printf("pos[i_loc]_old: %d | %d\n", pos[i_loc], k);
    MPI_Sendrecv_replace(pos+i_loc, 1, MPI_INT, k_j, 0/*tag*/,
      k_j, 0/*tag*/, comm, &st);
    //printf("pos[i_loc]_new: %d | %d\n", pos[i_loc], k);
  }
  else if(k == k_j){
    MPI_Status st;
    //printf("pos[j_loc]_old: %d | %d\n", pos[j_loc], k);
    MPI_Sendrecv_replace(pos+j_loc, 1, MPI_INT, k_i, 0/*tag*/,
      k_i, 0/*tag*/, comm, &st);
    //printf("pos[j_loc]_new: %d | %d\n", pos[j_loc], k);
  }
}
