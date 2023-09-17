#include "mpi_methods.h"

void matrix_mult_vector(double* a, double* b, double*c, int n, int m, int k,
  int p, MPI_Comm comm){

  int rows = get_rows(n, m, p, k);
  int max_rows = get_max_rows(n, m, p);
  int src = (k+1+p)%p;
  int dst = (k-1+p)%p;
  memset(c, 0, rows*m*sizeof(double));
  for(int l = 0; l < p; l++){
    int lk = (k+l)%p;
    int lk_rows = get_rows(n, m, p, lk);
    for(int lk_i = 0; lk_i < lk_rows; lk_i++){
      int lk_ig = l2g(n,m,p,lk,lk_i);
      for(int i = 0; i < rows; i++){
        //c_i += A_i,lk * b_lk
        int lk_m = (lk_i*m+m < n?m:n-lk_i*m);
        int i_m = (i*m+m < n?m:n-i*m);
        for(int ii = 0; ii<i_m; ii++){
          double s = 0;
          for(int jj = 0; jj < lk_m; jj++){
            s += a[i*m*n+ii*n+lk_ig*m+jj]*b[lk_i*m+jj];
          }
          c[i*m+ii] += s;
        }
      }
    }
    MPI_Status st;
    MPI_Sendrecv_replace(b, max_rows*m, MPI_DOUBLE, dst, 0/*tag*/,
      src, 0/*tag*/, comm, &st);
  }

}

double mpi_normMatr(double* a, int n, int m, int k, int p, MPI_Comm comm){
  int rows = get_rows(n, m, p, k);
  double norm_k = -1, norm = -1, t;
  for(int i = 0; i < rows; i++){
    t = 0;
    int h = (i*m+m<n?m:n-i*m);
    for(int i1 = 0; i1 < h; i1++){
      for(int j = 0; j < n; j++){
        t += fabs(a[i*m*n+i1*n+j]);
      }
      if(t > norm_k) {
        norm_k = t;
      }
    }
  }

  MPI_Allreduce(&norm_k, &norm, 1, MPI_DOUBLE, MPI_MAX, comm);
  return norm;
}


double mpi_norml1(double* b, int n, int m, int k, int p, MPI_Comm comm){
  int rows = get_rows(n, m, p, k);
  double norm_k = 0, norm = -1;
  for(int i = 0; i < rows; i++){
      int h = (i*m+m<n?m:n-i*m);
      for(int i1 = 0; i1 < h; i1++){
        norm_k += fabs(b[i*m+i1]);
      }
  }

  MPI_Allreduce(&norm_k, &norm, 1, MPI_DOUBLE, MPI_SUM, comm);
  return norm;
}
