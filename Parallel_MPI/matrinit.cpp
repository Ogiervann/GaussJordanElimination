#include "matrinit.h"
#include "basic_funcs.h"

void setId(double* a, int n){
  int i, j;
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      if(i != j){
        a[i*n+j] = 0;
      }
      else{
        a[i*n+j] = 1;
      }
    }
  }
}


void setPos(int* pos, int n){
  for(int i = 0; i < n; i++){
    pos[i] = i;
  }
}

void setPos_mpi(int* pos, int n, int m, int p, int k){
  int i;
  int rows = get_rows(n, m, p, k);
  for(i = 0; i < rows; i++){
    pos[i] = l2g(n, m, p, k, m*i)/m;
  }
}



double initElement(int n,int s, int i, int j){
  if(s == 1){
    return n - (i>j?i:j) + 1;
  }
  if(s == 2){
    return (i>j?i:j);
  }
  if(s == 3){
    return (i>j?i-j:j-i);
  }
  if(s==4){
    return 1/(i+j-1.0);
  }
  return 0;
}

void initMatrix(double* a, int n, int m, int p, int k, int s){
  int i_loc, i_loc_m , j_loc, i_glob, i_glob_0 , j_glob, rows;
  int h = m;
  rows = get_rows(n, m, p, k);
  //printf("%d\n", rows);
  for(i_loc_m = 0; i_loc_m < rows; i_loc_m++){
    i_glob_0 = l2g(n, m, p, k, m*i_loc_m);
    h = (i_glob_0+m<=n?m:n-i_glob_0);
    for(i_glob = i_glob_0; i_glob < i_glob_0 + h; i_glob++ ){
      for(j_loc = 0; j_loc < n; j_loc++){
        j_glob = j_loc;
        i_loc = g2l(n, m, p, k, i_glob);
        a[i_loc*n+j_loc] = initElement(n, s, i_glob+1, j_glob+1);
      }
    }
  }
}

int readArray(FILE *fp, double* a, int n){
  for(int i = 0; i < n; i++){
    if(fscanf(fp, "%lf", a+i) != 1)return 1;
  }
  return 0;
}

int printArray(double* a, int cols, int rows, int printed_rows, int max_print){
  if(printed_rows >= max_print) return 0;
  int p_cols = (cols > max_print?max_print:cols);
  int p_rows = (printed_rows+rows > max_print?max_print-printed_rows:rows);
  for(int i = 0; i < p_rows; i++){
    for(int j = 0; j < p_cols; j++){
      printf(" %10.3e", a[i*cols+j]);
    }
    printf("\n");
  }
  return p_rows;
}


int readMatrix(double* a, int n, int m, int p, int k, const char* name,
  double* buff/* размера блочной строки */, MPI_Comm comm){
    int main_k = 0; FILE *fp = nullptr; int err = 0;
    if(k == main_k){
      fp = fopen(name, "r"); if(fp == nullptr) err = 1;
    }
    MPI_Bcast(&err, 1, MPI_INT, main_k, comm);
    if(err) return err;
    memset(buff, 0, n*m*sizeof(double));
    int b, max_b = (n+m-1)/m;
    for(b = 0; b < max_b; b++){
      int owner = b%p;
      int rows = (b*m+m<=n?m:n-b*m);
      int b_loc = b/p;
      if(k == main_k){
        err += readArray(fp, buff, n*rows);
        if(owner == main_k){
          memcpy(a+b_loc*n*m, buff, n*rows*sizeof(double));
        }
        else{
          MPI_Send(buff, n*rows, MPI_DOUBLE, owner, 0/*tag*/, comm);

        }
      }
      else{
        if(owner == k){
          MPI_Status st;
          MPI_Recv(a+b_loc*n*m, n*rows, MPI_DOUBLE, main_k, 0/*tag*/, comm, &st);
        }

      }
    }
    if(k == main_k){
      fclose(fp); fp = nullptr;
    }
    MPI_Bcast(&err, 1, MPI_INT, main_k, comm);
    return err;
}

void printMatrix(double *a, int n, int m, int p, int k,
  double* buff/*n*m buff*/, int max_print, MPI_Comm comm){
  int main_k = 0;
  int b, max_b = (n+m-1)/m;
  int printed_rows = 0;
  for(b = 0; b < max_b; b++){
    int owner = b%p;
    int rows = (m<n-b*m?m:n-b*m);
    int l_loc = b/p;
    if(k == main_k){
      if(owner == main_k)
        printed_rows += printArray(a+l_loc*m*n, n, rows, printed_rows, max_print);
      else{
        MPI_Status st;
        MPI_Recv(buff, n*rows, MPI_DOUBLE, owner, 0/*tags*/, comm, &st);
        printed_rows += printArray(buff, n, rows, printed_rows, max_print);
      }
    }
    else{
      if(k == owner){
        MPI_Send(a+l_loc*n*m, n*rows, MPI_DOUBLE, main_k, 0/*tags*/, comm);
      }
    }
  }
}


void printVector(double *x, int n, int m, int p, int k,
  double* buff/* m buff */, int max_print, MPI_Comm comm){
  int main_k = 0;
  int b, max_b = (n+m-1)/m;
  int printed_rows = 0;
  for(b = 0; b < max_b; b++){
    int owner = b%p;
    int rows = (m<n-b*m?m:n-b*m);
    int l_loc = b/p;
    if(k == main_k){
      if(owner == main_k)
        printed_rows += printArray(x+l_loc*m, 1, rows, printed_rows, max_print);
      else{
        MPI_Status st;
        MPI_Recv(buff, rows, MPI_DOUBLE, owner, 0/*tags*/, comm, &st);
        printed_rows += printArray(buff, 1, rows, printed_rows, max_print);
      }
    }
    else{
      if(k == owner){
        MPI_Send(x+l_loc*m, rows, MPI_DOUBLE, main_k, 0/*tags*/, comm);
      }
    }
  }
}


void initVector(double* a, double* b, int n, int m, int k, int p){
  int i, j;
  int rows = get_rows(n, m, p, k);

  for(int ii = 0; ii < rows; ii++){
    int lrows = (ii*m+m < n?m:n-ii*m);
    for(i = 0; i < lrows; i++){
      b[ii*m+i] = 0;
      for(j = 0; j < (n+1)/2; j++){
        b[ii*m+i] += a[ii*m*n + i*n + 2 * j];
      }
    }
  }
}
