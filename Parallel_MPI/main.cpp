#include "mpi.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>

#include "linearsolve.h"
#define EPS 1e-14



int main(int argc, char* argv[]){
  int n, m, k, p, r, s, task = 17; int main_k = 0; int err = 0, err_loc = 0;
  char* filename = 0;
  double r1 = -1, r2 = -1;
  double t1 = -1, t2 = -1;
  //unsigned int start, end;


  MPI_Comm comm = MPI_COMM_WORLD;
  //MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &p);
  MPI_Comm_rank(comm, &k);


  double *a;
  double *b;
  double *x;
  double *buff;
  int* pos;
  int* mpos;
  double* block1;
  double* block2;
  double* block3;

  if(!((argc == 5 || argc == 6) && sscanf(argv[1], "%d", &n) == 1 &&
    sscanf(argv[2], "%d", &m) == 1 &&
    sscanf(argv[3], "%d", &r) == 1 && sscanf(argv[4], "%d", &s) == 1) ||
    (s == 0 && argc == 5)){
    printf("Usage mpirun -n p %s n m r s (file)\n", argv[0]);
    err_loc = 1;

      MPI_Allreduce(&err_loc, &err, 1, MPI_INT, MPI_MAX, comm);
      if(err) {
        MPI_Finalize();
        return 0;
      }
  }

  if(s == 0){
    filename = argv[5];
  }



//  int rows = get_rows(n, m, p, k);
  int max_rows = get_max_rows(n, m, p);
  a = new double[n*m*max_rows];
  b = new double[m*max_rows];
  x = new double[m*max_rows];
  buff = new double[m*n];
  pos = new int[max_rows];
  mpos = new int[m];
  block1 = new double[m*m];
  block2 = new double[m*m];
  block3 = new double[m*m];



  if(a == NULL || b == NULL || x == NULL || buff == NULL || pos == NULL ||
    mpos == NULL || block1 == NULL || block2 == NULL || block3 == NULL){
    if(a != NULL) delete[] a;
    if(b != NULL) delete[] b;
    if(x != NULL) delete[] x;
    if(buff != NULL) delete[] buff;
    if(pos != NULL) delete[] pos;
    if(mpos != NULL) delete[] mpos;
    if(block1 != NULL) delete[] block1;
    if(block2 != NULL) delete[] block2;
    if(block3 != NULL) delete[] block3;
    err_loc = 1;
  }



  MPI_Allreduce(&err_loc, &err, 1, MPI_INT, MPI_MAX, comm);


    if(err) {
      if(a != NULL) delete[] a;
      if(b != NULL) delete[] b;
      if(x != NULL) delete[] x;
      if(buff != NULL) delete[] buff;
      if(pos != NULL) delete[] pos;
      if(mpos != NULL) delete[] mpos;
      if(block1 != NULL) delete[] block1;
      if(block2 != NULL) delete[] block2;
      if(block3 != NULL) delete[] block3;
      MPI_Finalize();
      return 0;
    }



  memset(a, 0, m*max_rows*n*sizeof(double));
  memset(b, 0, m*max_rows*sizeof(double));
  memset(x, 0, m*max_rows*sizeof(double));
  memset(buff, 0, n*m*sizeof(double));
  memset(pos, 0, max_rows*sizeof(int));
  memset(mpos, 0, m*sizeof(int));
  memset(block1, 0, m*m*sizeof(double));
  memset(block2, 0, m*m*sizeof(double));
  memset(block3, 0, m*m*sizeof(double));


  setPos(mpos, m);
  setPos_mpi(pos, n, m, p, k);


  //printMatrix(a, n, m, p, k, buff, r, comm);


  if(s != 0){
    initMatrix(a, n, m, p, k, s);
  }
  else{
    err_loc = readMatrix(a, n, m, p, k, filename, buff, comm);
  }

  MPI_Allreduce(&err_loc, &err, 1, MPI_INT, MPI_MAX, comm);


  if(err) {

      delete[] a;
      delete[] b;
      delete[] x;
      delete[] buff;
      delete[] pos;
      delete[] mpos;
      delete[] block1;
      delete[] block2;
      delete[] block3;

    MPI_Finalize();
    return 0;
  }


  initVector(a, b, n, m, k, p);


  printMatrix(a, n, m, p, k, buff, r, comm);
  if(k == main_k)printf("\n");
  //printVector(b, n, m, p, k, buff, r, comm);
  double normA = mpi_normMatr(a, n, m, k, p, comm);

  t1 = get_full_time();
  err_loc = linearSolve(a, b, x, buff, pos, block1, block2, block3, mpos,
    normA * EPS, n, m, k, p, comm);
  t1 = get_full_time() - t1;
  MPI_Allreduce(&err_loc, &err, 1, MPI_INT, MPI_MAX, comm);

  if(err_loc == 0){
    if(s != 0){
      initMatrix(a, n, m, p, k, s);
    }
    else{
      err_loc = readMatrix(a, n, m, p, k, filename, buff, comm);
    }
    initVector(a, b, n, m, k, p);

/*
printMatrix(a, n, m, p, k, buff, r, comm);
if(k == main_k)printf("\n");
printVector(b, n, m, p, k, buff, r, comm);

if(k == main_k)printf("\n");
*/
    printVector(x, n, m, p, k, buff, r, comm);
    t2 = get_full_time();
    r1 = findr1(a, b, x, n, m, k, p, comm);
    r2 = findr2(x, n, m, k, p, comm);
    t2 = get_full_time() - t2;
/*
    if(k == main_k)printf("\n");
    printVector(b, n, m, p, k, buff, r, comm);
*/
  }
  if(k == main_k){

    printf (
      "%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d "
      "M = %d P = %d\n",
      argv[0], task, r1, r2, t1, t2, s, n, m, p);

  }
  if(a != NULL) delete[] a;
  if(b != NULL) delete[] b;
  if(x != NULL) delete[] x;
  if(buff != NULL) delete[] buff;
  if(pos != NULL) delete[] pos;
  if(mpos != NULL) delete[] mpos;
  if(block1 != NULL) delete[] block1;
  if(block2 != NULL) delete[] block2;
  if(block3 != NULL) delete[] block3;
  MPI_Finalize();
  return 0;
}
