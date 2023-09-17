#include "linearsolve.h"



Pos findMaxBlock(int q, double* a, double* block1, double* block2, int* mpos,
   double epsA, int n, int m, int k, int p, MPI_Comm comm, MPI_Op op_max_Pos,
   MPI_Datatype MY_MPI_POS){
    double norm, maxNorm = -1;
    int i, j, r, max_i = -1, max_j = -1;

    int max_b = n/m;
    int rows = get_rows(n, m, p, k);


    for(i = q/p + (k<q%p); i < rows && l2g(n, m, p, k, i*m)+m<=n; i++){
      for(j = q; j < max_b; j++){
        get_block(i, j, n, m, p, k, a, block1);
        //printf("q: %d i: %d i_glob: %d j: %d | %d \n", q, i, l2g(n, m, p, k, i*m)/m, j, k);
        //printMatrix_old(block1, m, m, m);
        //printf("| %d\n", k);
        r = matrixInverse(m, block1, block2, mpos, epsA);
        if(r == 0){

          //printMatrix_old(block2, m, m, m);
          norm = normMatr(block2, m);
          //printf("i_g: %d j: %d norm: %lf | %d\n", l2g(n,m,p,k,i*m)/m, j, norm, k);
          if(1/norm > 1/maxNorm){
            maxNorm = norm;
            max_i = l2g(n, m, p, k, i*m)/m;
            max_j = j;
          }
        }
      }
    }

    //printf("%d %d| %d\n", g2l(n, m, p, k_q, q*m)/m, max_b, k);

    if(q >= max_b){
      max_i = q;
      max_j = q;
      maxNorm = 1;
    }
    Pos pos_k;
    pos_k.norm = maxNorm;
    pos_k.i = max_i;
    pos_k.j = max_j;
    Pos pos;

    //printf("q: %d max_i: %d max_j: %d maxNorm: %lf | %d\n", q, max_i, max_j, maxNorm, k);

    MPI_Allreduce(&pos_k, &pos, 1, MY_MPI_POS, op_max_Pos, comm);



    return pos;
}


void returnToMonke(double* x, double* b, int *pos, double* block, int n, int m,
   int k, int p, MPI_Comm comm){
  int i, bl = n/m, rows = get_rows(n, m, p, k), i1, k1, k1_0, i1_0;
  MPI_Status st;
  for(i = 0; i < rows; i++){
    //i_0 = l2g(n, m, p, k, i*m)/m;
    //printf("i: %d i_0: %d| %d\n", i, i_0, k);
    i1 = pos[i];
    i1_0 = g2l(n, m, p, k1, i1*m)/m;
    //printf("i1: %d i1_0: %d k1: %d | %d\n", i1, i1_0, k1, k);

    if(k1 == k){
      get_vect_block(i, n, m, p, k, b, block);
      put_vect_block(i1_0, n, m, p, k, x, block);
    }
    else{
      MPI_Send(&k, 1, MPI_INT, k1, 0/*tags*/, comm);
      MPI_Send(pos+i, 1, MPI_INT, k1, 0/*tags*/, comm);
      MPI_Send(b+i*m, m, MPI_DOUBLE, k1, 0/*tags*/, comm);

      MPI_Recv(&k1_0, 1, MPI_INT, MPI_ANY_SOURCE, 0/*tags*/, comm, &st);
      MPI_Recv(&i1, 1, MPI_INT, k1_0, 0/*tags*/, comm, &st);
      i1_0 = g2l(n, m, p, k, i1*m)/m;
      MPI_Recv(x+i1_0*m, m, MPI_DOUBLE, k1_0, 0/*tags*/, comm, &st);
    }

  }



  if(rows > bl && k == get_k(n, m, p, bl)){
    get_vect_block(bl, n, m, p, k, b, block);
    put_vect_block(bl, n, m, p, k, x, block);
  }
}

int linearSolve(double* a, double* b, double* x, double* buff, int* pos,
   double* block1, double* block2, double *block3, int* mpos,
   double epsA, int n, int m, int k, int p, MPI_Comm comm){
  int i, j, q, h1, h2, main_k, err_k = 0, err = 0;
  //printf("Hello1 %d \n",  k);
  ///*
  int i1, j1;
  Pos block_pos;
  MPI_Op op_max_Pos;



  definePosOp(&op_max_Pos);
  MPI_Datatype MY_MPI_POS;
  definePosClass(&MY_MPI_POS);
  //*/
  int max_b = (n+m-1)/m;
  //if(k==0)printf("max_b: %d | %d\n", max_b, k);
  for(q = 0; q < max_b; q++){
    //if(k == 0)printf("Step %d\n", q);
    //MPI_Allreduce(&err_k, &err, 1, MPI_INT, MPI_MAX, comm);

    //printMatrix(a, n, m, p, k, buff, n, comm);
    //if(k == 0)printf("\n");
    //printVector(b, n, m, p, k, buff, n, comm);
    int q_loc = g2l(n, m, p, main_k, q*m)/m;
    ///*
    block_pos = findMaxBlock(q, a, block1, block2, mpos,
      epsA, n, m, k, p, comm, op_max_Pos, MY_MPI_POS);

    if(block_pos.i == -1){
      return -2;
    }
    i1 = block_pos.i;
    j1 = block_pos.j;

    //if(k == 0)printf("q: %d i1: %d j1: %d | %d\n", q, i1, j1, k);

    switchBlockLines(q, i1, a, b, block1, block3, n, m, k, p, comm);
    switchBlockColumns(q, j1, a, pos, block1, block3, n, m, k, p, comm);


    //printMatrix(a, n, m, p, k, buff, n, comm);
    //if(k == 0)printf("\n");
    //printVector(b, n, m, p, k, buff, n, comm);
    //MPI_Allreduce(&err_k, &err, 1, MPI_INT, MPI_MAX, comm);
    //*/
    int h = (q*m+m<n?m:n-q*m);
    //printf("q: %d q_loc: %d main_k: %d h: %d| %d\n", q, q_loc, main_k, h, k);
    if(k == main_k){
      get_block(q_loc, q, n, m, p, k, a, block1);

      err_k = matrixInverse(h, block1, block2, mpos, epsA);

      //printMatrix_old(block2, h, h, n);

      for(j = q+1; err_k == 0 && j < max_b; j++){
        //A_q,j
        h2 = (j * m + m < n?m:n-m*j);
        get_block(q_loc, j, n, m, p, k, a, block1);
        matrixMult(h, h, h2, block2, block1, block3);
        put_block(q_loc, j, n, m, p, k, a, block3);
      }
      //b_q
      if(err_k == 0){
        get_vect_block(q_loc, n, m, p, k, b, block1);
        matrixMult(h, h, 1, block2, block1, block3);
        put_vect_block(q_loc, n, m, p, k, b, block3);
        memcpy(buff, a+q_loc*n*m, n*m*sizeof(double));
      }
    }


    MPI_Bcast(&err_k, 1, MPI_INT, main_k, comm);

    if(err_k != 0){
      return -1;
    }

    MPI_Bcast(buff, n*h, MPI_DOUBLE, main_k, comm);

    int rows = get_rows(n, m, p, k);
    int i_glob;

    //printf("buff | %d\n", k);
    //printMatrix_old(buff, h, n, n);
    //MPI_Allreduce(&err_k, &err, 1, MPI_INT, MPI_MAX, comm);
    //printf("before %d | %d\n", q, k);
    for(j = q+1; j < max_b; j++){
      get_buff_block(j, h, n, m, buff, block3);
      h2 = (j*m+m<n?m:n-j*m);
      //printf("block_buff | %d\n", k);
      //printMatrix_old(block3, h, h2, n);
      //printf("\n");
      for(i = 0; i < rows; i++){
        if(k == main_k && i == q_loc)continue;
        //A_i,j
        i_glob = l2g(n, m, p, k, i*m)/m;
        h1 = (i_glob*m+m < n?m:n-i_glob*m);

        get_block(i, q, n, m, p, k, a, block1);
        matrixMult(h1, h, h2, block1, block3, block2);
        get_block(i, j, n, m, p, k, a, block1);
        matrSubstr_eq(h1, h2, block1, block2);
        put_block(i, j, n, m, p, k, a, block1);
      }
    }

    //printf("after %d | %d\n", q, k);

    if(k == main_k){
      get_vect_block(q_loc, n, m, p, k, b, block3);
    }
    MPI_Bcast(block3, h*h, MPI_DOUBLE, main_k, comm);
    for(i = 0; i < rows; i++){
      if(k == main_k && i==q_loc)continue;
      //b_i
      i_glob = l2g(n, m, p, k, i*m)/m;
      h1 = (i_glob*m+m < n?m:n-i_glob*m);
      get_block(i, q, n, m, p, k, a, block1);
      matrixMult(h1, h, 1, block1, block3, block2);
      get_vect_block(i, n, m, p, k, b, block1);
      vectSubstr_eq(h1, block1, block2);
      put_vect_block(i, n, m, p, k, b, block1);
    }

    MPI_Allreduce(&err_k, &err, 1, MPI_INT, MPI_MAX, comm);

  }


  //MPI_Allreduce(&err_k, &err, 1, MPI_INT, MPI_MAX, comm);

  //int rows_k = get_rows(n, m, p, k);
  //memcpy(x, b, m*rows_k*sizeof(double));
  //printVector(b, n, m, p, k, buff, n, comm);
  //for (i = 0; i < rows_k; i++) {
  //    printf("%d %d | %d\n", l2g(n, m, p, k, i*m)/m, pos[i], k);
  //}

  returnToMonke(x, b, pos, block1, n, m, k, p, comm);
  return 0;
}


double findr1(double* a, double* b, double* x,
  int n, int m, int k, int p, MPI_Comm comm){
  if(n > 11000){
    return 0;
  }
  double tmp = mpi_norml1(b, n, m, k, p, comm), r1 = 0;
  int rows = get_rows(n, m, p, k);

  int max_rows = get_max_rows(n, m, p);
  int src = (k+1+p)%p;
  int dst = (k-1+p)%p;
  for(int l = 0; l < p; l++){
    int lk = (k+l)%p;
    int lk_rows = get_rows(n, m, p, lk);
    for(int lk_i = 0; lk_i < lk_rows; lk_i++){
      for(int i = 0; i < rows; i++){
        //b_i -= A_i,lk * x_lk
        int lk_m = (l2g(n, m, p, lk, lk_i*m)+m<n?m:n-l2g(n, m, p, lk, lk_i*m));
        int i_m = (l2g(n, m, p, k, i*m)+m < n?m:n-l2g(n, m, p, k, i*m));
        //printf("lk: %d i: %d i_m: %d lk_i: %d lk_m: %d | %d\n", lk, i, i_m, lk_i, lk_m, k);
        for(int ii = 0; ii<i_m; ii++){
          double s = 0;
          for(int jj = 0; jj < lk_m; jj++){
            int lk_ig = l2g(n,m,p,lk,lk_i*m+jj);
            //printf("ii: %d jj: %d lk_ig: %d a[i*m*n+ii*n+lk_ig]: %lf x[lk_i*m+jj]: %lf | %d\n", ii, jj, lk_ig, a[i*m*n+ii*n+lk_ig], x[lk_i*m+jj], k);
            s += a[i*m*n+ii*n+lk_ig]*x[lk_i*m+jj];
          }
          //printf("b[i*m+ii]_old: %lf | %d\n", b[i*m+ii], k);
          b[i*m+ii] -= s;
          //printf("b[i*m+ii]_new: %lf | %d\n", b[i*m+ii], k);

        }
      }
    }
    MPI_Status st;
    MPI_Sendrecv_replace(x, max_rows*m, MPI_DOUBLE, dst, 0/*tag*/,
      src, 0/*tag*/, comm, &st);
  }

  r1 = mpi_norml1(b, n, m, k, p, comm);
  return r1/tmp;

}

double findr2(double* x, int n, int m, int k, int p, MPI_Comm comm){
  if(n > 11000){
    return 0;
  }
  int i, j, h;
  double res_k = 0, res = 0;
  int rows = get_rows(n, m, p, k);
  for(i = 0; i < rows; i++){
    int i0 = l2g(n, m, p, k, i*m)/m;
    h = (i0*m+m<n?m:n-i0*m);
    for(j = 0; j < h; j++){
      int i1 = (l2g(n, m, p, k, i*m+j)+1)%2;
      res_k += fabs(x[i*m+j] - i1);
    }
  }
  MPI_Allreduce(&res_k, &res, 1, MPI_DOUBLE, MPI_SUM, comm);
  return res;
}
