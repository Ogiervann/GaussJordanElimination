#include "linearsolve.h"


void returnToRegularMonke(int n, double* res, int* mpos){
  double k;
  int i = 0, s, pos;

  while(i < n){
    if(mpos[i] != i){
      pos = mpos[i];
      for(s = 0; s < n; s++){
        k = res[i*n+s];
        res[i*n+s] = res[pos*n+s];
        res[pos*n+s] = k;
      }
      k = mpos[pos];
      mpos[pos] = mpos[i];
      mpos[i] = k;
    }
    else{
      i++;
    }
  }
}

int matrixInverse(int n, double *a, double *res, int* mpos, double epsA){
  int i, j, p, maxL, maxC;
  double an, normA = normMatr(a, n), maxEl = 0;

  if(normA < epsA){
    return 1;
  }

  setId(res, n);
  setPos(mpos, n);

  for( p = 0; p < n; p++){
    maxEl = a[p*n+p];
    maxL = p;
    maxC = p;
    for(i = p; i < n; i++){
      for(j = p; j<n; j++){
        if(a[i*n+j] > maxEl){
          maxEl = a[i*n+j];
          maxL = i;
          maxC = j;
        }
      }
    }

    if(std::abs(maxEl) < epsA){
      return 1;
    }
    switchLines(p, maxL, n, a, res);
    switchColumns(p, maxC, n, a, mpos);



    an = 1.0/maxEl;

    for(j = 0; j < n; j++){
      if(j >= p+1){
        a[p*n+j] *= an;
      }
      res[p*n+j] *= an;
      for(i = 0; i < n; i++){
        if(i == p) continue;
        if(j >= p+1){
          a[i*n + j] -= a[i*n+p] * a[p*n+j];
        }
        res[i*n+j] -= a[i*n+p] * res[p*n+j];
      }
    }

  }



  returnToRegularMonke(n, res, mpos);
  return 0;
}




void findMaxBlock_byThread(int i0, int n, int m, int k, int l, double* a, double* block1,
  double* block2, int* mpos, maxBlockResults *mbres, double epsA, int t, int p){
  int i, j, r, maxpos = -1;
  double norm, maxNorm = -1;

  for(i = t + p*((i0-t+p-1)/p); i < k; i += p){
    for(j = i0; j < k; j++){
      get_block(i, j, n, m, k, l, a, block1);
      r = matrixInverse(m, block1, block2, mpos, epsA);
      if(r == 0){
        norm = normMatr(block2, m);
        if(norm < maxNorm || maxpos == -1){
          maxNorm = norm;
          maxpos = i*k + j;

        }
      }
    }
  }
  mbres->norm = maxNorm;
  mbres->pos = maxpos;
}

int findMaxBlock(int i0, int n, int m, int k, int l, double* a, double* block1,
  double* block2, int* mpos, maxBlockResults *mbres, double epsA, int t, int p){
    int res = -1;
    double maxNorm = -1;
    findMaxBlock_byThread(i0, n, m, k, l, a, block1, block2, mpos, &mbres[t], epsA, t, p);
    reduce_sum(p);
    for(int i = 0; i < p; i++){
      //if(t == 0) printf("%d: %d %lf\n", i, mbres[i].pos, mbres[i]. norm);
      if(1.0/(mbres[i].norm) > 1.0/maxNorm){
        maxNorm = mbres[i].norm;
        res = mbres[i].pos;
      }
    }
    return res;
}

void returnToMonke(int n, int m,  double* x, double* b, int *pos, double* block,
  int t, int p){
  int i, k = n/m;
  int l = n - k*m;
  for(i = t; i < k; i+=p){
    get_vect_block(i, m, k, l, b, block);
    put_vect_block(pos[i], m, k, l, x, block);
  }
  if( l != 0){
    get_vect_block(k, m, k, l, b, block);
    put_vect_block(k, m, k, l, x, block);
  }
  reduce_sum(p);
}

int linearSolve(int n, int m, double* a, double* b, double* x, int* pos,
   double* block1, double* block2, double *block3, int* mpos, maxBlockResults *mbres, double epsA, int t, int p){
  int i, j, q, k, l, block_pos, i1, j1, r, h1, h2;
  k = n/m;
  l = n - k*m;


  for(q = 0; q < k; q++){
    //if(t == 0){

    //printMatrix(a, n, n, n);
    //printf("\n");
    //printMatrix(b, 1, n, n);
  //}
  //reduce_sum(p);

    block_pos = findMaxBlock(q, n, m, k, l, a, block1, block2, mpos, mbres, epsA, t, p);

    if(block_pos == -1){
      return -2;
    }
    i1 = block_pos/k;
    j1 = block_pos - i1*k;

    //if(t == 0)printf("q: %d i: %d j: %d\n", q, i1, j1);

    switchBlockLines(q, i1, n, m, k, l, a, b, block1, block3, t, p);
    switchBlockColumns(q, j1, n, m, k, l, a, pos, block1, block3, t, p);

    get_block(q, q, n, m, k, l, a, block1);
    matrixInverse(m, block1, block2, mpos, epsA);


    if(q%p == t){
      for(j = q+1; j < k+1; j++){
        //A_q,j
        h2 = (j < k?m:l);
        get_block(q, j, n, m, k, l, a, block1);
        matrixMult(m, m, h2, block2, block1, block3);
        put_block(q, j, n, m, k, l, a, block3);
      }
      //b_q
      get_vect_block(q, m, k, l, b, block1);
      matrixMult(m, m, 1, block2, block1, block3);
      put_vect_block(q, m, k, l, b, block3);
    }

    reduce_sum(p);
    for(j = q+1; j < k+1; j++){
      get_block(q, j, n, m, k, l, a, block3);
      h2 = (j < k?m:l);
      for(i = t; i < k+1; i += p){
        if(i == q)continue;
        //A_i,j
        h1 = (i < k?m:l);

        get_block(i, q, n, m, k, l, a, block1);
        matrixMult(h1, m, h2, block1, block3, block2);
        get_block(i, j, n, m, k, l, a, block1);
        matrSubstr_eq(h1, h2, block1, block2);
        put_block(i, j, n, m, k, l, a, block1);
      }
    }


    get_vect_block(q, m, k, l, b, block3);
    for(i = t; i < k+1; i += p){
      if(i==q)continue;
      //b_i
      h1 = (i < k?m:l);
      get_block(i, q, n, m, k, l, a, block1);
      matrixMult(h1, m, 1, block1, block3, block2);
      get_vect_block(i, m, k, l, b, block1);
      vectSubstr_eq(h1, block1, block2);
      put_vect_block(i, m, k, l, b, block1);
    }
    reduce_sum(p);
  }
  if(l != 0 ){


    get_block(k, k, n, m, k, l, a, block1);
    r = matrixInverse(l, block1, block2, mpos, epsA);

    if(r != 0){
      return -2;
    }

    if(k%p == t){
      //b_k
      get_vect_block(k, m, k, l, b, block1);
      matrixMult(l, l, 1, block2, block1, block3);
      put_vect_block(k, m, k, l, b, block3);
    }

    reduce_sum(p);
    get_vect_block(k, m, k, l, b, block3);

    for(i = t; i < k; i += p){
      //b_i
      get_block(i, k, n, m, k, l, a, block1);
      matrixMult(m, l, 1, block1, block3, block2);
      get_vect_block(i, m, k, l, b, block1);
      vectSubstr_eq(m, block1, block2);
      put_vect_block(i, m, k, l, b, block1);
    }
    reduce_sum(p);
  }
/*
  if(t == 0){
    printMatrix(b, 1, n, n);
    for(i = 0; i < k; i++){
      printf("%d ", pos[i]);
    }
    printf("\n");
  }
*/
  returnToMonke(n, m, x, b, pos, block1, t, p);
  return 0;
}


double findr1(int n, int m, double* a, double* b, double* x, double* block1,
  double* block2, double* block3, int t, int p){
  if(n > 11000){
    return 0;
  }
  double tmp, r1 = 0;
  int i, j, k, l, s, q;
  k = n/m;
  l = n-m*k;

  //printMatrix(b, 1, n, n);

  tmp = norml1(b, n);
  //printf("tmp: %lf\n", tmp);
  //if(t == 0)printMatrix(b, 1, n, n);
  reduce_sum(p);

  for(i = t; i < k+1; i += p){
    for(j = 0; j < k+1; j++){
      get_block(i, j, n, m, k, l, a, block1);
      get_vect_block(j, m, k, l, x, block3);

      s = ( i!=k?m:l);
      q = (j!=k?m:l);

      matrixMult(s, q, 1, block1, block3, block2);
      get_vect_block(i, m, k, l, b, block1);

      vectSubstr_eq(s, block1, block2);

      put_vect_block(i, m, k, l, b, block1);
    }
  }

  reduce_sum(p);
  //if(t == 0)printMatrix(b, 1, n, n);
  r1 = norml1(b, n);
  return r1/tmp;

}

double findr2(int n, double* x){
  if(n > 11000){
    return 0;
  }
  int i;
  double res = 0;
  for(i = 0; i < n/2; i++){
    res += fabs(x[2*i] - 1);
  }
  return res;
}
