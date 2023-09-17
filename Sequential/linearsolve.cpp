#include "linearsolve.h"


void returnToRegularMonke(int n, double* res, int* mpos){
  double k;
  int i = 0, t, pos;

  while(i < n){
    if(mpos[i] != i){
      pos = mpos[i];
      for(t = 0; t < n; t++){
        k = res[i*n+t];
        res[i*n+t] = res[pos*n+t];
        res[pos*n+t] = k;
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
  //printf("%lf\n", normA);

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


int findMaxBlock(int i0, int n, int m, int k, int l, double* a, double* block1,
  double* block2, int* mpos, double epsA){
  int i, j, r, res = -1, flag = 0;
  double maxNorm = -1, norm;
  if(k == 1){
    return 0;
  }

  for(i = i0; i < k; i++){
    for(j = i0; j < k; j++){
      get_block(i, j, n, m, k, l, a, block1);
      r = matrixInverse(m, block1, block2, mpos, epsA);
      if(r == 0){
        norm = normMatr(block2, m);
        if(norm < maxNorm || flag == 0){
          maxNorm = norm;
          res = i*k+j;
          flag = 1;
        }
      }
    }
  }
  return res;
}


void returnToMonke(int n, int m,  double* x, double* b, int *pos, double* block){
  int i, k = n/m;
  int l = n - k*m;
  for(i = 0; i < k; i++){
    get_vect_block(i, m, k, l, b, block);
    put_vect_block(pos[i], m, k, l, x, block);
  }
  if( l != 0){
    get_vect_block(k, m, k, l, b, block);
    put_vect_block(k, m, k, l, x, block);
  }

}

int linearSolve(int n, int m, double* a, double* b, double* x, int* pos,
  double* block1, double* block2, double *block3, int* mpos, double epsA){
  int i, j, p, k, l, block_pos, i1, j1, r;
  k = n/m;
  l = n - k*m;


  for(p = 0; p < k; p++){
    block_pos = findMaxBlock(p, n, m, k, l, a, block1, block2, mpos, epsA);

    if(block_pos == -1){
      return 1;
    }
    i1 = block_pos/k;
    j1 = block_pos - i1*k;

    switchBlockLines(p, i1, n, m, k, l, a, b, block1, block3);
    switchBlockColumns(p, j1, n, m, k, l, a, pos, block1, block3);

    get_block(p, p, n, m, k, l, a, block1);
    matrixInverse(m, block1, block2, mpos, epsA);


    for(j = p+1; j < k; j++){
      //A_p,j
      get_block(p, j, n, m, k, l, a, block1);
      matrixMult(m, m, m, block2, block1, block3);
      put_block(p, j, n, m, k, l, a, block3);
    }

    if(l != 0){
      //A_p,k
      get_block(p, k, n, m, k, l, a, block1);
      matrixMult(m, m, l, block2, block1, block3);
      put_block(p, k, n, m, k, l, a, block3);
    }

    //b_p
    get_vect_block(p, m, k, l, b, block1);
    matrixMult(m, m, 1, block2, block1, block3);
    put_vect_block(p, m, k, l, b, block3);

    for(j = p+1; j < k; j++){
      get_block(p, j, n, m, k, l, a, block3);
      for(i = 0; i < k; i++){
        if(i == p)continue;
        //A_i,j
        get_block(i, p, n, m, k, l, a, block1);
        matrixMult(m, m, m, block1, block3, block2);
        get_block(i, j, n, m, k, l, a, block1);
        matrSubstr_eq(m, m, block1, block2);
        put_block(i, j, n, m, k, l, a, block1);
      }
      if(l != 0){
        //A_k,j
        get_block(k, p, n, m, k, l, a, block1);
        matrixMult(l, m, m, block1, block3, block2);
        get_block(k, j, n, m, k, l, a, block1);
        matrSubstr_eq(l, m, block1, block2);
        put_block(k, j, n, m, k, l, a, block1);
      }

    }
    if(l != 0){
      get_block(p, k, n, m, k, l, a, block3);
      for(i = 0; i < k; i++){
        if(i == p)continue;
        //A_i,k
        get_block(i, p, n, m, k, l, a, block1);
        matrixMult(m, m, l, block1, block3, block2);
        get_block(i, k, n, m, k, l, a, block1);
        matrSubstr_eq(m, l, block1, block2);
        put_block(i, k, n, m, k, l, a, block1);
      }
      //A_k,k
      get_block(k, p, n, m, k, l, a, block1);
      matrixMult(l, m, l, block1, block3, block2);
      get_block(k, k, n, m, k, l, a, block1);
      matrSubstr_eq(l, l, block1, block2);
      put_block(k, k, n, m, k, l, a, block1);
    }


    get_vect_block(p, m, k, l, b, block3);

    for(i = 0; i < k; i++){
      if(i==p)continue;
      //b_i
      get_block(i, p, n, m, k, l, a, block1);
      matrixMult(m, m, 1, block1, block3, block2);
      get_vect_block(i, m, k, l, b, block1);
      vectSubstr_eq(m, block1, block2);
      put_vect_block(i, m, k, l, b, block1);
    }
    if(l != 0){
      //b_k
      get_block(k, p, n, m, k, l, a, block1);
      matrixMult(l, m, 1, block1, block3, block2);
      get_vect_block(k, m, k, l, b, block1);
      vectSubstr_eq(l, block1, block2);
      put_vect_block(k, m, k, l, b, block1);
    }
  }
  if(l != 0){

    get_block(k, k, n, m, k, l, a, block1);
    r = matrixInverse(l, block1, block2, mpos, epsA);
    if(r != 0){
      return 1;
    }

    //b_k
    get_vect_block(k, m, k, l, b, block1);
    matrixMult(l, l, 1, block2, block1, block3);
    put_vect_block(k, m, k, l, b, block3);

    for(i = 0; i < k; i++){
      //b_i
      get_block(i, k, n, m, k, l, a, block1);
      matrixMult(m, l, 1, block1, block3, block2);
      get_vect_block(i, m, k, l, b, block1);
      vectSubstr_eq(m, block1, block2);
      put_vect_block(i, m, k, l, b, block1);
    }
  }
  returnToMonke(n, m, x, b, pos, block1);
  return 0;
}


int linearSolve1(int n, int m, double* a, double* b, double* x, int* pos,
  double* block1, double* block2, double *block3, int* mpos, double epsA){
  int i, j, p, k, l, block_pos, i1, j1, r, h1, h2;
  k = n/m;
  l = n - k*m;


  for(p = 0; p < k; p++){
    block_pos = findMaxBlock(p, n, m, k, l, a, block1, block2, mpos, epsA);

    if(block_pos == -1){
      return 1;
    }
    i1 = block_pos/k;
    j1 = block_pos - i1*k;

    switchBlockLines(p, i1, n, m, k, l, a, b, block1, block3);
    switchBlockColumns(p, j1, n, m, k, l, a, pos, block1, block3);

    get_block(p, p, n, m, k, l, a, block1);
    matrixInverse(m, block1, block2, mpos, epsA);


    for(j = p+1; j < k+1; j++){
      //A_p,j
      h2 = (j < k?m:l);

      get_block(p, j, n, m, k, l, a, block1);
      matrixMult(m, m, h2, block2, block1, block3);
      put_block(p, j, n, m, k, l, a, block3);
    }

    //b_p
    get_vect_block(p, m, k, l, b, block1);
    matrixMult(m, m, 1, block2, block1, block3);
    put_vect_block(p, m, k, l, b, block3);

    for(j = p+1; j < k+1; j++){
      get_block(p, j, n, m, k, l, a, block3);
      h2 = (j < k?m:l);
      for(i = 0; i < k+1; i++){
        if(i == p)continue;
        //A_i,j
        h1 = (i < k?m:l);

        get_block(i, p, n, m, k, l, a, block1);
        matrixMult(h1, m, h2, block1, block3, block2);
        get_block(i, j, n, m, k, l, a, block1);
        matrSubstr_eq(h1, h2, block1, block2);
        put_block(i, j, n, m, k, l, a, block1);
      }
    }


    get_vect_block(p, m, k, l, b, block3);

    for(i = 0; i < k; i++){
      if(i==p)continue;
      //b_i
      h1 = (i < k?m:l);
      get_block(i, p, n, m, k, l, a, block1);
      matrixMult(h1, m, 1, block1, block3, block2);
      get_vect_block(i, m, k, l, b, block1);
      vectSubstr_eq(h1, block1, block2);
      put_vect_block(i, m, k, l, b, block1);
    }
  }
  if(l != 0){

    get_block(k, k, n, m, k, l, a, block1);
    r = matrixInverse(l, block1, block2, mpos, epsA);
    if(r != 0){
      return 1;
    }

    //b_k
    get_vect_block(k, m, k, l, b, block1);
    matrixMult(l, l, 1, block2, block1, block3);
    put_vect_block(k, m, k, l, b, block3);

    for(i = 0; i < k; i++){
      //b_i
      get_block(i, k, n, m, k, l, a, block1);
      matrixMult(m, l, 1, block1, block3, block2);
      get_vect_block(i, m, k, l, b, block1);
      vectSubstr_eq(m, block1, block2);
      put_vect_block(i, m, k, l, b, block1);
    }
  }


  returnToMonke(n, m, x, b, pos, block1);

  //printMatrix(b, 1, n, n);
  //printMatrix(x, 1, n, n);

  return 0;
}


double findr1(int n, int m, double* a, double* b, double* x, double* block1, double* block2, double* block3){
  if(n > 11000){
    return 0;
  }
  double tmp, r1 = 0;
  int i, j, k, l, s, t;
  k = n/m;
  l = n-m*k;


  tmp = norml1(b, n);

  //printf("%lf\n", tmp);

  for(i = 0; i < k+1; i++){
    for(j = 0; j < k+1; j++){
      get_block(i, j, n, m, k, l, a, block1);
      get_vect_block(j, m, k, l, x, block3);

      s = ( i!=k?m:l);
      t = (j!=k?m:l);

      matrixMult(s, t, 1, block1, block3, block2);

      get_vect_block(i, m, k, l, b, block1);
      vectSubstr_eq(s, block1, block2);
      put_vect_block(i, m, k, l, b, block1);
    }
  }
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
