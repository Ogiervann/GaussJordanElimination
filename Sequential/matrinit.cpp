#include "matrinit.h"

void setZero(double *a, int n){
  int i;
  for(i = 0; i < n; i++){
    a[i] = 0;
  }
}

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
  int i;
  for(i = 0; i < n; i++){
    pos[i] = i;
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

int readMatrix(double* a, int n, char* filename){
  int i, j, res = 0;
  FILE * f = 0;

  f = fopen(filename, "r");

  if(f == 0){
    return 1;
  }


  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      res += fscanf(f, "%lf", &a[i*n+j]);
    }
  }
  fclose(f);
  if(res != n*n){
    return 1;
  }
  return 0;
}

int initMatrix(double* a, int n, int s, char* filename){
  int i, j, res = 0;
  if(s != 0){
    for(i = 0; i < n; i++){
      for(j =0; j < n; j++){
        a[i*n+j] = initElement(n, s, i+1, j+1);
      }
    }
  }
  else{
    res = readMatrix(a, n, filename);
  }
  return res;
}

void initVector(double* a, double* b, int n){
  int i, j;

  for(i = 0; i < n; i++){
    for(j = 0; j < (n+1)/2; j++){
      b[i] += a[i*n + 2 * j];
    }
  }
}


void printMatrix(double* a, int l, int n, int r){
  int i, j;

  for(i = 0; i < (l>r?r:l); i++){
    for(j = 0; j < (n>r?r:n); j++){
      printf(" %10.3e", a[i*n + j]);
    }
    printf("\n");
  }
}
