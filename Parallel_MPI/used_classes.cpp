#include "used_classes.h"


void definePosClass(MPI_Datatype *MY_MPI_POS){

  Pos x;
  MPI_Aint d1, d2, d3;

  d1 = (MPI_Aint)&x.norm - (MPI_Aint)&x;
  d2 = (MPI_Aint)&x.i - (MPI_Aint)&x;
  d3 = (MPI_Aint)&x.j - (MPI_Aint)&x;



  int count[3] = {1, 1, 1};
  MPI_Datatype types[3] = {MPI_DOUBLE, MPI_INT, MPI_INT};
  MPI_Aint disp[3] = {d1, d2, d3};


  MPI_Type_create_struct(3, count, disp, types, MY_MPI_POS);

  MPI_Type_commit(MY_MPI_POS);

}

void max_Pos(void *a, void *b, int* len, MPI_Datatype*){
  Pos* a1 = (Pos*) a;
  Pos* b1 = (Pos*) b;
  for(int i = 0; i < *len; i++){
    if(1.0/a1[i].norm > 1.0/b1[i].norm){
      b1[i].norm = a1[i].norm;
      b1[i].i = a1[i].i;
      b1[i].j = a1[i].j;
    }
  }
}


void definePosOp(MPI_Op *op){
  MPI_Op_create(max_Pos, 0, op);
}
