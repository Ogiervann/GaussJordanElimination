#ifndef USED_CLASSES_H
#define USED_CLASSES_H
#include "mpi.h"


class Pos{
public:
  double norm = -1;
  int i = -1;
  int j = -1;
};


void definePosClass(MPI_Datatype *MY_MPI_POS);
void definePosOp(MPI_Op *op);

#endif
