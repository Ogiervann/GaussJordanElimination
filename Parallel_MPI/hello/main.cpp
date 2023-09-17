#include "mpi.h"
#include <cstdio>
#include <cstring>
#define LEN 1234

int main(int argc, char* argv[]){
  int p, k;  const int tag = 0;

  MPI_Comm comm = MPI_COMM_WORLD;
  char buf[LEN];
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &p);
  MPI_Comm_rank(comm, &k);
  printf("wrong\n");
  snprintf(buf, LEN, "Hello from process %d", k);
  if(k != 0){
    MPI_Send(buf, strlen(buf)+1, MPI_CHAR, 0/*dest*/, tag, comm);
  }
  else{
    printf("%s\n", buf);
    for(int i = 1; i < p; i++){
      MPI_Recv(buf, LEN, MPI_CHAR, MPI_ANY_SOURCE, tag, comm, &status);
      printf("%s\n", buf);
    }
  }



  MPI_Finalize();

  return 0;
}
