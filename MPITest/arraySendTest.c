//collective.c
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main(int argc, char *argv[])
{
  int rank;
  int numprocs;

  MPI_Init(NULL,NULL);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

  int theArray[3];
  if (rank==0) {
    theArray[0]=10;
    theArray[1]=11;
    theArray[2]=12;
    //data,count,type,source,communicator
    MPI_Send(theArray,3,MPI_INT,1,0,MPI_COMM_WORLD);
    printf("Rank 0 has sent the array\n");
  }
 else {
    //receive arguments irrelevant
    MPI_Status Stat;
    MPI_Recv(theArray,3,MPI_INT, 0, 0, MPI_COMM_WORLD, &Stat);
      printf("Rank %d received array with elements: %d, %d, %d\n",rank,  theArray[0], theArray[1], theArray[2]);
  }

  //free(theArray);

  MPI_Finalize();
}
