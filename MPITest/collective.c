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

  int* theArray = (int*) malloc(3*sizeof(int));
  if (rank==0) {
    theArray[0]=10;
    theArray[1]=11;
    theArray[2]=12;
  }
  //data,count,type,source,communicator
  MPI_Bcast(theArray,3,MPI_INT,0,MPI_COMM_WORLD);
  printf("Rank %d has the array %d %d %d\n",
      rank,theArray[0],theArray[1],theArray[2]);
  free(theArray);

  if (rank==0) {
    int* allNums=(int*) malloc(numprocs*sizeof(int));
    //sendData,count,type,receivePtr,countPer,type,destination,communicator
    MPI_Gather(&rank,1,MPI_INT,allNums,1,MPI_INT,0,MPI_COMM_WORLD);
    printf("Got ");
    for (int i =0; i < numprocs; i++)
      printf("%d ",allNums[i]);
    printf("\n");
    free(allNums);
  } else {
    //receive arguments irrelevant
    MPI_Gather(&rank,1,MPI_INT,NULL,0,MPI_INT,0,MPI_COMM_WORLD);
  }

  int sumOfRanks;
  //inData,answer,countPer,type,operation,destination,communicator
  MPI_Reduce(&rank,&sumOfRanks,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
  if (rank==0)
    printf("Sum of the ranks is %d\n",sumOfRanks);


  MPI_Finalize();
}
