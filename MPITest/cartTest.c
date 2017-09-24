#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

MPI_Comm cartcomm;


int main(int argc, char const *argv[]) {
  int ndims = 2; int dims[2]= {0,0};
//  int ndims = 2; int dims = {0, 0}; int periods = {1,1}; int reorder = 1;
  int coords[ndims];

  int rank, numprocs;
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  // get information on how to partition depending on desired structure

  MPI_Dims_create(numprocs, ndims, dims); // # of nodes, #dims 2D, 3D.., dims -> integer array of size ndims specifying the number of nodes  in each dimension
  //printf("dims: %d, %d\n", dims[0], dims[1]);

  // create cartesian partitioning with given structure
  MPI_Comm cartcomm; // create a communicator
  int periods[] = {1, 1}; int reorder = 1;
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &cartcomm);


  // get information about the topology of communicator
  int ndimMesh;
  MPI_Cartdim_get(cartcomm, &ndimMesh);
  //printf("cartesian dimention of mesh: %d\n", ndimMesh );

  MPI_Cart_coords(cartcomm, rank, ndims, coords);

  printf("I am rank: %d. My coordinates in grid: %d, %d\n", rank, coords[0], coords[1] );

  if (rank == 0){
    int direction = 0; // 1 = x direction 0 = y direction
    int displacement = -1; // 1 = shift positive direction, -1 = shift negative direction
    int source; int dest;
    MPI_Cart_shift(cartcomm, direction, displacement, &source, &dest);
    printf("FROM 0: source = %d, dest = %d\n", source, dest);
  }

  MPI_Request reqs[2];
  MPI_Status stats[2];

  reqs[0] = MPI_REQUEST_NULL;
  reqs[1] = MPI_REQUEST_NULL;

  MPI_Waitall(2,reqs,stats);

/*
int MPI_Cart_create ( MPI_Comm comm_old, int ndims, int *dims, int *periods,
                            int reorder, MPI_Comm *comm_cart )
*/

  MPI_Finalize();
  return 0;
}
