#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <stddef.h>
#include <math.h>

#include "RPS_MPI.h"

void initialize();
void initialize_petri();
void exchange_borders();
void iterate_CA();
void gather_petri();
void create_types();

void fill_border_row_data();
void fill_border_col_data();

void populate_local_col();
void populate_local_row();

void print_local_area();

void allocate_memory();
void free_memory();

void attack_and_update();
void send_and_receive();

int index_of();

int rank; int size;
// The dimensions of the processor grid. Same for every process
int proc_x_dim; int proc_y_dim;
// The location of a process in the process grid. Unique for every process
int my_x_dim; int my_y_dim;
int north, south, east, west;
// The dimensions for the process local area
int local_x_dim; int local_y_dim;
int side_buffer_size = 1;
MPI_Comm cart_comm;
// some datatypes, useful for sending data with somewhat less primitive semantics
MPI_Datatype border_row_t; MPI_Datatype border_col_t;
MPI_Datatype local_area_t; MPI_Datatype mpi_cell_t;

cell* local_current; cell* local_next;

cell* row_to_send; cell* col_to_send;
cell* row_to_recv; cell* col_to_recv;

cell* collected_grid;

int main(int argc, char** argv){
  srand(1234);

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  initialize();

  create_types();
  allocate_memory();

  for(int it = 0; it < ITERATIONS; it++){
    exchange_borders();
    iterate_CA();
  }

  MPI_Barrier(cart_comm);

  gather_petri();


  // WTF AM I DOING?
  // cell** bufferpointer;
  // *bufferpointer = collected_grid;

  if(rank==0){
   //make_bmp(bufferpointer, 0);

  }

  free_memory(); // NOTE CRASH SEEMS TO HAPPEN HERE. for numproc = 1, 4, 6, 9 , 12
  MPI_Finalize();
  exit(0);
}

void initialize(){

  int dims[2];
  dims[0] = proc_x_dim;
  dims[1] = proc_y_dim;

  int periods[2]; // we set these to 0 because we are not interested in wrap-around
  periods[0] = 0;
  periods[1] = 0;

  int coords[2];
  coords[0] = my_x_dim;
  coords[1] = my_y_dim;

  MPI_Dims_create(size, 2, dims);
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);
  MPI_Cart_coords(cart_comm, rank, 2, coords);
  MPI_Cart_shift(cart_comm, 1, 1, &north, &south);
  MPI_Cart_shift(cart_comm, 0, 1, &west, &east);

  proc_x_dim = dims[0];
  proc_y_dim = dims[1];

  my_x_dim = coords[0];
  my_y_dim = coords[1];

  if( my_x_dim < (IMG_X % proc_x_dim)  ){
    local_x_dim = (IMG_X / proc_x_dim) + side_buffer_size*2 +1;
  }
  else{
    local_x_dim = IMG_X / proc_x_dim + side_buffer_size*2;
  }

  if( my_y_dim < (IMG_Y % proc_y_dim)){
    local_y_dim = IMG_Y / proc_y_dim + side_buffer_size*2 +1;
  }
  else{
    local_y_dim = IMG_Y / proc_y_dim + side_buffer_size*2;
  }
    // printf("RANK %d procxdim: %d, procydim: %d, IMAGE: %d %d. My coords: %d %d, mydims: %d %d , rests: %d, %d. exta line? %d, %d\n",
    //  rank, proc_x_dim, proc_y_dim, IMG_X, IMG_Y, my_x_dim, my_y_dim, local_x_dim, local_y_dim,
    // (IMG_X % proc_x_dim), (IMG_Y % proc_y_dim), my_x_dim < (IMG_X % proc_x_dim), my_y_dim < (IMG_Y % proc_y_dim) );
    //
    // printf("rank: %d north %d , south %d west %d east %d\n", rank, north, south, west, east);

  // Free'd free_memory()
  local_current = calloc(local_x_dim * local_y_dim , sizeof(cell));
  local_next = calloc(local_x_dim * local_y_dim , sizeof(cell));

  // Randomly perturb local area
  for(int row = side_buffer_size; row < local_y_dim - side_buffer_size; row++){
    for(int col = side_buffer_size; col < local_x_dim - side_buffer_size; col++){
      int rand_color = rand() % 4;
      int rand_str = 1;
      if(rand_color > 0){
        rand_str = 1 + rand() % 6;
      }
      local_current[index_of(row,col)].color = rand_color;
      local_current[index_of(row,col)].strength = rand_str;

    }
  }

}
void create_types(){
  // cell type
  const int    nitems=2;
  int          blocklengths[2] = {1,1};
  MPI_Datatype types[2] = {MPI_INT, MPI_INT};
  MPI_Aint     offsets[2];
  offsets[0] = offsetof(cell, color);
  offsets[1] = offsetof(cell, strength);

  // ARE THESE TYPES UNIQUE TO THE PROCSS? CAUSE THEY ARE MADE DIFFERENTLY AT DIFFERENT PROCSSES

  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_cell_t);
  MPI_Type_commit(&mpi_cell_t);
  // A message for a local petri-dish
  MPI_Type_contiguous((local_x_dim - 2*side_buffer_size) * (local_y_dim - 2*side_buffer_size), mpi_cell_t, &local_area_t);
  MPI_Type_commit(&local_area_t);

  MPI_Type_contiguous((local_x_dim-2*side_buffer_size), mpi_cell_t, &border_row_t);
  MPI_Type_commit(&border_row_t);
  MPI_Type_contiguous((local_y_dim-2*side_buffer_size), mpi_cell_t, &border_col_t);
  MPI_Type_commit(&border_col_t);
}
void allocate_memory(){
  row_to_send = (cell*) calloc(local_x_dim-2*side_buffer_size, sizeof(cell));
  row_to_recv = (cell*) calloc(local_x_dim-2*side_buffer_size, sizeof(cell));
  col_to_send = (cell*) calloc(local_y_dim-2*side_buffer_size, sizeof(cell));
  col_to_recv = (cell*) calloc(local_y_dim-2*side_buffer_size, sizeof(cell));
}
void free_memory(){
  free(local_current);
  free(local_next);

  free(row_to_send);
  free(row_to_recv);

  free(col_to_send);
  free(col_to_recv);

  free(collected_grid);

  MPI_Type_free(&border_row_t); MPI_Type_free(&border_col_t);
  MPI_Type_free(&local_area_t); MPI_Type_free(&mpi_cell_t);

}

void print_local_area(){
  // DEBUG-PRINT
  printf("COLORS:\n");
    for(int row = side_buffer_size; row < local_y_dim - side_buffer_size; row++){
      for(int col = side_buffer_size; col < local_x_dim - side_buffer_size; col++){
        printf("%d ", local_current[index_of(row,col)].color);
      }
      printf("\n");
    }
    printf("STRENGTHS:\n");
    for(int row = side_buffer_size; row < local_y_dim - side_buffer_size; row++){
      for(int col = side_buffer_size; col < local_x_dim - side_buffer_size; col++){
        printf("%d ", local_current[index_of(row,col)].color);
      }
      printf("\n");
    }
}
void exchange_borders(){ // CRASHES WHEN NUMPROC = 6,8,10,12,14, when we have 3*2, 4*2 4*3 or 7*2 dims

  fill_border_row_data(side_buffer_size, row_to_send);
  send_and_receive(north, south, row_to_send, row_to_recv, 1);
  if(south != -1)
    populate_local_row(local_y_dim - 1, row_to_recv);

  fill_border_row_data(local_y_dim - side_buffer_size - 1, row_to_send);
  send_and_receive(south, north, row_to_send, row_to_recv, 1);
  if(north != -1)
    populate_local_row(0, row_to_recv);

  fill_border_col_data(side_buffer_size, col_to_send);
  send_and_receive(west, east, col_to_send, col_to_recv, 0);
  if(east != -1)
    populate_local_col(local_x_dim - 1, col_to_recv);

  fill_border_col_data(local_x_dim - side_buffer_size - 1, col_to_send);
  send_and_receive(east, west, col_to_send, col_to_recv, 0);
  if(west != -1)
    populate_local_col(0, col_to_recv);

}
void fill_border_row_data(int row_index, cell* border_data){
  for(int jj = side_buffer_size; jj < local_x_dim - side_buffer_size; jj++){
    border_data[jj - side_buffer_size].strength = local_current[index_of(row_index,jj)].strength;
    border_data[jj - side_buffer_size].color = local_current[index_of(row_index,jj)].color;

  }
}
void fill_border_col_data(int col_index, cell* border_data){
  for(int ii = side_buffer_size; ii < local_y_dim - side_buffer_size; ii++){
    border_data[ii - side_buffer_size].strength = local_current[index_of(ii,col_index)].strength;
    border_data[ii - side_buffer_size].color = local_current[index_of(ii,col_index)].color;

  }
}
void populate_local_row(int row_index, cell* row_data){
  for(int jj = side_buffer_size; jj < local_x_dim - side_buffer_size; jj++){
    local_current[index_of(row_index,jj)] = row_data[jj - side_buffer_size];
  }
}
void populate_local_col(int col_index, cell* col_data){
  for(int ii = side_buffer_size; ii < local_y_dim - side_buffer_size; ii++){
    local_current[index_of(ii,col_index)] = col_data[ii - side_buffer_size];
  }
}
void send_and_receive(int send_direction, int receive_direction, cell* border_data_send, cell* border_data_recv, int is_row_type){
  MPI_Request reqs[2]; MPI_Status stats[2];
  if(is_row_type){
    if(receive_direction != -1){
      MPI_Irecv(border_data_recv, 1, border_row_t, receive_direction, receive_direction, cart_comm, &reqs[0]);
    }
    else{
      reqs[0] = MPI_REQUEST_NULL;
    }
    if(send_direction != -1){
      MPI_Isend(border_data_send, 1, border_row_t, send_direction, rank, cart_comm, &reqs[1]);
    }
    else{
      reqs[1] = MPI_REQUEST_NULL;
    }
  }
  else{
    if(receive_direction != -1){
      MPI_Irecv(border_data_recv, 1, border_col_t, receive_direction, receive_direction, cart_comm, &reqs[0]);
    }
    else{
      reqs[0] = MPI_REQUEST_NULL;
    }
    if(send_direction != -1){
      MPI_Isend(border_data_send, 1, border_col_t, send_direction, rank, cart_comm, &reqs[1]);
    }
    else{
      reqs[1] = MPI_REQUEST_NULL;
    }
  }
    // printf("RANK: %d\n", rank);
    MPI_Waitall(2, reqs, stats);
}
void iterate_CA(){ // PROBLEM!!!!! SINCE BORDER CELLS ACT INDEPENDENT OF EACH OTHER IN 2 PROCESSES, THEY GET TO ATTACK 2 TIMES!, ONCE IN EACH PROCESS :(

  // Iterate local area. Since we have sidebuffer, we dont need to worry about which direction to choose. RIGHT NOW YOU CAN CHOOSE YOURSELF. FIX??
  int random_xdir; int random_ydir; int random_directions[] = {-1,0,1};
  for(int ii = side_buffer_size; ii < local_y_dim - side_buffer_size; ii++){
    for(int jj = side_buffer_size; jj < local_x_dim - side_buffer_size; jj++){
      // find direction of cell to attempt to attack
      random_xdir = random_directions[rand() % 3];
      random_ydir = random_directions[rand() % 3];
      attack_and_update(index_of(ii,jj), index_of(ii + random_ydir, jj + random_xdir));
    }
  }

  // swap local areas
  cell* temp;
  temp = local_current;
  local_current = local_next;
  local_next =  temp;
}
void attack_and_update(int source_cell_index, int target_cell_index){
  if(local_current[source_cell_index].color == WHITE){ // SOURCE CELL IS WHITE. TRANSFORM NEXT STATE COLOR TO TARGET COLOR
    local_next[source_cell_index].color = local_current[target_cell_index].color;
    local_current[source_cell_index].strength = 1;
  }
  else{ // CELL HAS COLOR DIFFERENT THAN WHITE
    int source_cell_result = WINNER_TABLE[local_current[source_cell_index].color][local_current[target_cell_index].color];
    local_current[source_cell_index].strength += source_cell_result;

    if(local_current[source_cell_index].strength < 0){
      local_next[source_cell_index].color = WHITE;
      local_next[source_cell_index].strength = 1; // UNNECASSARY??????
    }
    else{
      if(local_current[source_cell_index].strength > 6){ // make sure strength is capped at 6
        local_current[source_cell_index].strength = 6;
      }
      local_next[source_cell_index].color = local_current[source_cell_index].color; // transfer color to next state
    }
  }
}

void gather_petri(){ // image size must be wholy divided by number of processors in each direction. I have not fucking yet implemented for processes to send different fucking size local fucking areas
  //TODO: Gather the final petri for process rank 0

  if(rank == 0){
    collected_grid = malloc(IMG_X*IMG_Y*sizeof(cell));
    cell* data_to_recv = malloc((local_x_dim - 2*side_buffer_size)*(local_y_dim - 2*side_buffer_size)*sizeof(cell));

    int coordinates[2];
    MPI_Status stat;
    // printf("receive elements %d\n", local_x_dim - 2*side_buffer_size)*(local_y_dim - 2*side_buffer_size);
    printf("Receiving stuff\n");

    for(int p = 1; p < size; p++){
      MPI_Cart_coords(cart_comm, p, 2, coordinates);
      // MPI_Recv(&recv_value, 1, MPI_INT, p, p, cart_comm, &stat);
      MPI_Recv(data_to_recv, 1, local_area_t, p, p, cart_comm, &stat);
      // printf("Rank number %d is at x= %d, y = %d\n", p, coordinates[0], coordinates[1]);

      int row_length = (local_x_dim - 2*side_buffer_size);
      int col_length = (local_y_dim - 2*side_buffer_size);

      for(int ii = 0; ii < local_y_dim - 2*side_buffer_size; ii++){
        for(int jj = 0; jj < local_x_dim - 2*side_buffer_size; jj++){
          // printf("%d ",coordinates[0]*row_length + jj + IMG_X*coordinates[1]*col_length + IMG_X*ii);
          // printf("%d \n", data_to_recv[ii*row_length + jj].color);
          collected_grid[coordinates[0]*row_length + jj + IMG_X*coordinates[1]*col_length + IMG_X*ii].strength = data_to_recv[ii*row_length + jj].strength;
          collected_grid[coordinates[0]*row_length + jj + IMG_X*coordinates[1]*col_length + IMG_X*ii].color = data_to_recv[ii*row_length + jj].color;
        }
        // printf("\n");
      }
      printf("Received from %d\n", p);
    }
    printf("Received all\n");
    free(data_to_recv);
  }
  else{
    cell* data_to_send = malloc((local_x_dim - 2*side_buffer_size) * (local_y_dim - 2*side_buffer_size) * sizeof(cell));
    // int data_to_send = rank;
    int row_length = (local_x_dim - 2*side_buffer_size);
    for(int ii = 0; ii < local_y_dim - 2*side_buffer_size; ii++){
      for(int jj = 0; jj < local_x_dim - 2*side_buffer_size; jj++){
        data_to_send[ii*row_length + jj].strength = local_current[index_of(ii + side_buffer_size, jj + side_buffer_size)].strength;
        data_to_send[ii*row_length + jj].color = local_current[index_of(ii + side_buffer_size, jj + side_buffer_size)].color;
      }
    }

  //  MPI_Send(&data_to_send, 1, MPI_INT, 0, rank, cart_comm);
   MPI_Send(data_to_send, 1, local_area_t, 0, rank, cart_comm);
    free(data_to_send);
  }

}

int index_of(int row, int col){
  return row*local_x_dim + col;
}
