#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <stddef.h>
#include <math.h>

#include "RPS_MPI.h"


void initialize();

void exchange_borders();
void send_and_receive();
void gather_petri();

void create_types();

void fill_border_row_data();
void fill_border_col_data();

void populate_local_col();
void populate_local_row();

void print_local_area();
void print_collected_area();

void allocate_memory_border_exchange();
void free_memory();

void iterate_CA();
void attack_and_update();

int index_of();



int rank; int size;
// The dimensions of the processor grid. Same for every process
int proc_x_dim; int proc_y_dim;
// The location of a process in the process grid. Unique for every process
int my_x_dim; int my_y_dim;
// Variables to store the rank of neighboring processes
int north, south, east, west;
// The dimensions for the process's local area
int local_x_dim; int local_y_dim;

MPI_Comm cart_comm;

// some datatypes, useful for sending data with somewhat less primitive semantics
MPI_Datatype border_row_t; MPI_Datatype border_col_t;
MPI_Datatype local_area_t; MPI_Datatype mpi_cell_t;

cell* local_current; cell* local_next;

// Variables for sending and receiving data
cell* row_to_send;  cell* col_to_send;
cell* row_to_recv;  cell* col_to_recv;

// variable to store the collected cell grid for all processes
cell* collected_area;

// Table used to determine the result of one cell targeting another
const int WINNER_TABLE[4][4] = {  {0, -1, -1, -1},
                                  {0, 0, -1, 1},
                                  {0, 1, 0, -1},
                                  {0, -1, 1, 0} };


int main(int argc, char** argv){
  srand(time(NULL));

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  initialize();

  create_types();
  allocate_memory_border_exchange();

  for(int it = 0; it < ITERATIONS; it++){
    exchange_borders();
    iterate_CA();
  }

  gather_petri();

  if(rank==0){
    cell** img = alloc_img(collected_area, 0);
    make_bmp(img);
  }

  free_memory();
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
    local_x_dim = (IMG_X / proc_x_dim) + BORDER_SIZE*2 +1;
  }
  else{
    local_x_dim = IMG_X / proc_x_dim + BORDER_SIZE*2;
  }

  if( my_y_dim < (IMG_Y % proc_y_dim)){
    local_y_dim = IMG_Y / proc_y_dim + BORDER_SIZE*2 +1;
  }
  else{
    local_y_dim = IMG_Y / proc_y_dim + BORDER_SIZE*2;
  }

  // Free'd free_memory()
  local_current = calloc(local_x_dim * local_y_dim , sizeof(cell));
  local_next = calloc(local_x_dim * local_y_dim , sizeof(cell));

  // Randomly perturb local area
  int rand_color; int rand_str;
  for(int row = BORDER_SIZE; row < local_y_dim - BORDER_SIZE; row++){
    for(int col = BORDER_SIZE; col < local_x_dim - BORDER_SIZE; col++){
      rand_color = rand() % 4;
      rand_str = 1;
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

  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_cell_t);
  MPI_Type_commit(&mpi_cell_t);
  // A message for a local petri-dish
  MPI_Type_contiguous((local_x_dim - 2*BORDER_SIZE) * (local_y_dim - 2*BORDER_SIZE), mpi_cell_t, &local_area_t);
  MPI_Type_commit(&local_area_t);

  MPI_Type_contiguous((local_x_dim-2*BORDER_SIZE), mpi_cell_t, &border_row_t);
  MPI_Type_commit(&border_row_t);
  MPI_Type_contiguous((local_y_dim-2*BORDER_SIZE), mpi_cell_t, &border_col_t);
  MPI_Type_commit(&border_col_t);
}

void allocate_memory_border_exchange(){ // free'd in free_memory
  row_to_send = (cell*) calloc(local_x_dim-2*BORDER_SIZE, sizeof(cell));
  row_to_recv = (cell*) calloc(local_x_dim-2*BORDER_SIZE, sizeof(cell));
  col_to_send = (cell*) calloc(local_y_dim-2*BORDER_SIZE, sizeof(cell));
  col_to_recv = (cell*) calloc(local_y_dim-2*BORDER_SIZE, sizeof(cell));
}

void free_memory(){
  free(local_current);
  free(local_next);

  free(row_to_send);
  free(row_to_recv);

  free(col_to_send);
  free(col_to_recv);

  free(collected_area);

  MPI_Type_free(&border_row_t); MPI_Type_free(&border_col_t);
  MPI_Type_free(&local_area_t); MPI_Type_free(&mpi_cell_t);
}

// Prints the grid cells and received border exhange values for the individual processes
void print_local_area(){
  printf("COLORS:\n");
    for(int row = 0; row < local_y_dim ; row++){
      if(row == BORDER_SIZE || row == local_y_dim-1){
        printf("\n");
      }
      for(int col = 0; col < local_x_dim ; col++){
        if(col == BORDER_SIZE || col == local_x_dim-1){
          printf("  ");
        }
        printf("%d ", local_current[index_of(row,col)].color);
      }

      printf("\n");
    }
    printf("STRENGTHS:\n");
    for(int row = 0; row < local_y_dim ; row++){
      if(row == BORDER_SIZE || row == local_y_dim-1){
        printf("\n");
      }
      for(int col = 0; col < local_x_dim ; col++){
        if(col == BORDER_SIZE || col == local_x_dim-1){
          printf("  ");
        }
        printf("%d ", local_current[index_of(row,col)].strength);
      }
      printf("\n");
    }
    printf("\n");
}

// Prints the collected_area cell values. Must be used after gather_petri
void print_collected_area(){
  printf("\n");
  printf("Color\n");
  for(int i =0; i < IMG_Y; i++){
    if(i % (local_y_dim -2*BORDER_SIZE) == 0){
      printf("\n");
    }
    for(int j = 0; j < IMG_X; j++){
      if(j % (local_x_dim -2*BORDER_SIZE) == 0){
        printf("   ");
      }
      printf("%d ", collected_area[i*IMG_X + j].color);
    }
    printf("\n");

  }
  printf("\n");

/*
  printf("\n");
  printf("Strength\n");
  for(int i =0; i < IMG_Y; i++){
    if(i % (local_y_dim -2*BORDER_SIZE) == 0){
      printf("\n");
    }
    for(int j = 0; j < IMG_X; j++){
      if(j % (local_x_dim -2*BORDER_SIZE) == 0){
        printf("   ");
      }
      printf("%d ", collected_area[i*IMG_X + j].strength);
      // printf("%d ", i*IMG_X +j);
    }
    printf("\n");

  }
  printf("\n");
  */

}
// Exhange border with neighboring processes
void exchange_borders(){

  fill_border_row_data(BORDER_SIZE, row_to_send);
  send_and_receive(north, south, row_to_send, row_to_recv, 1);
  if(south != -1)
    populate_local_row(local_y_dim - 1, row_to_recv);

  fill_border_row_data(local_y_dim - BORDER_SIZE - 1, row_to_send);
  send_and_receive(south, north, row_to_send, row_to_recv, 1);
  if(north != -1)
    populate_local_row(0, row_to_recv);

  fill_border_col_data(BORDER_SIZE, col_to_send);
  send_and_receive(west, east, col_to_send, col_to_recv, 0);
  if(east != -1)
    populate_local_col(local_x_dim - 1, col_to_recv);

  fill_border_col_data(local_x_dim - BORDER_SIZE - 1, col_to_send);
  send_and_receive(east, west, col_to_send, col_to_recv, 0);
  if(west != -1)
    populate_local_col(0, col_to_recv);

}
// Fill row from local data into border_data
void fill_border_row_data(int row_index, cell* border_data){
  for(int jj = BORDER_SIZE; jj < local_x_dim - BORDER_SIZE; jj++){
    border_data[jj - BORDER_SIZE].strength = local_current[index_of(row_index,jj)].strength;
    border_data[jj - BORDER_SIZE].color = local_current[index_of(row_index,jj)].color;
  }
}
// Fill col from local data into border_data
void fill_border_col_data(int col_index, cell* border_data){
  for(int ii = BORDER_SIZE; ii < local_y_dim - BORDER_SIZE; ii++){
    border_data[ii - BORDER_SIZE].strength = local_current[index_of(ii,col_index)].strength;
    border_data[ii - BORDER_SIZE].color = local_current[index_of(ii,col_index)].color;
  }
}
// Takes received data and fills a row in local area
void populate_local_row(int row_index, cell* row_data){
  for(int jj = BORDER_SIZE; jj < local_x_dim - BORDER_SIZE; jj++){
    local_current[index_of(row_index,jj)] = row_data[jj - BORDER_SIZE];
  }
}
// Takes received data and fills a col in local area
void populate_local_col(int col_index, cell* col_data){
  for(int ii = BORDER_SIZE; ii < local_y_dim - BORDER_SIZE; ii++){
    local_current[index_of(ii,col_index)] = col_data[ii - BORDER_SIZE];
  }
}
// Performs the sending and receiving of data in a direction. Used in exchange_borders
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

    MPI_Waitall(2, reqs, stats);
}

// Performs iteration of CA
void iterate_CA(){
  // Iterate local area. Since we have sidebuffer, we dont need to worry about which direction to choose. RIGHT NOW YOU CAN CHOOSE YOURSELF. FIX??
  int random_xdir; int random_ydir; int random_directions[] = {-1,0,1};
  for(int ii = BORDER_SIZE; ii < local_y_dim - BORDER_SIZE; ii++){
    for(int jj = BORDER_SIZE; jj < local_x_dim - BORDER_SIZE; jj++){
      // find direction of cell to attempt to attack
      do{
        random_xdir = random_directions[rand() % 3];
        random_ydir = random_directions[rand() % 3];
      } while( random_xdir == 0 && random_ydir == 0);
      attack_and_update(index_of(ii,jj), index_of(ii + random_ydir, jj + random_xdir));
    }
  }

  // swap local areas
  cell* temp;
  temp = local_current;
  local_current = local_next;
  local_next =  temp;
}

// Performs the needed logic when one cell targets another. Used by iterate_CA
void attack_and_update(int source_cell_index, int target_cell_index){
  if(local_current[source_cell_index].color == WHITE){ // SOURCE CELL IS WHITE. TRANSFORM NEXT STATE COLOR TO TARGET COLOR
    local_next[source_cell_index].color = local_current[target_cell_index].color;
    local_next[source_cell_index].strength = 1;
  }
  else{ // CELL HAS COLOR DIFFERENT THAN WHITE
    int source_cell_result = WINNER_TABLE[local_current[source_cell_index].color][local_current[target_cell_index].color];
    local_current[source_cell_index].strength += source_cell_result;

    if(local_current[source_cell_index].strength <= 0){
      local_next[source_cell_index].color = WHITE;
      local_next[source_cell_index].strength = 1; // UNNECASSARY??????
    }
    else{
      if(local_current[source_cell_index].strength > 6){ // make sure strength is capped at 6
        local_current[source_cell_index].strength = 6;
      }
      local_next[source_cell_index].color = local_current[source_cell_index].color; // transfer color to next state
      local_next[source_cell_index].strength = local_current[source_cell_index].strength; // transfer color to next state

    }
  }
}

// Gather all local data to process 0. This is only function to implement differently for being able to have
// image sizes that is not dependant on number of processes, but since we got informatino we did'nt have to
// I got to lazy to implement this. NOTE Also requires some trixing wiith the commited MPI datatypes!!
void gather_petri(){

  if(rank == 0){
    // Memory to hold data of final collected area of cells. Free'd in free_memory
    collected_area = malloc(IMG_X*IMG_Y*sizeof(cell));

    // Free'd at bottom of if branch
    cell* data_to_recv = malloc((local_x_dim - 2*BORDER_SIZE)*(local_y_dim - 2*BORDER_SIZE)*sizeof(cell));

    int coordinates[2];
    MPI_Status stat;

    int row_length = (local_x_dim - 2*BORDER_SIZE);
    int col_length = (local_y_dim - 2*BORDER_SIZE);

    // Add data from rank 0 to final image
    for(int ii = 0; ii < row_length; ii++){
      for(int jj = 0; jj < col_length; jj++){
        collected_area[ii*IMG_X + jj].strength = local_current[index_of(BORDER_SIZE + ii, BORDER_SIZE + jj)].strength;
        collected_area[ii*IMG_X + jj].color = local_current[index_of(BORDER_SIZE + ii, BORDER_SIZE + jj)].color;
      }
    }

    // Receive and add data from all other ranks, and add them to final image
    for(int p = 1; p < size; p++){
      MPI_Cart_coords(cart_comm, p, 2, coordinates);
      MPI_Recv(data_to_recv, 1, local_area_t, p, p, cart_comm, &stat);

      // variables used for indexing
      int data_row = 0;
      int data_col = 0;
      // Place received data in final image
      for(int ii = coordinates[1]*row_length; ii < coordinates[1]*row_length + row_length; ii++){
        data_col = 0;
        for(int jj = coordinates[0]*col_length; jj < coordinates[0]*col_length + col_length; jj++){
          collected_area[ii*IMG_X + jj].strength = data_to_recv[data_row*row_length + data_col].strength;
          collected_area[ii*IMG_X + jj].color = data_to_recv[data_row*row_length + data_col].color;
          data_col++;
        }
        data_row++;
      }
    }
    free(data_to_recv);
  }
  else{
    // Allogate memory for dta to send. Freed at bot of else branch
    cell* data_to_send = malloc((local_x_dim - 2*BORDER_SIZE) * (local_y_dim - 2*BORDER_SIZE) * sizeof(cell));

    // Add data from local area to buffer that will be sent to process 0
    int row_length = (local_x_dim - 2*BORDER_SIZE);
    for(int ii = 0; ii < local_y_dim - 2*BORDER_SIZE; ii++){
      for(int jj = 0; jj < local_x_dim - 2*BORDER_SIZE; jj++){
        data_to_send[ii*row_length + jj].strength = local_current[index_of(ii + BORDER_SIZE, jj + BORDER_SIZE)].strength;
        data_to_send[ii*row_length + jj].color = local_current[index_of(ii + BORDER_SIZE, jj + BORDER_SIZE)].color;
      }
    }
   MPI_Send(data_to_send, 1, local_area_t, 0, rank, cart_comm);
   free(data_to_send);
  }
}

// Used for accessing correct element in
int index_of(int row, int col){
  return row*local_x_dim + col;
}
