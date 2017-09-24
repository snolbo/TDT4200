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

void let_them_fight();
void update_color_from_strengths();
void send_and_receive();

int index_of();

int rank;
int size;

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

cell* local_A; cell* local_B;

cell* local_next; cell* local_current;

cell* row_to_send; cell* col_to_send;
cell* row_to_recv; cell* col_to_recv;

MPI_Request reqs[2]; MPI_Status stats[2];

int main(int argc, char** argv){

  srand(1234);

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  initialize();
  create_types();
  // Create memory for data to send and receive. Free'd in bottom of main
  allocate_memory();

  printf("RANK: %d\n",rank);


  for(int it = 0; it < ITERATIONS; it++){
    exchange_borders();
    iterate_CA();
    if(rank == 0)
      printf("RANK: %d, Iteration %d of %d C \n",rank, it+1, ITERATIONS);
  }

  MPI_Barrier(cart_comm);
  printf("Error happened after for loop has finished\n");
  // NOTE CRASH HAPPENS AFTER THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!! AL PROCSS AND IMAGE COMBINATIONS TESTEED
  // NOTE SOMETIMES I GET DOUBLE FREE OR CORRUPTION ERROR

  if(rank==0){
    // TODO: Write the petri to an image
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
  MPI_Cart_shift(cart_comm, 0, 1, &north, &south);
  MPI_Cart_shift(cart_comm, 1, 1, &west, &east);


  proc_x_dim = dims[0];
  proc_y_dim = dims[1];

  // NOTE THIS IS WRONG, IT WORK IN A ROW MAJOR FASHION, NOT BY CARTESISNA COORDINATES, same with above
  my_x_dim = coords[0];
  my_y_dim = coords[1];
  // Dimentions of local area, including send and receive side buffers
  if( my_y_dim +1 == proc_y_dim){ // NOTE tprev NOTE
    local_x_dim = IMG_X / proc_x_dim + side_buffer_size*2 +1;
  }
  else{
    local_x_dim = IMG_X / proc_x_dim + side_buffer_size*2;
  }

  if( my_x_dim +1  == proc_x_dim){ // NOTE PREV PREV NOTE
    local_y_dim = IMG_Y / proc_y_dim + side_buffer_size*2 +1;
  }
  else{
    local_y_dim = IMG_Y / proc_y_dim + side_buffer_size*2;
  }

  // Free'd in bottom of main
  local_A = calloc(local_x_dim * local_y_dim , sizeof(cell));
  local_B = calloc(local_x_dim * local_y_dim , sizeof(cell));

  if(rank == 8){
    printf("xdim: %d, ydim: %d, IMAGE: %d %d. My coords: %d %d, local dims: %d %d \n", dims[0], dims[1], IMG_X, IMG_Y, my_x_dim, my_y_dim, local_x_dim, local_y_dim);
  }
  if(rank == 6){
    printf("xdim: %d, ydim: %d, IMAGE: %d %d. My coords: %d %d, local dims: %d %d \n", dims[0], dims[1], IMG_X, IMG_Y, my_x_dim, my_y_dim, local_x_dim, local_y_dim);
  }
  if(rank == 0){
    printf("xdim: %d, ydim: %d, IMAGE: %d %d. My coords: %d %d, local dims: %d %d \n", dims[0], dims[1], IMG_X, IMG_Y, my_x_dim, my_y_dim, local_x_dim, local_y_dim);
  }
  if(rank == 2){
    printf("xdim: %d, ydim: %d, IMAGE: %d %d. My coords: %d %d, local dims: %d %d \n", dims[0], dims[1], IMG_X, IMG_Y, my_x_dim, my_y_dim, local_x_dim, local_y_dim);
  }
  if(rank == 4){
    printf("Rank: %d, NORTH: %d, SOUTH: %d, WEST: %d, EAST: %d\n", rank, north, south, west, east);
  }
  if(rank == 5){
    printf("xdim: %d, ydim: %d, IMAGE: %d %d. My coords: %d %d, local dims: %d %d \n", dims[0], dims[1], IMG_X, IMG_Y, my_x_dim, my_y_dim, local_x_dim, local_y_dim);
  }

  // Randomly perturb local area
  for(int row = side_buffer_size; row < local_x_dim - side_buffer_size; row++){
    for(int col = side_buffer_size; col < local_y_dim - side_buffer_size; col++){
      int rand_color = rand() % 4;
      int rand_str = 1;
      if(rand_color > 0){
        rand_str = 1 + rand() % 6;
      }
      local_A[index_of(row,col)].color = rand_color;
      local_A[index_of(row,col)].strength = rand_str;
    }
  }
  local_current = local_A;
  local_next = local_B;
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
  MPI_Type_contiguous((local_x_dim - side_buffer_size) * (local_y_dim - side_buffer_size), mpi_cell_t, &local_area_t);
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
void free_memory(){ // THIS CAUSES CRASHES!!!!!!
  free(local_A); free(local_B); local_A = NULL; local_B = NULL; local_current = NULL; local_next = NULL;
  free(row_to_send); free(col_to_send); free(row_to_recv); free(col_to_recv);
  MPI_Type_free(&border_row_t); MPI_Type_free(&border_col_t);
  MPI_Type_free(&local_area_t); MPI_Type_free(&mpi_cell_t);
}

void print_local_area(){
  // DEBUG-PRINT
  printf("COLORS:\n");
    for(int row = 0; row < local_y_dim; row++){
      for(int col = 0; col < local_x_dim; col++){
        printf("%d ", local_current[index_of(row,col)].color);
      }
      printf("\n");
    }
    printf("STRENGTHS:\n");
    for(int row = 0; row < local_y_dim; row++){
      for(int col = 0; col < local_x_dim; col++){
        printf("%d ", local_current[index_of(row,col)].strength);
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
    // Reuse Request so it dont have to be initialized and destroyed all the time
    reqs[0] = NULL;
    reqs[1] = NULL;
}
void iterate_CA(){ // PROBLEM!!!!! SINCE BORDER CELLS ACT INDEPENDENT OF EACH OTHER IN 2 PROCESSES, THEY GET TO ATTACK 2 TIMES!, ONCE IN EACH PROCESS :(


  // Iterate local area. Since we have sidebuffer, we dont need to worry about which direction to choose. RIGHT NOW YOU CAN CHOOSE YOURSELF. FIX??
  int random_xdir; int random_ydir; int random_directions[] = {-1,0,1};
  for(int ii = side_buffer_size; ii < local_y_dim - side_buffer_size; ii++){
    for(int jj = side_buffer_size; jj < local_x_dim - side_buffer_size; jj++){
      // find direction of cell to attempt to attack
      random_xdir = random_directions[rand() % 3];
      random_ydir = random_directions[rand() % 3];
      // random_xdir = 0;
      // random_ydir = -1;
      let_them_fight(index_of(ii,jj), index_of(ii + random_ydir, jj + random_xdir));
    }
  }

  // TODO: OVERLLOOK AND MAKE SURE THIS IS CORRECT. Written very late at night
  // Iterate row sidebuffers
  for(int jj = side_buffer_size; jj < local_x_dim - side_buffer_size; jj++){
    // must have random_ydir == 1 for top row to attack and -1 for bottom row. else, no effect occurs
    random_ydir = random_directions[rand() % 3];
    // random_xdir = 0;
    // random_ydir = -1;
    if(random_ydir == 1){
      random_xdir = random_directions[rand() % 3];
      // random_xdir = 0;
      // random_ydir = -1;
      let_them_fight(index_of(0,jj), index_of(random_ydir, jj + random_xdir));
    }

    //bottom row
    random_ydir = random_directions[rand() % 3];
    // random_xdir = 0;
    // random_ydir = -1;
    if(random_ydir == -1){
      random_xdir = random_directions[rand() % 3];
      // random_xdir = 0;
      // random_ydir = -1;
      let_them_fight(index_of(local_y_dim -1, jj), index_of(local_y_dim-1 + random_ydir, jj + random_xdir));
    }
  }

  // Iterate column sidebuffers
  for(int ii = side_buffer_size; ii < local_y_dim - side_buffer_size; ii++){
    random_xdir = random_directions[rand() % 3];
    // random_xdir = 0;
    // random_ydir = -1;
    if(random_xdir == 1){
      random_ydir = random_directions[rand() % 3];
      // random_xdir = 0;
      // random_ydir = -1;
      let_them_fight(index_of(ii, 0), index_of(ii +random_ydir, random_xdir));
    }

    // right ow
    random_xdir = random_directions[rand() % 3];
    // random_xdir = 0;
    // random_ydir = -1;
    if(random_xdir == -1){
      random_ydir = random_directions[rand() % 3];
      // random_xdir = 0;
      // random_ydir = -1;
      let_them_fight(index_of(ii, local_x_dim -1), index_of(ii + random_ydir, local_x_dim -1 + random_xdir));
    }
  }

  update_color_from_strengths();



  // swap local areas
  // cell* temp;
  // temp = local_A;
  // local_A = local_next;
  // local_next = temp;
  if(local_current == local_A){
    local_current = local_B;
    local_next = local_A;
  }
  else{
    local_current = local_A;
    local_next = local_B;
  }
}
void let_them_fight(int source_cell_index, int target_cell_index){
  //printf("cell1: %d, cell2 %d\n", local_A[source_cell_index].color, local_A[target_cell_index].color);
  if(local_current[source_cell_index].color == WHITE){ // SOURCE CELL IS WHITE. TRANSFORM NEXT STATE COLOR TO TARGET COLOR
    local_next[source_cell_index].color = local_current[target_cell_index].color;
    local_current[source_cell_index].strength = 1;
  }
  else{ // BOTH CELL HAVE COLOR DIFFERENT THAN WHITE
    local_next[source_cell_index].color = local_current[source_cell_index].color; // transfer color to next state
    // Find result of source attack target
    int source_cell_result = WINNER_TABLE[local_current[source_cell_index].color][local_current[target_cell_index].color];
    // Calculate strenghts from result,  and update current state with new strengths.
    local_current[source_cell_index].strength += source_cell_result;
    local_current[target_cell_index].strength -= source_cell_result;
    if(local_current[source_cell_index].strength > 6){
      local_current[source_cell_index].strength = 6;
    }
    else if(local_current[target_cell_index].strength > 6){
      local_current[target_cell_index].strength = 6;
    }
  }
}
void update_color_from_strengths(){
  for(int ii = side_buffer_size; ii < local_y_dim - side_buffer_size; ii++){
    for(int jj = side_buffer_size; jj < local_x_dim - side_buffer_size; jj++){
      // Handle the color for next state for the source cell
      if(local_current[index_of(ii, jj)].strength <= 0){ // check ifhis cell is defeated this iteration, Then change cell to WHITE for next state
        local_next[index_of(ii, jj)].color = WHITE;
        local_next[index_of(ii, jj)].strength = 1; // UNNECASSARY??????
      }
      else{ // transfer cell to next state as it surviced this round
        local_next[index_of(ii, jj)].strength = local_current[index_of(ii, jj)].strength; // We have already transfered colors
      }

    }
  }

}
void gather_petri(){
  //TODO: Gather the final petri for process rank 0

  if(rank == 0){

  }
  else{
    cell* data_to_send = malloc((local_x_dim - side_buffer_size) * (local_y_dim - side_buffer_size) * sizeof(cell));
    for(int ii = 0; ii < local_y_dim - 2*side_buffer_size; ii++){
      for(int jj = 0; jj < local_x_dim - 2*side_buffer_size; jj++){
        data_to_send[ii*(local_x_dim - side_buffer_size) + jj].strength = local_current[index_of(ii + side_buffer_size, jj + side_buffer_size)].strength;
        data_to_send[ii*(local_x_dim - side_buffer_size) + jj].color = local_current[index_of(ii + side_buffer_size, jj + side_buffer_size)].color;
      }
    }

    // DO I NEED TO MAKE CONTIGUOUS DATATYPE FOR GRID OF LOCAL_AREA_T AS RECV BUFFER TYPE? IS BEST WAY TO COORDINATE THIS BY RECEIVEING FROM ONE process
    // AT A TIME WITH RECV AND THEN PLACING THE GRID IN THE CORRECT PALCE?
    free(data_to_send);
  }


}

int index_of(int row, int col){
  return row*local_x_dim + col;
}
