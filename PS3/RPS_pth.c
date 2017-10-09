#include "RPS.h"
#include <time.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

void *iterate_image_thread(void *tid);
void swap_petris();
void initialize_petri();
void barrier_init();
void barrier_wait();

int NUM_THREADS;

pthread_mutex_t mutex;
pthread_cond_t cond_var;
pthread_t* threads;
pthread_attr_t attr;
int counter;

cell* petri_A;
cell* petri_B;

// Datatype used to store data about thread work distribution
typedef struct{
  int len;
  int rest;
} LOADBALANCEDATA;

LOADBALANCEDATA load_balance_data;



int main(int argc, char** argv){
  // Handle input arguments
  if(argc!=2) {
    puts("Input number of threads");
    return 0;
  }
  NUM_THREADS = strtol(argv[1], NULL, 0);

  // Allocate memory for array to hold threads created
  threads = calloc(NUM_THREADS, sizeof(pthread_t));

  srand(time(NULL));
  printf("running %d iterations\n",ITERATIONS);

  // Allocate memory for petri dishes
  petri_A = calloc(IMG_X*IMG_Y, sizeof(cell));
  petri_B = calloc(IMG_X*IMG_Y, sizeof(cell));


  // CAlculate load distribution data and store it in load_balance_data
  int rows_per_thread = (IMG_Y) / NUM_THREADS;
  load_balance_data.len = rows_per_thread;
  load_balance_data.rest = IMG_Y - rows_per_thread*NUM_THREADS; // These rows are dumped on last thread

  // Initialize pertri with random colors
  initialize_petri();

  // Initialize pthread variables used for barrier
  pthread_mutex_init(&mutex, NULL);
  pthread_cond_init (&cond_var, NULL);


  barrier_init();


  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  // Create threads and pass in function to calculate. i is passed as id used for load distribution
  for(long i = 0; i < NUM_THREADS; i++){
    pthread_create(&threads[i], &attr, iterate_image_thread, (void *)i);
  }

  // Wait for all threads to complete and join them back into main thread
  for(int i=0;i<NUM_THREADS;i++) {
    pthread_join(threads[i], NULL);
  }

  // free pthread variables allocated in main
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&mutex);
  pthread_cond_destroy(&cond_var);

  // Create bmp visualization of result
  char filename[] = "RPS_pthread";
  make_bmp(petri_A, filename);

  // free memory allocated in main
  free(petri_A);
  free(petri_B);
  free(threads);
}


void swap_petris(){
  cell* tmp1 = petri_A;
  petri_A = petri_B;
  petri_B = tmp1;
}

void initialize_petri(){
  int seed = rand();

  // Seed some CAs
  for(int ii = 0; ii < 100; ii++){
    int rx = rand() % (IMG_X - 1);
    int ry = rand() % (IMG_Y - 1);
    int rt = rand() % 4;
    petri_A[TRANS(rx,ry)].color = rt;
    petri_A[TRANS(rx,ry)].strength = 1;
  }
}


void *iterate_image_thread(void *tid){
  // CALCULATING WHAT PARTS OF IMAGE TO PROCESS
  int id = (long)tid;
  int y_start = id*load_balance_data.len;
  int y_end = y_start + load_balance_data.len;
  // IGNORE OUTER BORDERS OF IMAGE
  if(id == 0){ // Dont want to calculate from top row
    y_start++;
  }
  if(id == NUM_THREADS -1){ // Dont want to calculate from button row
    y_end-=2;
    y_end+= load_balance_data.rest; // Last thread take extra lines, Is this ok or is it not good load balancing?????
  }

  for(int ii = 0; ii < ITERATIONS; ii++){
    // IMAGE ITERATION
    int seed = rand();
    seed = (seed * 0x5DEECE66DL + 0xBL) & 0xFFFFFFFFFFFFL;
    for(int yy = y_start; yy < y_end; yy++){
      for(int xx = 1; xx < IMG_X - 2; xx++){
        petri_B[TRANS(xx,yy)] = next_cell(xx, yy, petri_A, (seed % 8) + 8*(seed < 8));
      }
    }
    // SYNCRONIZE BETWEEN ITERATIONS AND SWAP PETRIS
    barrier_wait();
    if(id == 0){ swap_petris();}
    barrier_wait();
    } // end for ITERATIONS
  }

//// Implementation of barrier using mutex
void barrier_init(){
  counter = 0;
}

void barrier_wait(){
  pthread_mutex_lock(&mutex);
  counter++;
  if(counter == NUM_THREADS){
    counter = 0;
    pthread_cond_broadcast(&cond_var);
  }
  else{
    while(pthread_cond_wait(&cond_var, &mutex) != 0);
  }
  pthread_mutex_unlock(&mutex);
}
