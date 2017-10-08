#include "RPS.h"
#include <time.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>


int NUM_THREADS;

pthread_barrier_t iteration_barrier;
pthread_t* threads;
pthread_attr_t attr;

cell* petri_A;
cell* petri_B;

typedef struct{
  int len;
  int rest;
} LOADBALANCEDATA;

LOADBALANCEDATA load_balance_data;



void swap_petris();
void *iterate_image_thread(void *spesification);
void initialize_petri();


int main(int argc, char** argv){


  if(argc!=2) {
    puts("Input number of threads");
    return 0;
  }
  NUM_THREADS = strtol(argv[1], NULL, 0); // OK???????????????????????????????????+
  threads = calloc(NUM_THREADS, sizeof(pthread_t));  // OK???????????



  srand(time(NULL));

  printf("running %d iterations\n",ITERATIONS);

  // ALLOCATE MEMORY
  petri_A = calloc(IMG_X*IMG_Y, sizeof(cell));
  petri_B = calloc(IMG_X*IMG_Y, sizeof(cell));


  // CALE CALUCATE AND STORE DATA FOR LOAD BALANCING AND WORK DISTRIBUTION
  int rows_per_thread = (IMG_Y) / NUM_THREADS;
  load_balance_data.len = rows_per_thread;
  load_balance_data.rest = IMG_Y - rows_per_thread*NUM_THREADS; // holds number lines extra that needs to be procssed.
  printf("Rows per thread to process: %d\n", rows_per_thread);
  printf("Extra rows dumped to last thread: %d\n", load_balance_data.rest);

  // Initialize pertri with random colors
  initialize_petri();

  pthread_barrier_init(&iteration_barrier, NULL, NUM_THREADS);
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for(int i = 0; i < NUM_THREADS; i++){
    pthread_create(&threads[i], &attr, iterate_image_thread, (void *)i);
  }

  // Wait or all threads to complete
  for(int i=0;i<NUM_THREADS;i++) {
    pthread_join(threads[i], NULL);
  }
  printf("Joined all threads\n");

  pthread_attr_destroy(&attr);
  pthread_barrier_destroy(&iteration_barrier);


  char filename[] = "RPS_pthread";
  make_bmp(petri_A, filename);


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
  int id = (int)tid;
  int y_start = id*load_balance_data.len;
  int y_end = y_start + load_balance_data.len;

  // IGNORE OUTER BORDERS OF IMAGE
  if(id == 0){ // Dont want to calculate from top row
    y_start++;
  }
  if(id == NUM_THREADS -1){ // DDont want to calculate from button row
    y_end-=2;
    y_end+= load_balance_data.rest; // Last thread take extra lines, Is this ok or is it not good load balancing?????
  }


  // // DEBUG PRINT
  // printf("Thread: %d, ystart row %d yend %d, len: %d\n", id, y_start, y_end, len);
  // pthread_barrier_wait(&iteration_barrier);

  for(int ii = 0; ii < ITERATIONS; ii++){
    // IMAGE ITERATION
    int seed = rand();
    seed = (seed * 0x5DEECE66DL + 0xBL) & 0xFFFFFFFFFFFFL;
    for(int xx = 1; xx < IMG_X - 2; xx++){
      for(int yy = y_start; yy < y_end; yy++){
        petri_B[TRANS(xx,yy)] = next_cell(xx, yy, petri_A, (seed % 8) + 8*(seed < 8));
      }
    }
    // SYNCRONIZE BETWEEN ITERATIONS AND SWAP PETRIS
    pthread_barrier_wait(&iteration_barrier); // BAD IMPLEMENTAION?
    if(id == 0){ // BAD IMPLEMENTATION???????
      swap_petris();
      //PROGRESS PRINT
      if(ii % 100 == 0){
        printf("Thread %d on iteration %d\n", id, ii);
      }
    }
    pthread_barrier_wait(&iteration_barrier); // BAD IMPLEMENTAION?

  } // end for ITERATIONS


}
