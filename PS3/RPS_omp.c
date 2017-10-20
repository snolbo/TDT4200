#include "RPS.h"
#include <time.h>
#include <omp.h>

void swap_petris();

cell* petri_A;
cell* petri_B;
int NUM_THREADS;

int main(int argc, char** argv){

  // Handle input data
  if(argc!=2) {
    puts("Input number of threads");
    return 0;
  }
  NUM_THREADS = strtol(argv[1], NULL, 0);

  srand(time(NULL));
  printf("running %d iterations\n",ITERATIONS);

  // Allocate memory for petri dishes
  petri_A = calloc(IMG_X*IMG_Y, sizeof(cell));
  petri_B = calloc(IMG_X*IMG_Y, sizeof(cell));

  // Initialize petri dish
  int seed = rand();
  // Seed some CAs
  for(int ii = 0; ii < 100; ii++){
    int rx = rand() % (IMG_X - 1);
    int ry = rand() % (IMG_Y - 1);
    int rt = rand() % 4;

    petri_A[TRANS(rx,ry)].color = rt;
    petri_A[TRANS(rx,ry)].strength = 1;
  }


    for(int ii = 0; ii < ITERATIONS; ii++){
      seed = (seed * 0x5DEECE66DL + 0xBL) & 0xFFFFFFFFFFFFL;
      // Iterations are parallelized
      #pragma omp parallel for num_threads(NUM_THREADS) schedule(static) collapse(2)
      for(int yy = 1; yy < IMG_Y - 2; yy++){
        int seed = rand();
        for(int xx = 1; xx < IMG_X - 2; xx++){
          petri_B[TRANS(xx,yy)] = next_cell(xx, yy, petri_A, (seed % 8) + 8*(seed < 8));
        }
      }
      //Swap petris between each iteration
      swap_petris();
    }

  // Make visualization
  char filename[] = "RPS_omp";
  make_bmp(petri_A, filename);

  // free memory allocated in main
  free(petri_A);
  free(petri_B);
}

void swap_petris(){
  cell* tmp1 = petri_A;
  petri_A = petri_B;
  petri_B = tmp1;
}
