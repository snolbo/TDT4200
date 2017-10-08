#include "RPS.h"
#include <time.h>
#include <omp.h>

void swap_petris();

cell* petri_A;
cell* petri_B;
int NUM_THREADS;

int main(int argc, char** argv){


  if(argc!=2) {
    puts("Input number of threads");
    return 0;
  }
  NUM_THREADS = strtol(argv[1], NULL, 0); // OK???????????????????????????????????+


  printf("running %d iterations\n",ITERATIONS);

  srand(time(NULL));
  petri_A = calloc(IMG_X*IMG_Y, sizeof(cell));
  petri_B = calloc(IMG_X*IMG_Y, sizeof(cell));

  int seed = rand();

  // Seed some CAs
  for(int ii = 0; ii < 100; ii++){
    int rx = rand() % (IMG_X - 1);
    int ry = rand() % (IMG_Y - 1);
    int rt = rand() % 4;

    petri_A[TRANS(rx,ry)].color = rt;
    petri_A[TRANS(rx,ry)].strength = 1;
  }


  // Main loop
  for(int ii = 0; ii < ITERATIONS; ii++){
    int seed = rand();
    seed = (seed * 0x5DEECE66DL + 0xBL) & 0xFFFFFFFFFFFFL;

    int chunk = ceil((IMG_X-2)*(IMG_Y-2)/NUM_THREADS); // will chunk size and static cause overlap of iterations????
    // Parallelized iterations of petri
    #pragma omp parallel for shared(petri_A, petri_B) schedule(static, chunk) collapse(2)
    for(int xx = 1; xx < IMG_X - 2; xx++){
      for(int yy = 1; yy < IMG_Y - 2; yy++){
        petri_B[TRANS(xx,yy)] = next_cell(xx, yy, petri_A, (seed % 8) + 8*(seed < 8));
      }
    }

    // PROGRESS PRINT
    if(ii % 100 == 0)
      printf("Iteration %d\n", ii);

    swap_petris();
  }

  char filename[] = "RPS_omp";
  make_bmp(petri_A, filename);

  free(petri_A);
  free(petri_B);
  // WHY DOES VALGRIND SAY THAT I HAVE NOT FREED ALL MEMORY??????


}






void swap_petris(){
  cell* tmp1 = petri_A;
  petri_A = petri_B;
  petri_B = tmp1;
}
