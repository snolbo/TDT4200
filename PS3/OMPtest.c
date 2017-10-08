#include <omp.h>
#include <stdio.h>
#include <time.h>

int  a, b, i, tid;
float x;


void loopty_loop(){
  a=2;
  for(int i = 0; i < 10000; i++){
    a+=a +1;
  }
}

int main(int argc, char *argv[]) {


   int   i, n, chunk;
   float a[100], b[100], result;

   /* Some initializations */
   n = 100;
   chunk = 10;
   result = 0.0;
   for (i=0; i < n; i++) {
     a[i] = i * 1.0;
     b[i] = i * 2.0;
     }
int thread_id;
#pragma omp parallel shared(a,b,n) private(thread_id)
{
  thread_id = omp_get_thread_num();
  for(int j = 0; j < 20; j++){
    int seed = rand();
    printf("%d\n", seed);

    #pragma omp parallel for private(i) schedule(static,chunk) reduction(+:result)
    for (i=0; i < n; i++){
      result = result + (a[i] * b[i]);
      loopty_loop();
    }
    #pragma omp barrier
    // printf("Hello from thread %d on iteration %d \n",thread_id, j);
  }

}

}
