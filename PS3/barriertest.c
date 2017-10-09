#include <time.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

pthread_mutex_t mutex;
pthread_cond_t cond_var;
pthread_attr_t attr;

int counter;



#define NUM_THREADS 8
pthread_t threads[NUM_THREADS];




void barrier_initself(){
  counter = 0;
}

void barrier_waitself(){
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



void* hello(void *tid){
  for(int i = 0; i < 3; i++){
    printf("Helllo from %li, iteration %d\n", (long)tid, i);
    barrier_waitself();
  }

}




int main(int argc, char const *argv[]) {


  pthread_mutex_init(&mutex, NULL);
  pthread_cond_init (&cond_var, NULL);

  barrier_initself();


  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);



  for(int i = 0; i < NUM_THREADS; i++){
    pthread_create(&threads[i], &attr, &hello, (void *)i);
  }


    pthread_attr_destroy(&attr);
    pthread_mutex_destroy(&mutex);
    pthread_cond_destroy(&cond_var);


    // Wait or all threads to complete
    for(int i=0;i<NUM_THREADS;i++) {
      pthread_join(threads[i], NULL);
    }
    printf("Joined all threads\n");

  return 0;
}
