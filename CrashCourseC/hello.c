// preprocessor commands
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"


// compiling from terminal: gcc <filename>.c (optional) -o <desiredExecutableName>
// <desiredExecutableName> ( to run .exe)

// segmentation fault: tried to access memory that was not yours.. Shame on you

//prototypes, compilers must know this. Implementation irrelevant
double addTwo(int a, double b);
void assign(int* aPtr);
void fillArray(int* aPtr, int len);

// main function
int main(int argc, char* argv[]){
  printf("Hello, world!\n");

  //printing different types
  int anInt = 4;
  double aDoub = 3.2;
  char aChar = 'c';
  printf("An int %d, a double %lf, and a char %c\n", anInt, aDoub, aChar);

  //if statements
  if(0==0){
    printf("Math still works\n");
  }
  else if( 1== 1){
    printf("Math kinds works?\n");
  }
  else{
    printf("Universe is broken \n");
  }

  // for loops
  printf("Watch me count\n");
  for(int i = 0; i < 10; i++){
    printf("%d\n", i);
  }

  // variables and function calls
  int ay = 3;
  double bee = 5.2;
  double see = addTwo(ay, bee);
  printf("%d+%lf=%lf\n", ay, bee, see);



  //command line arguments. arguments passed when running compiled program gets put in char array
  for(int i = 0; i < argc; i++){
    printf("Command line argument %d: %s\n",i, argv[i] );
  }


// strings and static arrays. C dont have strings
  char str[6];
  str[0] = 'h';
  str[1] = 'e';
  str[2] = 'l';
  str[3] = 'l';
  str[4] = 'o';
  str[5] = 0;
  printf("The string is %s\n", str);
  printf("The pointer value is %p\n", str);



  // pointers
  int* ayPtr = &ay;
  printf("The pointer of ayPtr is %p and ay is %d\n",ayPtr, ay );
  ay = 2;
  printf("The pointer of ayPtr is %p and ay is %d\n",ayPtr, ay );
  *ayPtr = 5;
  printf("The pointer of ayPtr is %p and ay is %d\n",ayPtr, ay );
  assign(ayPtr);
  printf("The pointer of ayPtr is %p and ay is %d\n",ayPtr, ay );


  // pointers and functions
  int dee = 2;
  assign(&dee);
  printf("Dee is %d\n",dee );

  // arrays and functions
  int length = 3;
  int someInts[length];
  fillArray(someInts, length); // arrays are really passed with the pointer to first element
  printf("The array constains: ");
  for(int i = 0; i < length; i++){
    printf("%d ", someInts );
  }
  printf("\n");


  // dynamic arraysint
  // int* anArray = (int*) malloc(5*sizeof(int)); // elements have values from previous memory
  int* anArray = (int*) calloc(5, sizeof(int)); // allocates all elements to zero
  // fillArray(anArray, 5);
  printf("The array contains ");
  for(int i = 0; i < length; i++){
    printf("%d ", anArray[i]);
  }
  printf("\n" );
  free(anArray);

  // invalid memory access
  int* anotherArray = (int*) malloc(50*sizeof(int));
  fillArray(anotherArray, 50);
  for(int i = 0; 1==1; i++){
    printf("index %d contains %d\n",i, anotherArray[i] );
  }
  free(anotherArray);



} // end main



// function
double addTwo(int a, double b){
  return a+b;
}

// function and pointer dereferencing
void assign(int* aPtr){
  *aPtr = 10;
}

// function and arrays. aPtr is array
void fillArray(int* aPtr, int len){
  for(int i = 0; i < len; i++){
    aPtr[i] = i;
  }
}
