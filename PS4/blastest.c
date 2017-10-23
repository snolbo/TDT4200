/* -*- Mode:C; Coding:us-ascii-unix; fill-column:132 -*- */
/* ****************************************************************************************************************************** */
/**
   @file      blas2C.c
   @author    Mitch Richling <https://www.mitchr.me/>
   @Copyright Copyright 1997 by Mitch Richling.  All rights reserved.
   @brief     Demonstrate several cblas (level 1) functions. @EOL
   @Keywords  blas cblas C fortran numerical linear algebra vector matrix gemv ger
   @Std       C89

   This is a simple program intended to illistrate how to make use of #gemv and #ger blas routines (as implimented in the cblas).

*/

/* ------------------------------------------------------------------------------------------------------------------------------ */

#include <cblas.h>
#include <stdlib.h>             /* Standard Lib    ISOC  */
#include <math.h>               /* Math stuff      ISOC  */
#include <stdio.h>              /* I/O lib         ISOC  */

void print_matrix(double* mat, int size, int col_size){
  int col = 0;
  for(int i = 0; i < size; i++){
    if(col == col_size){
      col = 0;
      printf("\n");
    }
    printf("%f ", mat[i]);
    col++;
  }
  printf("\n");
}

void print_vector(double *vector, int size){
  for(int i = 0; i < size; i++){
    printf("%f   ", vector[i]);
  }
  printf("\n");
}

int main(int argc, char **argv) {
  double a[5*5] = {  1, 2, 3, 4, 5,
                     6, 7, 8, 9,10,
                    11,12,13,14,15,
                    16,17,18,19,20,
                    21, 22, 23, 24, 25
                  };
  double x[5]  = {1,1,1,1,1};
  double z[5] = {2, 2, 2, 2, 2};
  int lenX = 5;
  int lenY = 5;
  int lda = lenY;
  int alpha = 1;
  int beta = 0;
  int incX = 1;
  int incY = 1;

  double y[lenY];

  // printMatrix(CblasRowMajor, 4, 5, a, 8, 3, NULL, NULL, NULL, NULL, NULL, "              a = ");
  // printVector(5, x, 8, 3, NULL, NULL, NULL, "              x = ");

           /* row_order      transform     lenY lenX alpha  a  lda  X  incX  beta  Y, incY */
  cblas_dgemv(CblasColMajor, CblasNoTrans, lenY,   lenX,   alpha,     a,   lda, x, incX,    beta,    y, incY);
  // printVector(4, y, 8, 3, NULL, NULL, NULL, "       y<-1.0*a*xT+0.0*y= ");
  // print_matrix(a, 20, 5);
  print_vector(y, lenY);
  cblas_daxpy(lenX, alpha, y, incX, a, lenX);
  print_matrix(a, 25, 5);



for( int i = 0; i < 10; ++i){
  printf("%d\n", i);
}


  return 0;
} /* end func main */
