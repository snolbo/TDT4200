#ifndef JULIA_H
#define JULIA_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bitmap.h"

#define XSIZE 10000
#define YSIZE 10000


#define MAXITER 255

// note that we are extra careful with preprocessor macros. Adding parenthesises is never the
// wrong choice.
 #define PIXEL(i,j) ((i)+(j)*XSIZE)

// In order to work with complex numbers we need a datatype to hold a real and imaginary part
typedef struct {
	double real;
	double imag;
} complex_t;

// It's probably a good idea to add some functions for working on imaginary numbers
// although this is not necessary
complex_t square_complex(complex_t c);
complex_t add_complex(complex_t a, complex_t b);
complex_t add_real(complex_t a, int b);
complex_t add_imag(complex_t a, int b);
double absolute_value(complex_t a);




#endif
