#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "julia_mpi.h"
#include "mpi.h"

double x_start=-2.01;
double x_end=1;
double yupper;
double ylower;

double ycenter=1e-6;
double step;

int pixel[XSIZE*YSIZE];


// I suggest you implement these, however you can do fine without them if you'd rather operate
// on your complex number directly.
complex_t square_complex(complex_t c){
	// ???
	complex_t squared;
	squared.real = (pow(c.real, 2) - pow(c.imag, 2));
	squared.imag = 2 * c.real * c.imag;
	return squared;
}

complex_t add_complex(complex_t a, complex_t b){
  // ???
	complex_t sum;
	sum.real = a.real + b.real;
	sum.imag = a.imag + b.imag;
	return sum;
}

complex_t add_real(complex_t a, int b){
  // ???
	complex_t sum;
	sum.real = a.real + b;
	return sum;
}

complex_t add_imag(complex_t a, int b) {
	// ???
	complex_t sum;
	sum.imag = a.imag + b;
	return sum;
}

double absolute_value(complex_t a) {
	return sqrt(pow(a.real, 2) + pow(a.imag, 2));
}

// add julia_c input arg here?
void calculate(complex_t julia_C, double lines2Process, int numproc) {

	int row;
	for(int i = 0; i < lines2Process, i++){
		for(int j = 0; j < XSIZE, j++){
			row = rank + i * numproc;
			complex_t c, z, temp;
			int iter=0;

	    // find our starting complex number c
			c.real = (x_start + step*row);
			c.imag = (ylower + step*j);

	    // our starting z is c
			z = c;

			while(z.real*z.real + z.imag*z.imag < 4) {
				z = add_complex(square_complex(z), julia_C); // ??????????????
				if(++iter==MAXITER) break;
			}
			pixel[PIXEL(i,j)]=iter;

		}
	}
}


int main(int argc,char **argv) {
	if(argc==1) {
		puts("Usage: JULIA\n");
		puts("Input real and imaginary part. ex: ./julia 0.0 -0.8");
		return 0;
	}

	int rank;
	int numproc;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numproc);

	double lines2Process = floor(ceil(YSIZE/numproc));
	int remainder = YSIZE % numproc;

	if(rank < remainder){
		lines2Process++;
	}


	/* Calculate the range in the y-axis such that we preserve the
	   aspect ratio */
	step=(x_end-x_start)/XSIZE;
	yupper=ycenter+(step*YSIZE)/2;
	ylower=ycenter-(step*YSIZE)/2;

  // Unlike the mandelbrot set where C is the coordinate being iterated, the
  // julia C is the same for all points and can be chosed arbitrarily
  complex_t julia_C;

  // Get the command line args
  julia_C.real = strtod(argv[1], NULL);
  julia_C.imag = strtod(argv[2], NULL);

	calculate(julia_C, lines2Process, numproc);

	MPI_Barrier(MPI_COMM_WORLD); // hold so all pixel are set before process 0 save image
  /* create nice image from iteration counts. take care to create it upside
     down (bmp format) */
	if(rank == 0){
		unsigned char *buffer=calloc(XSIZE*YSIZE*3,1);
		for(int i=0;i<XSIZE;i++) {
			for(int j=0;j<YSIZE;j++) {
				int p=((YSIZE-j-1)*XSIZE+i)*3;
				fancycolour(buffer+p,pixel[PIXEL(i,j)]);
			}
		}
		/* write image to disk */
		savebmp("julia.bmp",buffer,XSIZE,YSIZE);
	}
	MPI_Finalize();
}
