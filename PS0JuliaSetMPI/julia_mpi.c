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
void calculate(complex_t julia_C, int rank, int workingProcessors) {
	int lines2Process = (int) floor(YSIZE / workingProcessors);
	int remainder = YSIZE % workingProcessors;
	if(rank < remainder){ // lower ranks takes extra work, no extra pay though...
		lines2Process++;
	}

	int rowData[XSIZE];
	rank = rank - 1; // rank 0 handles results, and i like 0-indexing sooooo.....
	int row;
	for(int i = 0; i < lines2Process; i++){
		for(int j = 0; j < XSIZE; j++){

			row = rank + i * workingProcessors; // cyclic distribution of workload
			complex_t c, z, temp;
			int iter=0;

	    // find our starting complex number c
			c.real = (x_start + step*row);
			c.imag = (ylower + step*j);

	    // our starting z is c
			z = c;

			while(z.real*z.real + z.imag*z.imag < 4) {
				z = add_complex(square_complex(z), julia_C);
				if(++iter==MAXITER) break;
			}
			rowData[j] = iter;
		}
		MPI_Send(rowData, XSIZE, MPI_INT, 0, 0, MPI_COMM_WORLD); // sending one row of results, BLOCKING!! // save results and send collected or send each row?
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

	int workingProcessors = numproc - 1; // working processors. Process 0 receive result and organize only.

	// Idea is to laod balance using cyclic distribution of rows that need caclulating, since C is row major, and laods rows in cache. Does this matter?

	// lowest ranks gets more work! Obviously not rank 0 though

	if(rank == 0){
		int pixel[XSIZE*YSIZE];
		unsigned char *buffer=calloc(XSIZE*YSIZE*3,1);
		MPI_Status Status;
		int rowData[XSIZE];


		// Receive and store data to pixel
		for(int i = 0; i < YSIZE; i++){
			int senderRank = (i+1) % workingProcessors; // get rank of worker we expect data from
			MPI_Recv(&rowData, XSIZE, MPI_INT, senderRank, 0, MPI_COMM_WORLD, &Status);

			/* create nice image from iteration counts. take care to create it upside
				 down (bmp format) */
			for(int j=0;j<XSIZE;j++) {
				int p=((YSIZE-j-1)*XSIZE+i)*3;
				fancycolour(buffer+p,pixel[PIXEL(i,j)]);
			}
		}

		/* write image to disk */
		savebmp("julia.bmp",buffer,XSIZE,YSIZE);
	} // END RANK 0
	else{
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

		calculate(julia_C, rank, workingProcessors);

	} // END WORKERS




	MPI_Finalize();
}
