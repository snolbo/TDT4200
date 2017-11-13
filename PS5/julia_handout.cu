#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <cuda.h>
#include <math.h>


#define XSIZE 2560
#define YSIZE 2048

#define MAXITER 255
#define PIXEL(i,j) ((i)+(j)*XSIZE)


// Same as PS1
typedef unsigned char uchar;
typedef struct {
	float real;
  float imag;
} complex_t;

// implement these
void calculate_cuda(float x_start, float ylower, float step);
__global__
void julia_kernel(int* pixel_device, float x_start, float ylower, float step);

// utilities
void output_bmp();
double walltime();
void calculate_serial();



float x_start=-2.01;
float x_end=1;
float yupper;
float ylower;
float ycenter=1e-6;
float step;

complex_t julia_num;


int pixel_host[XSIZE*YSIZE];
int pixel[XSIZE*YSIZE];

double walltime() {
    static struct timeval t;
    gettimeofday(&t, NULL);
    return (t.tv_sec + 1e-6 * t.tv_usec);
}


// Set up the cuda memory transfers, launch your kernel and extract the finished image
void calculate_cuda(float x_start, float ylower, float step){
	int* pixel_device;
	cudaMalloc(&pixel_device, XSIZE*YSIZE*sizeof(int));
	size_t threads_per_block_dim = 32;
	// Assumin that XSIZE and YSIZE is dividable by 32
	dim3 gridBlock(XSIZE/threads_per_block_dim, YSIZE/threads_per_block_dim);
	dim3 threadBlock(threads_per_block_dim, threads_per_block_dim);
	julia_kernel<<<gridBlock, threadBlock>>>(pixel_device, x_start, ylower, step);
	cudaMemcpy(pixel_host, pixel_device, XSIZE*YSIZE*sizeof(int), cudaMemcpyDeviceToHost);
  cudaFree(pixel_device);

}


// Implement the kernel responsible for iterating a single pixel
__global__
void julia_kernel(int* pixel_device, float x_start, float ylower, float step){
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	complex_t c, z, temp;
	int iter = 0;
	c.real = (x_start + step * i);
	c.imag = (ylower + step * j);
	z = c;
	while(z.real*z.real + z.imag*z.imag < 4){
		temp.real = z.real*z.real - z.imag*z.imag + c.real;
		temp.imag = 2.0*z.real*z.imag + c.imag;
		z = temp;
		if(++iter==MAXITER) break;
	}
	pixel_device[PIXEL(i,j)] = iter;
}

int main(int argc, char **argv) {

	if(argc==1) {
		puts("Usage: JULIA\n");
		puts("Input real and imaginary part. ex: ./julia 0.0 -0.8");
		return 0;
	}

  julia_num.real = strtod(argv[1], NULL);
  julia_num.imag = strtod(argv[2], NULL);

  /* Calculate the range in the y - axis such that we preserve the aspect ratio */
  step = (x_end - x_start)/XSIZE;
  yupper = ycenter + (step * YSIZE)/2;
  ylower = ycenter - (step * YSIZE)/2;



  printf("Calculating with the serial implementation...\n");
  double start_serial = walltime();
  calculate_serial();
  double end_serial = walltime();
  printf("Computation complete. It took %7.3f ms\n\n\n", end_serial - start_serial);


  printf("Checking GPU(s)\n");

  int n_devices;
  cudaGetDeviceCount(&n_devices);
  printf("Number of CUDA devices: %d\n", n_devices);
  cudaDeviceProp device_prop;
  cudaGetDeviceProperties(&device_prop, 0);
  printf("CUDA device name 1: %s\n" , device_prop.name);

  if((n_devices < 1) || (n_devices > 2)){
    printf("You're either on more than 2 GPUs, or something is broken\n");
    printf("Exiting");
    exit(0);
  }

  printf("Calculating with CUDA...\n");
  double start_gpu = walltime();
  calculate_cuda(x_start, ylower, step);
  double end_gpu = walltime();
  printf("Computation complete. It took %7.3f ms\n", end_gpu - start_gpu);


  output_bmp();

  return 0;
}



//////////////////////////////////////////
//////////////////////////////////////////
//////////////////////////////////////////
////// UTILITIES, ALREADY IMPLEMENTED
complex_t add_complex(complex_t a, complex_t b){
  complex_t temp;
  temp.real = a.real + b.real;
  temp.imag = a.imag + b.imag;
  return temp;
}

complex_t add_real(complex_t a, int b){
  complex_t temp;
  temp.real = a.real + b;
  return temp;
}

complex_t square_complex(complex_t c){
  complex_t temp;
  temp.real = c.real*c.real - (c.imag*c.imag);
  temp.imag = 2*c.imag*c.real;
  return temp;
}


void savebmp(char *name,uchar *buffer,int x,int y);
void fancycolour(uchar *p,int iter);

void output_bmp(){
  unsigned char* img_buffer = (unsigned char*)calloc(XSIZE*YSIZE*3, 1);
  for(int ii = 0; ii < XSIZE; ii++){
    for(int jj = 0; jj < YSIZE; jj++){
      int p=((YSIZE-jj-1)*XSIZE+ii)*3;
      fancycolour(img_buffer+p,pixel_host[PIXEL(ii,jj)]);
    }
  }

  char filename[20] = "julia.bmp";
  savebmp(filename, img_buffer, XSIZE, YSIZE);
  free(img_buffer);
}



/* save 24-bits bmp file, buffer must be in bmp format: upside-down */
void savebmp(char *name,uchar *buffer,int x,int y) {
	FILE *f=fopen(name,"wb");
	if(!f) {
		printf("Error writing image to disk.\n");
		return;
	}
	unsigned int size=x*y*3+54;
	uchar header[54]={'B','M',size&255,(size>>8)&255,(size>>16)&255,size>>24,0,
                    0,0,0,54,0,0,0,40,0,0,0,x&255,x>>8,0,0,y&255,y>>8,0,0,1,0,24,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	fwrite(header,1,54,f);
	fwrite(buffer,1,XSIZE*YSIZE*3,f);
	fclose(f);
}


/* given iteration number, set a colour */
void fancycolour(uchar *p,int iter) {
	if(iter==MAXITER);
	else if(iter<8) { p[0]=128+iter*16; p[1]=p[2]=0; }
	else if(iter<24) { p[0]=255; p[1]=p[2]=(iter-8)*16; }
	else if(iter<160) { p[0]=p[1]=255-(iter-24)*2; p[2]=255; }
	else { p[0]=p[1]=(iter-160)*2; p[2]=255-(iter-160)*2; }
}


void calculate_serial() {
	for(int i=0;i<XSIZE;i++) {
		for(int j=0;j<YSIZE;j++) {

			/* Calculate the number of iterations until divergence for each pixel.
			   If divergence never happens, return MAXITER */
			complex_t c;
      complex_t z;
      complex_t temp;
			int iter=0;

      // find our starting complex number c
			c.real = (x_start + step*i);
			c.imag = (ylower + step*j);

      // our starting z is c
			z = c;

      // iterate until we escape
			while(z.real*z.real + z.imag*z.imag < 4) {
        temp.real = (z.real*z.real) - (z.imag*z.imag);
        temp.imag = 2*z.real*z.imag;

        temp.real += julia_num.real;
        temp.imag += julia_num.imag;

				z = temp;
				if(++iter==MAXITER) break;
			}
			pixel[PIXEL(i,j)]=iter;
		}
	}
}
