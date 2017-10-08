CFLAGS+= -std=c99 -g -O3

# by uncommenting this line the preprocessor will see #ifdef DEBUG as true
# CFLAGS+= -DDEBUG
objects = CA.o bitmap.o RPS.o utils.o

SHELL := /bin/bash

serial :
	gcc -std=c99 -g -O3 CA.h CA.c bitmap.h bitmap.c RPS_serial.c RPS.h -o serialRPS.out

omp :
	gcc -std=c99 -g -O3 -fopenmp CA.h CA.c bitmap.h bitmap.c RPS_omp.c RPS.h -o RPS_omp -lm

pth : # IS IT OK TO USE gnu99 to get pthread_barrier_t???????
	gcc -std=gnu99 -g -O3 -pthread CA.h CA.c bitmap.h bitmap.c RPS_pth.c RPS.h -o RPS_pthread

remake : clean all

# We add .PHONY when a target doesn't actually create any output. In this case we just run a shell
# command, removing all object files, i.e files ending on .o
# the * syntax means [anything].o
.PHONY : clean
clean :
	rm -f *.o && rm -f *.gch && rm -f data/*.bmp && rm -f *.mp4

# Finally, the test target. Builds the 'all' target, then runs the test script on the output
.PHONY : video
video : clean all
	./myProgram && ffmpeg -framerate 60 -i data/CA-%000d.bmp -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4

.PHONY : run
run : clean omp
	./meme.out


omptest:
	gcc -std=c99 -g -O3 -fopenmp OMPtest.c -o omptest.out

pthtest:
	gcc -std=c99 -g -O3 -pthread pthtest.c -o pthtest.out -lm
