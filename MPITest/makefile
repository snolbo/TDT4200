CFLAGS+= -std=c99

# by uncommenting this line the preprocessor will see #ifdef DEBUG as true
# CFLAGS+= -DDEBUG
objects = arraySendTest.o
objects2 = cartTest.o
SHELL := /bin/bash
# When running "make" in your terminal, the first target will be chosen.
# The syntax of make is basically:
# name : [stuff I need to make this target]

# In this case the dependencies are easy to figure out, so I will not elaborate further
all : $(objects)
	mpicc $(objects) -o arrayTest

# In this target [stuff I need to make this target] is two other targets, namely clean and all
# This command simply runs the clean target, then the all target, thus recompiling everything.
remake : clean all

# We add .PHONY when a target doesn't actually create any output. In this case we just run a shell
# command, removing all object files, i.e files ending on .o
# the * syntax means [anything].o
.PHONY : clean
clean :
	rm -f *.o && rm -f *.gch

# Finally, the test target. Builds the 'all' target, then runs the test script on the output
#.PHONY : test
#test : all
#	./test.sh myProgram
.PHONY : arrayTest
arrayTest : all
	mpirun -np 2 ./arrayTest

cartMake : $(objects2)
	mpicc $(objects2) -o cartTest
