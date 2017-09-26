#ifndef RPS_MPI_H
#define RPS_MPI_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "CA.h"
#include "bitmap.h"
#include "RPS.h"

#define WHITE   0
#define ROCK    1
#define PAPER   2
#define SCISSOR 3

const int WINNER_TABLE[4][4] = {  {0, -1, -1, -1},
                                  {0, 0, -1, 1},
                                  {0, 1, 0, -1},
                                  {0, -1, 1, 0} };

// seems to crash for some values when size does not match up with numproc. Think is has to do with definition of local dim in initialize
#define IMG_X 18
#define IMG_Y 18



// Each cell is updated based on neighbors of distance 1
#define BORDER_SIZE 1

// How many iterations?
#define ITERATIONS 70
#endif
