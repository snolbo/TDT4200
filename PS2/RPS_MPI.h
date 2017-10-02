#ifndef RPS_MPI_H
#define RPS_MPI_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

typedef struct {
  int color;
  int strength;
} cell;

#include "CA.h"
#include "bitmap.h"

#define WHITE   0
#define ROCK    1
#define PAPER   2
#define SCISSOR 3

#define IMG_X 600
#define IMG_Y 600

#define BORDER_SIZE 1

#define ITERATIONS 1000
#endif
