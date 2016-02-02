#ifndef DYNAMIC_ARR_H
#define DYNAMIC_ARR_H

#include <stdlib.h>

//Functions for dynamically allocating arrays at runtime (on the heap) that have indices for each dimension
double*** allocateArrayD3(int,int,int);
double** allocateArrayD2(int,int);
int*** allocateArrayI3(int,int,int);
int** allocateArrayI2(int,int);

#endif
