#ifndef PEAK_FIND
#define PEAK_FIND

#include <stdlib.h>
#include <stdio.h>

#define BIG_NUMBER      1E10

typedef struct
{
  int windowSize,searchMin,searchMax;
}peak_search_par; //parameters used when finding a peak

typedef struct
{
  double centroid,width;
}peak_fit_par; //fit parameters


//Functions for dynamically allocating arrays at runtime (on the heap) that have indices for each dimension
peak_fit_par findPeak(const peak_search_par*,const double*,const int);

#endif
