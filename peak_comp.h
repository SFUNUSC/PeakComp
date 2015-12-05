#include <stdlib.h>
#include <stdio.h>

#include "gnuplot_i.h"
#include "lin_eq_solver.h"

#define S32K        32768
#define NSPECT      100
#define NSIMDATA    8
#define BIG_NUMBER  99999999999999

//forward declarations
void compareSpectra();
void computeLinearBackground(int);
void plotSpectra();

int expHist[NSPECT][S32K],simHist[NSIMDATA][NSPECT][S32K];
char str[256];
double chisq;
double spectChisq[NSPECT];
double expInt, simInt[NSIMDATA];//integral of the simulated and experimental data over the specified channel range
double scaleFactor[NSIMDATA][NSPECT];//factor to scale a given simulated sprectrum by
double scaledSimHist[NSIMDATA][NSPECT][S32K];
int binsSkipped;
int numBinsUsed;
int i,j,k,l;
//parameters used by background addition algorithm
long double bgA[NSPECT],bgB[NSPECT];//linear background parameters
//gnuplot interface
gnuplot_ctrl *handle;

lin_eq_type linEq;
