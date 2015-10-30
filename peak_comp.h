#include <stdlib.h>
#include <stdio.h>

#include "gnuplot_i.h"

#define S32K   32768
#define NSPECT 100
#define BIG_NUMBER 99999999999999

//forward declarations
void compareSpectra();
void computeLinearBackground();
void plotSpectra();

int expHist[NSPECT][S32K],simHist[NSPECT][S32K];
char str[256];
double chisq;
double expInt, simInt;//integral of the simulated and experimental data over the specified channel range
double scaleFactor[NSPECT];//factor to scale a given simulated sprectrum by
double scaledSimHist[NSPECT][S32K];
int binsSkipped;
int numBinsUsed;
int i,j;
//parameters used by background addition algorithm
long double m_sum,s_sum,ss_sum,ms_sum,mi_sum,si_sum,i_sum,ii_sum,sum1;//sums for determinants
long double detA[NSPECT], detAi[NSPECT][3];//determinants for Cramer's rule soln
long double bgA[NSPECT],bgB[NSPECT];//linear background parameters
//gnuplot interface
gnuplot_ctrl *handle ;
