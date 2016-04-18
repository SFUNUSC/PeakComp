#ifndef TS_H
#define TS_H

#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "dynamic_arrays.h"
#include "gnuplot_i.h"
#include "lin_eq_solver.h"
#include "peak_find.h"

#define S32K        32768
#define NSPECT      100
#define NSIMDATA    MAX_DIM-3 //3 parameters reserved for background fitting

//structures
typedef struct
{
  int spectrum[NSPECT];//spectrum indices in the .mca file(s) to compare
  int startCh[NSPECT],endCh[NSPECT],fixBG[NSPECT],numSpectra,endSpectrum,maxNumCh,numSimData,numFittedSimData;
  double fixedBGPar[NSPECT][3];
  int addBackground;//0=no,1=lin background,2=quadratic background
  int fitAddBackground[NSPECT];
  int peakSearch;//0=no,1=yes
  int peakSearchWidth;//width of the window to set around the peak found by peak search, in channels 
  int plotOutput;//0=no,1=yes,2=detailed
  int saveOutput;//0=no,1=yes
  int verbose;//0=normal,-1=only output chisq value, no plots or other stats
  char expDataName[256],simDataName[NSIMDATA][256],fittedSimDataName[NSIMDATA][256];//filenames for the simulated and experiment data
  int simDataFixedAmp[NSIMDATA];//bool specifying whether amplitude of each set of simulated data is fixed
  double simDataFixedAmpValue[NSIMDATA];//value at which amplitude is fixed for each set of simulated data
}par; //parameters for peak comparison (from parameter file)

typedef struct
{
  double scaleFactor[NSIMDATA][NSPECT];//factor to scale a given simulated sprectrum by
  long double bgA[NSPECT],bgB[NSPECT],bgC[NSPECT];//background parameters (y = A + B*x + C*x*x)
}fitpar; //fit parameters

typedef struct
{
  int expHist[NSPECT][S32K];
  int fittedExpHist[NSPECT][S32K];
  int simHist[NSIMDATA][NSPECT][S32K];
  int fittedSimHist[NSIMDATA][NSPECT][S32K];
}data;

typedef struct
{
  double scaledSimHist[NSIMDATA][NSPECT][S32K];
  double bgHist[NSPECT][S32K];
}fitdata;

//global variables (umad compsci profs?)
gnuplot_ctrl *handle;
int plotOpen;//1 if plots are being displayed, 0 otherwise

#endif

