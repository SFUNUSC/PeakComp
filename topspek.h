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
  int commonScaling;//0=disabled,1=all spectra will have the same scaling
  int verbose;//0=normal,-1=only output chisq value, no plots or other stats
  char expDataName[256],simDataName[NSIMDATA][256],fittedSimDataName[NSIMDATA][256];//filenames for the simulated and experiment data
  int simDataFixedAmp[NSIMDATA];//bool specifying whether amplitude of each set of simulated data is fixed,1=fixed scaline,2=relative scaling
  double simDataFixedAmpValue[NSIMDATA][NSPECT];//value at which amplitude is fixed for each set of simulated data and each spectrum
  int simDataCommonScaling[NSIMDATA];//bool specifying whether scaling is common to each spectrum
  double channelScaling;//value to scale all channel values specified in the parameter file by (useful for looking at the same data with different contraction factors)
  int forcePositiveS;//1=scaling factors in fit will be forced to be positive
}par; //parameters for peak comparison (from parameter file)

typedef struct
{
  double scaleFactor[NSIMDATA][NSPECT];//factor to scale a given simulated sprectrum by
  long double bgA[NSPECT],bgB[NSPECT],bgC[NSPECT];//background parameters (y = A + B*x + C*x*x)
}fitpar; //fit parameters

typedef struct
{
  long double m_sum,s_sum[NSIMDATA],ss_sum[NSIMDATA][NSIMDATA],ms_sum[NSIMDATA],
              mi_sum,mii_sum,si_sum[NSIMDATA],sii_sum[NSIMDATA],i_sum,ii_sum,
              iii_sum,iiii_sum,sum1; //sums needed to construct system of equations
}fitsum; //sums used in the fitting routine

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

