#include <signal.h>
#include <stdlib.h>
#include <stdio.h>

#include "gnuplot_i.h"
#include "lin_eq_solver.h"

#define S32K        32768
#define NSPECT      100
#define NSIMDATA    MAX_DIM-2 //2 parameters reserved for background fitting

//forward declarations
void compareSpectra(int);
void computeBackgroundandScaling(int,int);
void plotSpectra();
void saveSpectra();
void sigint_cleanup();

//global variables (umad compsci profs?)
int expHist[NSPECT][S32K],fittedExpHist[NSPECT][S32K],simHist[NSIMDATA][NSPECT][S32K],fittedSimHist[NSIMDATA][NSPECT][S32K];
double scaleFactor[NSIMDATA][NSPECT], fittedScaleFactor[NSIMDATA][NSPECT];//factor to scale a given simulated sprectrum by
double scaledSimHist[NSIMDATA][NSPECT][S32K];
int i,j,k,l;
long double bgA[NSPECT],bgB[NSPECT];//linear background parameters
gnuplot_ctrl *handle;
int plotOpen;//1 if plots are being displayed, 0 otherwise
int spectrum[NSPECT],startCh[NSPECT],endCh[NSPECT],numSpectra,endSpectrum,maxNumCh,numSimData,numFittedSimData;
int addBackground;//0=no,1=constant background
int plotOutput;//0=no,1=yes,2=detailed
int plotStyle;//0=lin-lin, 1=log-lin
int saveOutput;//0=no,1=yes
char expDataName[256],simDataName[NSIMDATA][256],fittedSimDataName[NSIMDATA][256];//filenames for the simulated and experiment data
int simDataFixedAmp[NSIMDATA];//bool specifying whether amplitude of each set of simulated data is fixed
double simDataFixedAmpValue[NSIMDATA];//value at which amplitude is fixed for each set of simulated data
