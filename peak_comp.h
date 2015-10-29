#include <stdlib.h>
#include <stdio.h>

#define S32K   32768
#define NSPECT 100
#define BIG_NUMBER 99999999999999

//forward declarations
void compareSpectra();

int expHist[NSPECT][S32K],simHist[NSPECT][S32K];
char str[256];
double chisq;
double expInt, simInt;//integral of the simulated and experimental data over the specified channel range
double scaleFactor;//factor to scale data by
double scaledSimHist[NSPECT][S32K];
int binsSkipped;
int numBinsUsed;
int i,j,k;
//parameters used by background addition algorithm
int origSimHist[NSPECT][S32K];
double minchisq;
int bestBackground, bestBinsSkipped;
