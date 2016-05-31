//definitions
#include "topspek.h"
//functions and logic
#include "chisq.c"
#include "fitter.c"
#include "plotter.c"
#include "peak_window.c"
#include "read_parameters.c"
#include "read_data.c"
#include "save_data.c"

int main(int argc, char *argv[])
{
 
  //set up handler to take action upon SIGINT (CTRL-C command)
  struct sigaction sigIntHandler;
  sigIntHandler.sa_handler = sigint_cleanup;
  sigaction(SIGINT, &sigIntHandler, NULL);

  if(argc!=2)
    {
      printf("\ntopspek parameter_file\n----------------------\n\n");
      printf("Compares the .mca spectra designated in the parameter file specified and generates cool statistics.\n\n");
      exit(-1);
    }

  int i,j,k;
  
  //allocate structures
  par *p=(par*)malloc(sizeof(par));
  fitpar *fp=(fitpar*)malloc(sizeof(fitpar));
  data *d=(data*)malloc(sizeof(data));
  fitdata *fd=(fitdata*)malloc(sizeof(fitdata));

  readParFile(argv[1],p); //grab data from the parameter file (see read_config.c)
  
  //check that the number of spectra being compared is fine
  if(p->endSpectrum>=NSPECT)
    {
      printf("ERROR: A spectrum number specified in the parameter file is larger than the maximum value of %i.  Reduce it or increase NSPECT in peak_comp.h and recompile.\n",NSPECT);
      exit(-1);
    }
    

  //read in the data files
  readDataFile(p->expDataName,p->endSpectrum+1,d->expHist);
  readDataFile(p->expDataName,p->endSpectrum+1,d->fittedExpHist);
  for (i=0;i<p->numSimData;i++) 
    {
      readDataFile(p->simDataName[i],p->endSpectrum+1,d->simHist[i]);
      readDataFile(p->simDataName[i],p->endSpectrum+1,d->fittedSimHist[i]);
    }
    
  //generate data for fitting
  for (j=0;j<=p->endSpectrum;j++)
    for (i=0;i<p->numSpectra;i++)
      if(p->spectrum[i]==j)
        {
          if(p->fixBG[i]!=0)
            for (k=0;k<S32K;k++)
              d->fittedExpHist[j][k]=d->expHist[j][k] - p->fixedBGPar[i][0] - p->fixedBGPar[i][1]*k - p->fixedBGPar[i][2]*k*k;
        }
    
  int fi=0;//index for simulated data to be fitted
  for (i=0;i<p->numSimData;i++)
    {
          
      //determine whether simulated data is fitted and read into histograms for fitting as needed
      if(p->simDataFixedAmp[i]==0)
        {
          for (j=0;j<=p->endSpectrum;j++)
            for (k=0;k<S32K;k++)
              d->fittedSimHist[fi][j][k]=d->simHist[i][j][k];
          fi++;
        }
      else if(p->simDataFixedAmp[i]==2)//data scaling is fixed relative to the previous fitted dataset
        {
          if(fi>0)
            for (j=0;j<=p->endSpectrum;j++)
              for (k=0;k<S32K;k++)
                d->fittedSimHist[fi-1][j][k]+=p->simDataFixedAmpValue[i]*d->simHist[i][j][k];//add this data to the data that it is scaled relative to
        }
      else if(p->simDataFixedAmp[i]==1)//data scaling is fixed to a specified value (will not be fitted)
        for (j=0;j<=p->endSpectrum;j++)
          for (k=0;k<S32K;k++)
            d->fittedExpHist[j][k]-=p->simDataFixedAmpValue[i]*d->simHist[i][j][k];
    }
    
  if(p->verbose>=0) printf("Spectra read in...\n");
  
  if(p->peakSearch==1)
    findFittingWindow(p,d);//find peaks and shift fitting windows (see peak_window.c)

  computeBackgroundandScaling(p,d,fp);//get background coefficients and scaling factors (see fitter.c)
   
  applyBackgroundandScaling(p,fp,d,fd);//apply fit parameters to data (see fitter.c)

  compareSpectra(p,d,fd);//generate chisq stats (see chisq.c)
  
  if((p->plotOutput>=1)&&(p->verbose>=0))
    plotSpectra(p,d,fd);//see plotter.c
  
  if(p->saveOutput==1)
    saveSpectra(p,fd);//see save_data.c
  
  //free structures
  free(p);
  free(fp);
  free(d);
  free(fd);

  return 0; //great success
}
