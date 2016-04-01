//definitions
#include "topspek.h"
//functions and logic
#include "chisq.c"
#include "fitter.c"
#include "plotter.c"
#include "peak_window.c"
#include "read_parameters.c"
#include "read_mca.c"
#include "save_data.c"

int main(int argc, char *argv[])
{
 
  //set up handler to take action upon SIGINT (CTRL-C command)
  struct sigaction sigIntHandler;
  sigIntHandler.sa_handler = sigint_cleanup;
  sigaction(SIGINT, &sigIntHandler, NULL);

  if(argc!=2)
    {
      printf("\ntopspek parameter_file\n");
      printf("Compares the .mca spectra designated in the parameter file specified and generates cool statistics.\n\n");
      exit(-1);
    }

  int i,j,k;
  
  //allocate structures
  pc_par *par=(pc_par*)malloc(sizeof(pc_par));
  fit_par *fpar=(fit_par*)malloc(sizeof(fit_par));
  histdata *data=(histdata*)malloc(sizeof(histdata));
  fitteddata *fdata=(fitteddata*)malloc(sizeof(fitteddata));

  readParFile(argv[1],par); //grab data from the parameter file (see read_config.c)
  
  //check that the number of spectra being compared is fine
  if(par->endSpectrum>=NSPECT)
    {
      printf("ERROR: A spectrum number specified in the parameter file is larger than the maximum value of %i.  Reduce it or increase NSPECT in peak_comp.h and recompile.\n",NSPECT);
      exit(-1);
    }
    

  //read in the .mca files
  readMCA(par->expDataName,par->endSpectrum+1,data->expHist);
  readMCA(par->expDataName,par->endSpectrum+1,data->fittedExpHist);
  for (i=0;i<par->numSimData;i++) 
    {
      readMCA(par->simDataName[i],par->endSpectrum+1,data->simHist[i]);
      readMCA(par->simDataName[i],par->endSpectrum+1,data->fittedSimHist[i]);
    }
    
  //generate data for fitting
  for (j=0;j<=par->endSpectrum;j++)
    for (i=0;i<par->numSpectra;i++)
      if(par->spectrum[i]==j)
        {
          if(par->fixBG[i]!=0)
            for (k=0;k<S32K;k++)
              data->fittedExpHist[j][k]=data->expHist[j][k] - par->fixedBGPar[i][0] - par->fixedBGPar[i][1]*k - par->fixedBGPar[i][2]*k*k;
        }
    
  int fi=0;//index for simulated data to be fitted
  for (i=0;i<par->numSimData;i++)
    {
          
      //determine whether simulated data is fitted and read into histograms for fitting as needed
      if(par->simDataFixedAmp[i]==0)
        {
          for (j=0;j<=par->endSpectrum;j++)
            for (k=0;k<S32K;k++)
              data->fittedSimHist[fi][j][k]=data->simHist[i][j][k];
          fi++;
        }
      else if(par->simDataFixedAmp[i]==2)//data scaling is fixed relative to the previous fitted dataset
        {
          if(fi>0)
            for (j=0;j<=par->endSpectrum;j++)
              for (k=0;k<S32K;k++)
                data->fittedSimHist[fi-1][j][k]+=par->simDataFixedAmpValue[i]*data->simHist[i][j][k];//add this data to the data that it is scaled relative to
        }
      else if(par->simDataFixedAmp[i]==1)//data scaling is fixed to a specified value (will not be fitted)
        for (j=0;j<=par->endSpectrum;j++)
          for (k=0;k<S32K;k++)
            data->fittedExpHist[j][k]-=par->simDataFixedAmpValue[i]*data->simHist[i][j][k];
    }
  if(par->verbose>=0) printf("Spectra read in...\n");
  
  if(par->peakSearch==1)
    findFittingWindow(par,data);//find peaks and shift fitting windows (see peak_window.c)

  computeBackgroundandScaling(par,fpar,data);//get background coefficients and scaling factors (see fitter.c)
   
  applyBackgroundandScaling(par,fpar,data,fdata);//apply fit parameters to data (see fitter.c)

  compareSpectra(par,data,fdata);//generate chisq stats (see chisq.c)
  
  if((par->plotOutput>=1)&&(par->verbose>=0))
    plotSpectra(par,data,fdata);//see plotter.c
  
  if(par->saveOutput==1)
    saveSpectra(par,fdata);//see save_data.c
  
  //free structures
  free(par);
  free(fpar);
  free(data);
  free(fdata);

  return 0; //great success
}
