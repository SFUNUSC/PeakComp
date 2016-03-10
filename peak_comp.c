#include "peak_comp.h"

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
  printf("\n");

  int i,j,k;
  
  //allocate structures
  pc_par *parameters=(pc_par*)malloc(sizeof(pc_par));
  fit_par *fparameters=(fit_par*)malloc(sizeof(fit_par));
  histdata *data=(histdata*)malloc(sizeof(histdata));
  fitteddata *fdata=(fitteddata*)malloc(sizeof(fitteddata));

  readConfigFile(argv[1],parameters); //grab data from the parameter file (see read_config.c)
  
  //check that the number of spectra being compared is fine
  if(parameters->endSpectrum>=NSPECT)
    {
      printf("ERROR: A spectrum number specified in the parameter file is larger than the maximum value of %i.  Reduce it or increase NSPECT in peak_comp.h and recompile.\n",NSPECT);
      exit(-1);
    }
    

  //read in the .mca files
  readMCA(parameters->expDataName,parameters->endSpectrum+1,data->expHist);
  readMCA(parameters->expDataName,parameters->endSpectrum+1,data->fittedExpHist);
  for (i=0;i<parameters->numSimData;i++) 
    {
      readMCA(parameters->simDataName[i],parameters->endSpectrum+1,data->simHist[i]);
      readMCA(parameters->simDataName[i],parameters->endSpectrum+1,data->fittedSimHist[i]);
    }
    
  //generate data for fitting
  for (j=0;j<=parameters->endSpectrum;j++)
    for (i=0;i<parameters->numSpectra;i++)
      if(parameters->spectrum[i]==j)
        {
          if(parameters->fixBG[i]!=0)
            for (k=0;k<S32K;k++)
              data->fittedExpHist[j][k]=data->expHist[j][k] - parameters->fixedBGPar[i][0] - parameters->fixedBGPar[i][1]*k - parameters->fixedBGPar[i][2]*k*k;
        }
    
  int fi=0;//index for simulated data to be fitted
  for (i=0;i<parameters->numSimData;i++)
    {
          
      //determine whether simulated data is fitted and read into histograms for fitting as needed
      if(parameters->simDataFixedAmp[i]==0)
        {
          for (j=0;j<=parameters->endSpectrum;j++)
            for (k=0;k<S32K;k++)
              data->fittedSimHist[fi][j][k]=data->simHist[i][j][k];
          fi++;
        }
      else if(parameters->simDataFixedAmp[i]==2)//data scaling is fixed relative to the previous fitted dataset
        {
          if(fi>0)
            for (j=0;j<=parameters->endSpectrum;j++)
              for (k=0;k<S32K;k++)
                data->fittedSimHist[fi-1][j][k]+=parameters->simDataFixedAmpValue[i]*data->simHist[i][j][k];//add this data to the data that it is scaled relative to
        }
      else if(parameters->simDataFixedAmp[i]==1)//data scaling is fixed to a specified value (will not be fitted)
        for (j=0;j<=parameters->endSpectrum;j++)
          for (k=0;k<S32K;k++)
            data->fittedExpHist[j][k]-=parameters->simDataFixedAmpValue[i]*data->simHist[i][j][k];
    }
  printf("Spectra read in...\n");
  
  if(parameters->peakSearch==1)
    findFittingWindow(parameters,data);//find peaks and shift fitting windows (see peak_window.c)

  computeBackgroundandScaling(parameters,fparameters,data);//get background coefficients and scaling factors (see fitter.c)
   
  applyBackgroundandScaling(parameters,fparameters,data,fdata);//apply fit parameters to data (see fitter.c)

  compareSpectra(parameters,data,fdata);//generate chisq stats (see chisq.c)
  
  if(parameters->plotOutput>=1)
    plotSpectra(parameters,data,fdata);//see plotter.c
  
  if(parameters->saveOutput==1)
    saveSpectra(parameters,fdata);//see save_data.c
  
  //free structures
  free(parameters);
  free(fparameters);
  free(data);
  free(fdata);

  return 0; //great success
}
