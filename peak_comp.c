#include "peak_comp.h"

int main(int argc, char *argv[])
{
 
  //set up handler to take action upon SIGINT (CTRL-C command)
  struct sigaction sigIntHandler;
  sigIntHandler.sa_handler = sigint_cleanup;
  sigaction(SIGINT, &sigIntHandler, NULL);

  if(argc!=2)
    {
      printf("\npeak_comp parameter_file\n");
      printf("Compares the .mca spectra designated in the parameter file specified and generates cool statistics.\n\n");
      exit(-1);
    }
  printf("\n");

  FILE *expData,*simData[NSIMDATA];
  int i,j,k;
  pc_par parameters;
  fit_par fparameters;
  
  //initialize values
  memset(expHist,0,sizeof(expHist));
  memset(fittedExpHist,0,sizeof(fittedExpHist));
  memset(simHist,0,sizeof(simHist));
  memset(fittedSimHist,0,sizeof(fittedSimHist));
  memset(scaledSimHist,0,sizeof(scaledSimHist));

  readConfigFile(argv[1],&parameters); //grab data from the parameter file (see read_config.c)
  
  //check that the number of spectra being compared is fine
  if(parameters.endSpectrum>=NSPECT)
    {
      printf("ERROR: A spectrum number specified in the parameter file is larger than the maximum value of %i.  Reduce it or increase NSPECT in peak_comp.h and recompile.\n",NSPECT);
      exit(-1);
    }

  //read in the .mca files
  if((expData=fopen(parameters.expDataName,"r"))==NULL)
    {
      printf("ERROR: Cannot open the experiment data file %s!\n",parameters.expDataName);
      exit(-1);
    }
  for (i=0;i<parameters.numSimData;i++)  
    if((simData[i]=fopen(parameters.simDataName[i],"r"))==NULL)
      {
        printf("ERROR: Cannot open the simulated data file %s!\n",parameters.simDataName[i]);
        exit(-1);
      }
  for (i=0;i<=parameters.endSpectrum;i++)
    {
      if(fread(expHist[i],S32K*sizeof(int),1,expData)!=1)
        {
          printf("ERROR: Error reading file %s!\n",parameters.expDataName);
          printf("Verify that the format and number of spectra in the file are correct.\n");
          exit(-1);
        }
      for (j=0;j<S32K;j++)
        fittedExpHist[i][j]=expHist[i][j];
    }
  //generate data for fitting
  for (j=0;j<=parameters.endSpectrum;j++)
    for (i=0;i<parameters.numSpectra;i++)
      if(parameters.spectrum[i]==j)
        {
          if(parameters.fixBG[i]!=0)
            for (k=0;k<S32K;k++)
              fittedExpHist[j][k]=expHist[j][k] - parameters.fixedBGPar[i][0] - parameters.fixedBGPar[i][1]*k - parameters.fixedBGPar[i][2]*k*k;
        }
    
  int fi=0;//index for simulated data to be fitted
  for (i=0;i<parameters.numSimData;i++)
    {
      //read all simulated data in
      for (j=0;j<=parameters.endSpectrum;j++)
        if(fread(simHist[i][j],S32K*sizeof(int),1,simData[i])!=1)
          {
            printf("ERROR: Error reading file %s!\n",parameters.simDataName[i]);
            printf("Verify that the format and number of spectra in the file are correct.\n");
            exit(-1);
          }
          
      //determine whether simulated data is fitted and read into histograms for fitting as needed
      if(parameters.simDataFixedAmp[i]==0)
        {
          for (j=0;j<=parameters.endSpectrum;j++)
            for (k=0;k<S32K;k++)
              fittedSimHist[fi][j][k]=simHist[i][j][k];
          fi++;
        }
      else if(parameters.simDataFixedAmp[i]==2)//data scaling is fixed relative to the previous fitted dataset
        {
          if(fi>0)
            for (j=0;j<=parameters.endSpectrum;j++)
              for (k=0;k<S32K;k++)
                fittedSimHist[fi-1][j][k]+=parameters.simDataFixedAmpValue[i]*simHist[i][j][k];//add this data to the data that it is scaled relative to
        }
      else if(parameters.simDataFixedAmp[i]==1)//data scaling is fixed to a specified value (will not be fitted)
        for (j=0;j<=parameters.endSpectrum;j++)
          for (k=0;k<S32K;k++)
            fittedExpHist[j][k]-=parameters.simDataFixedAmpValue[i]*simHist[i][j][k];
    }
  fclose(expData);
  for (i=0;i<parameters.numSimData;i++) 
    fclose(simData[i]);
  printf("Spectra read in...\n");
  
  if(parameters.peakSearch==1)
    findFittingWindow(&parameters);//find peaks and shift fitting windows (see peak_window.c)

  computeBackgroundandScaling(&parameters,&fparameters);//get background coefficients and scaling factors (see fitter.c)
      
  //scale simulated data
  for (i=0;i<parameters.numSimData;i++)
    for (j=0;j<parameters.numSpectra;j++)
      for (k=0;k<S32K;k++)
        scaledSimHist[i][parameters.spectrum[j]][k]=fparameters.scaleFactor[i][parameters.spectrum[j]]*simHist[i][parameters.spectrum[j]][k];
  //generate background data
  for (i=0;i<parameters.numSpectra;i++)
    for (j=0;j<S32K;j++)
      if(parameters.addBackground>=1)
        bgHist[i][j]=fparameters.bgA[parameters.spectrum[i]] + fparameters.bgB[parameters.spectrum[i]]*j + fparameters.bgC[parameters.spectrum[i]]*j*j;

  compareSpectra(&parameters);//generate chisq stats (see chisq.c)
  
  if(parameters.plotOutput>=1)
    plotSpectra(&parameters);//see plotter.c
  
  if(parameters.saveOutput==1)
    saveSpectra(&parameters);//see save_data.c

  return 0; //great success
}
