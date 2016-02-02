#include "peak_comp.h"
#include <string.h>

void readConfigFile(const char * fileName, pc_par * par) 
{
  FILE *config;
  char str[256],str1[256],str2[256];
  int index=0;
  par->numSpectra=0;
  par->endSpectrum=0;
  par->maxNumCh=0;
  par->numSimData=0;
  par->numFittedSimData=0;
  if((config=fopen(fileName,"r"))==NULL)
    {
      printf("ERROR: Cannot open the config file %s!\n",fileName);
      exit(-1);
    }
  while(!(feof(config)))//go until the end of file is reached
    {
      if(fgets(str,256,config)!=NULL)
        {
        
          if(par->numSimData<NSIMDATA)
            if(sscanf(str,"%i %i %i",&par->spectrum[index],&par->startCh[index],&par->endCh[index])!=3) //no spectrum and channel data
              if(sscanf(str,"%s %s %lf",par->simDataName[par->numSimData],str1,&par->simDataFixedAmpValue[par->numSimData])==3) //simulated dataset info
                {
                  if(strcmp(str1,"yes")==0)
                    par->simDataFixedAmp[par->numSimData]=1;
                  else if(strcmp(str1,"rel")==0)
                    par->simDataFixedAmp[par->numSimData]=2;
                  else
                    {
                      par->simDataFixedAmp[par->numSimData]=0;
                      strcpy(par->fittedSimDataName[par->numFittedSimData],par->simDataName[par->numSimData]);
                      par->numFittedSimData++;
                    }
                  par->numSimData++;
                }
              
          if(index<NSPECT)
            if(sscanf(str,"%i %i %i",&par->spectrum[index],&par->startCh[index],&par->endCh[index])==3) //spectrum and channel data
              {
                if(par->spectrum[index]>par->endSpectrum)
                  par->endSpectrum=par->spectrum[index];
                if((par->endCh[index]-par->startCh[index]+1)>par->maxNumCh)
                  par->maxNumCh=par->endCh[index]-par->startCh[index]+1;
                index++;
                par->numSpectra++;
              }
              
          if(sscanf(str,"%s %s",str1,str2)==2) //single parameter data
            {
              if(strcmp(str1,"EXPERIMENT_DATA")==0)
                strcpy(par->expDataName,str2);
              if(strcmp(str1,"ADD_BACKGROUND")==0)
                {
                  if(strcmp(str2,"yes")==0)
                    par->addBackground=1;
                  else
                    par->addBackground=0;
                }
              if(strcmp(str1,"PLOT_OUTPUT")==0)
                {
                  if(strcmp(str2,"yes")==0)
                    par->plotOutput=1;
                  else if(strcmp(str2,"detailed")==0)
                    par->plotOutput=2;
                  else
                    par->plotOutput=0;
                }
              if(strcmp(str1,"PLOT_STYLE")==0)
                {
                  if(strcmp(str2,"lin")==0)
                    par->plotStyle=0;
                  else if(strcmp(str2,"log")==0)
                    par->plotStyle=1;
                  else
                    par->plotStyle=0;
                }
              if(strcmp(str1,"SAVE_OUTPUT")==0)
                {
                  if(strcmp(str2,"yes")==0)
                    par->saveOutput=1;
                  else
                    par->saveOutput=0;
                }
            }
          
          if(sscanf(str,"%s %s",str1,str2)==1) //listing of simulated data
            {
              if(strcmp(str1,"<---END_OF_PARAMETERS--->")==0)
                break;
              else if(strcmp(str1,"SIMULATED_DATA")!=0)
                if(par->numSimData<NSIMDATA)
                  {
                    strcpy(par->simDataName[par->numSimData],str1);
                    par->simDataFixedAmp[par->numSimData]=0;
                    par->simDataFixedAmpValue[par->numSimData]=1.;
                    par->numFittedSimData++;
                    par->numSimData++;
                  }
            }
        }
    }
  fclose(config);
  
  //print parameters read from the file
  printf("Taking experiment data from file: %s\n",par->expDataName);
  for(index=0;index<par->numSimData;index++)
    {
      printf("Taking simulated data from file (%i of %i): %s\n",index+1,par->numSimData,par->simDataName[index]);
      if(par->simDataFixedAmp[index]==1)
        printf("Fixing scaling factor for this data to %lf\n",par->simDataFixedAmpValue[index]);
      if(par->simDataFixedAmp[index]==2)
        printf("Fixing scaling factor for this data to a factor of %lf relative to the last fitted data.\n",par->simDataFixedAmpValue[index]);
    }
  for(index=0;index<par->numSpectra;index++)
    printf("Will compare spectrum %i from channels %i to %i.\n",par->spectrum[index],par->startCh[index],par->endCh[index]);
  if(par->addBackground==0)
    printf("Will not add background to simulated data.\n");
  if(par->addBackground==1)
    printf("Will add a linear background to simulated data.\n");
  if(par->plotOutput==0)
    printf("Will not plot output data.\n");
  if(par->plotOutput==1)
    printf("Will plot output data.\n");
  if(par->plotOutput==2)
    printf("Will plot detailed output data.\n");
  if(par->plotOutput>0)
    if(par->plotStyle==1)
      printf("Will plot using logarithmic y-axis.\n");
  if(par->saveOutput==0)
    printf("Will not save fitted simulation data.\n");
  if(par->saveOutput==1)
    {
      if(par->addBackground==0)
        {
          if(par->numSimData>1)
            printf("Will save fitted simulation data to files fit_sim0.mca through fit_sim%i.mca.\n",par->numSimData-1);
          else
            printf("Will save fitted simulation data to file fit_sim0.mca.\n");
        }
      else if(par->addBackground==1)
        {
          if(par->numSimData>1)
            printf("Will save fitted simulation data to files fit_background.mca and fit_sim0.mca through fit_sim%i.mca.\n",par->numSimData-1);
          else
            printf("Will save fitted simulation data to files fit_background.mca and fit_sim0.mca.\n");
        }
    }

  
  printf("Finished reading parameter file...\n");
  
}
