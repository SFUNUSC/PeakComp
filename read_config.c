#include "peak_comp.h"
#include <string.h>

FILE *config;
int spectrum[NSPECT],startCh[NSPECT],endCh[NSPECT],numSpectra,endSpectrum,maxNumCh,numSimData;
int addBackground;//0=no,1=constant background
int plotOutput;//0=no,1=yes,2=detailed
char expDataName[256],simDataName[NSIMDATA][256];//filenames for the simulated and experiment data
int simDataFixedAmp[NSIMDATA];//bool specifying whether amplitude of each set of simulated data is fixed
double simDataFixedAmpValue[NSIMDATA];//value at which amplitude is fixed for each set of simulated data
char str[256],str1[256],str2[256];


void readConfigFile(const char * fileName) 
{
  int index=0;
  numSpectra=0;
  endSpectrum=0;
  maxNumCh=0;
  numSimData=0;
  if((config=fopen(fileName,"r"))==NULL)
    {
      printf("ERROR: Cannot open the config file %s!\n",fileName);
      exit(-1);
    }
  while(!(feof(config)))//go until the end of file is reached
    {
      if(fgets(str,256,config)!=NULL)
        {
        
          if(numSimData<NSIMDATA)
            if(sscanf(str,"%i %i %i",&spectrum[index],&startCh[index],&endCh[index])!=3) //no spectrum and channel data
              if(sscanf(str,"%s %s %lf",simDataName[numSimData],str1,&simDataFixedAmpValue[numSimData])==3) //simulated dataset info
                {
                  if(strcmp(str1,"yes")==0)
                    simDataFixedAmp[numSimData]=1;
                  else
                    simDataFixedAmp[numSimData]=0;
                  numSimData++;
                }
              
          if(index<NSPECT)
            if(sscanf(str,"%i %i %i",&spectrum[index],&startCh[index],&endCh[index])==3) //spectrum and channel data
              {
                if(spectrum[index]>endSpectrum)
                  endSpectrum=spectrum[index];
                if((endCh[index]-startCh[index]+1)>maxNumCh)
                  maxNumCh=endCh[index]-startCh[index]+1;
                index++;
                numSpectra++;
              }
              
          if(sscanf(str,"%s %s",str1,str2)==2) //single parameter data
            {
              if(strcmp(str1,"EXPERIMENT_DATA")==0)
                strcpy(expDataName,str2);
              if(strcmp(str1,"ADD_BACKGROUND")==0)
                {
                  if(strcmp(str2,"yes")==0)
                    addBackground=1;
                  else
                    addBackground=0;
                }
              if(strcmp(str1,"PLOT_OUTPUT")==0)
                {
                  if(strcmp(str2,"yes")==0)
                    plotOutput=1;
                  else if(strcmp(str2,"detailed")==0)
                    plotOutput=2;
                  else
                    plotOutput=0;
                }
            }
          
          if(sscanf(str,"%s %s",str1,str2)==1) //listing of simulated data
            {
              if(strcmp(str1,"<---END_OF_PARAMETERS--->")==0)
                break;
              else if(strcmp(str1,"SIMULATED_DATA")!=0)
                if(numSimData<NSIMDATA)
                  {
                    strcpy(simDataName[numSimData],str1);
                    simDataFixedAmp[numSimData]=0;
                    simDataFixedAmpValue[numSimData]=1;
                    numSimData++;
                  }
            }
        }
    }
  fclose(config);
  
  printf("Taking experiment data from file: %s\n",expDataName);
  for(index=0;index<numSimData;index++)
    {
      printf("Taking simulated data from file (%i of %i): %s\n",index+1,numSimData,simDataName[index]);
      if(simDataFixedAmp[index]==1)
        printf("Fixing scaling factor for this data to %lf\n",simDataFixedAmpValue[index]);
    }
  for(index=0;index<numSpectra;index++)
    printf("Will compare spectrum %i from channels %i to %i.\n",spectrum[index],startCh[index],endCh[index]);
  if(addBackground==0)
    printf("Will not add background to simulated data.\n");
  if(addBackground==1)
    printf("Will add a linear background to simulated data.\n");
  if(plotOutput==0)
    printf("Will not plot output data.\n");
  if(plotOutput==1)
    printf("Will plot output data.\n");
  if(plotOutput==2)
    printf("Will plot detailed output data.\n");
  
  printf("Finished reading parameter file...\n");
  
}
