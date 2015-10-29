#include <string.h>

FILE *config;
int startSp,endSp,startCh,endCh;
int addBackground;//0=no,1=constant background
char expDataName[256],simDataName[256];//filenames for the simulated and experiment data
char str1[256],str2[256];

void readConfigFile(const char * fileName) 
{

  if((config=fopen(fileName,"r"))==NULL)
    {
      printf("ERROR: Cannot open the config file %s!\n",fileName);
      exit(-1);
    }
  while(fscanf(config,"%s %s",str1,str2)!=EOF)
    {
      //printf("%s %s.\n",str1,str2);
      if(strcmp(str1,"EXPERIMENT_DATA")==0)
        {
          strcpy(expDataName,str2);
          printf("Experiment data filename is %s.\n",expDataName);
        }
      if(strcmp(str1,"SIMULATED_DATA")==0)
        {
          strcpy(simDataName,str2);
          printf("Simulated data filename is %s.\n",simDataName);
        }
      if(strcmp(str1,"START_SPECTRUM")==0)
        {
          startSp=atoi(str2);
          printf("Comparison starts at spectrum %i.\n",startSp);
        }
      if(strcmp(str1,"END_SPECTRUM")==0)
        {
          endSp=atoi(str2);
          printf("Comparison ends at spectrum %i.\n",endSp);
        }
      if(strcmp(str1,"START_CHANNEL")==0)
        {
          startCh=atoi(str2);
          printf("Comparison starts at channel %i.\n",startCh);
        }
      if(strcmp(str1,"END_CHANNEL")==0)
        {
          endCh=atoi(str2);
          printf("Comparison ends at channel %i.\n",endCh);
        }
      if(strcmp(str1,"ADD_BACKGROUND")==0)
        {
          if(strcmp(str2,"yes")==0)
            {
              addBackground=1;
              printf("Will add a linear background to simulated data.\n");
            }
          else
            {
              addBackground=0;
              printf("Will not add background to simulated data.\n");
            }           

        }
    }
  fclose(config);
    
  printf("Finished reading configuration file...\n");
  
}
