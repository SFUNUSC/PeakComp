#include <string.h>

FILE *config;
int startSp,endSp,startCh,endCh;
int addBackground;//0=no,1=constant background
int plotOutput;//0=no,1=yes
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
        strcpy(expDataName,str2);
      if(strcmp(str1,"SIMULATED_DATA")==0)
        strcpy(simDataName,str2);
      if(strcmp(str1,"START_SPECTRUM")==0)
        startSp=atoi(str2);
      if(strcmp(str1,"END_SPECTRUM")==0)
        endSp=atoi(str2);
      if(strcmp(str1,"START_CHANNEL")==0)
        startCh=atoi(str2);
      if(strcmp(str1,"END_CHANNEL")==0)
        endCh=atoi(str2);
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
          else
            plotOutput=0;
        }
    }
  fclose(config);
  
  printf("Experiment data filename is %s.\n",expDataName);
  printf("Simulated data filename is %s.\n",simDataName);
  printf("Comparison starts at spectrum %i.\n",startSp);
  printf("Comparison ends at spectrum %i.\n",endSp);
  printf("Comparison starts at channel %i.\n",startCh);
  printf("Comparison ends at channel %i.\n",endCh);
  if(addBackground==0)
    printf("Will not add background to simulated data.\n");
  if(addBackground==1)
    printf("Will add a linear background to simulated data.\n");
  if(plotOutput==0)
    printf("Will not plot output data.\n");
  if(plotOutput==1)
    printf("Will plot output data.\n");
  
  printf("Finished reading configuration file...\n");
  
}
