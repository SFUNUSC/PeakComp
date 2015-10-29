#include "peak_comp.h"
#include "read_config.c"

int main(int argc, char *argv[])
{

  FILE *expData,*simData;

  if(argc!=2)
    {
      printf("\npeak_comp config_file\n");
      printf("Compares the .mca spectra designated in the configuration file specified and generates cool statistics.\n\n");
      exit(-1);
    }
  printf("\n");
  
  //initialize values
  addBackground=0;
  for (i=0;i<NSPECT;i++)
    for (j=0;j<S32K;j++)
      {
        expHist[i][j]=0;
        simHist[i][j]=0;
        scaledSimHist[i][j]=0.;
      }

  readConfigFile(argv[1]); //grab data from the config file

  //check that the number of spectra being compared is fine
  if(endSp>=NSPECT)
    {
      printf("ERROR: END_SPECTRUM in the config file is too large!  Reduce it to %i or lower, or increase NSPECT in peak_comp.h and recompile.\n",NSPECT);
      exit(-1);
    }

  //read in the .mca files
  if((expData=fopen(expDataName,"r"))==NULL)
    {
      printf("ERROR: Cannot open the experiment data file %s!\n",expDataName);
      exit(-1);
    }
  if((simData=fopen(simDataName,"r"))==NULL)
    {
      printf("ERROR: Cannot open the simulated data file %s!\n",simDataName);
      exit(-1);
    }
  for (i=0;i<=endSp;i++)
    if(fread(expHist[i],S32K*sizeof(int),1,expData)!=1)
      {
        printf("ERROR: Error reading file %s!\n",expDataName);
        printf("Verify that the format and number of spectra in the file are correct.\n");
        exit(-1);
      }
  for (i=0;i<=endSp;i++)
    if(fread(simHist[i],S32K*sizeof(int),1,simData)!=1)
      {
        printf("ERROR: Error reading file %s!\n",simDataName);
        printf("Verify that the format and number of spectra in the file are correct.\n");
        exit(-1);
      }
  fclose(expData);
  fclose(simData);
  printf("Spectra read in...\n");
  
  if(addBackground==0)//no background addition
    {
      compareSpectra();
    }
  else if(addBackground==1)//constant background addition
    {
      //get the original histogram values
      for (i=startSp;i<=endSp;i++)
        for (j=startCh;j<=endCh;j++)
          origSimHist[i][j]=simHist[i][j];
    
      minchisq=BIG_NUMBER;
      for (k=minBackground;k<=maxBackground;k++)
        {
          printf("Constant background of %i count(s) added to simulated data...\n",k);
          //add a count to each channel
          for (i=startSp;i<=endSp;i++)
            for (j=startCh;j<=endCh;j++)
              simHist[i][j]=origSimHist[i][j]+k;
          //do the chisq analysis
          compareSpectra();
          printf("chisq using this background: %f\n\n",chisq);
          //determine whether this result the best so far
          if(chisq<minchisq)
            {
              minchisq=chisq;
              bestBackground=k;
              bestBinsSkipped=binsSkipped;
            }
        }
      
      chisq=minchisq;
      binsSkipped=bestBinsSkipped;
      printf("\nchisq is minimized for a constant background of %i count(s).\n",bestBackground);
      printf("Data using this background follows...\n");
      
    }
      
  //print output
  printf("\nCOMPARISON DATA\n---------------\n");
  printf("chisq: %f\n",chisq);
  printf("number of bins: %i\n",numBinsUsed);
  printf("chisq / number of bins: %f\n",chisq/(numBinsUsed));

  return(0); //great success
}

//function compares spectra and gets chisq and other stats
void compareSpectra()
{

  //initialize values
  binsSkipped=0;
  numBinsUsed=0;
  
  //calculate integrals of data in each spectrum and scale the experiment and sim to each other
  printf("Calculating integrals of data over all spectra and channels...\n");
  expInt=0.;
  simInt=0.;
  for (i=startSp;i<=endSp;i++)
    for (j=startCh;j<=endCh;j++)
      {
        expInt+=(double)expHist[i][j];
        simInt+=(double)simHist[i][j];
      }
  scaleFactor=expInt/simInt;
  printf("Experiment: %1.0f, Simulated: %1.0f\n",expInt,simInt);
  printf("Scaling simulated data by a factor of: %f\n",scaleFactor);
  for (i=startSp;i<=endSp;i++)
    for (j=startCh;j<=endCh;j++)
      {
        scaledSimHist[i][j]=(scaleFactor*simHist[i][j]);
      }
  
  //compute chisq for data in the spectra
  chisq=0;
  for (i=startSp;i<=endSp;i++)
    for (j=startCh;j<=endCh;j++)
      {
        if(expHist[i][j]!=0)//avoid dividing by zero
          chisq+=((scaledSimHist[i][j]-expHist[i][j])*(scaledSimHist[i][j]-expHist[i][j]))/((double)expHist[i][j]);
        else
          binsSkipped++;
      }
      
  //print warnings
  if(binsSkipped>0)
    {
      printf("Warning: some of the bins in the experiment data have values of zero.  These have been skipped when calculating chisq.\n");
      printf("Bins skipped: %i.\n",binsSkipped);
    }

  numBinsUsed = (endCh-startCh+1)*(endSp-startSp+1) - binsSkipped;
}
