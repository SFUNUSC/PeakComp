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
            scaledSimHist[i][j]=scaleFactor*simHist[i][j];
          }
    
      compareSpectra();
    }
  else if(addBackground==1)//linear background addition
    {
    
      computeLinearBackground();//get background coefficients and scaling factor
      
      //scale simulated data
      for (i=startSp;i<=endSp;i++)
        for (j=startCh;j<=endCh;j++)
          scaledSimHist[i][j]=scaleFactor*simHist[i][j];
      //add background to simulated data
      for (i=startSp;i<=endSp;i++)
        for (j=startCh;j<=endCh;j++)
          scaledSimHist[i][j]=scaledSimHist[i][j] + bgA + bgB*j;
      
      
      compareSpectra();
    }
      
  //print output
  printf("\nCOMPARISON DATA\n---------------\n");
  printf("chisq: %f\n",chisq);
  printf("number of bins: %i\n",numBinsUsed);
  printf("chisq / number of bins: %f\n",chisq/(numBinsUsed));
  
  if(plotOutput==1)
    plotSpectra();

  return(0); //great success
}

//function compares spectra and gets chisq and other stats
void compareSpectra()
{

  //initialize values
  binsSkipped=0;
  numBinsUsed=0;
  
  //compute chisq for data in the spectra
  chisq=0;
  for (i=startSp;i<=endSp;i++)
    for (j=startCh;j<=endCh;j++)
      {
        if(expHist[i][j]!=0)//avoid dividing by zero
          chisq+=((expHist[i][j]-scaledSimHist[i][j])*(expHist[i][j]-scaledSimHist[i][j]))/((double)expHist[i][j]);
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

//function computes background coefficients and scaling factor
//using an analytic solution to Cramer's rule minimizing chisq 
//for the expression:
//chisq=sum_i[(meas_i - scaleFactor*sim_i - A - B*i)^2 / meas_i]
void computeLinearBackground()
{
  //get sums
  m_sum=0.;
  s_sum=0.;
  ss_sum=0.;
  ms_sum=0.;
  mi_sum=0.;
  si_sum=0.;
  i_sum=0.;
  ii_sum=0.;
  sum1=0.;
  for (i=startSp;i<=endSp;i++)
    for (j=startCh;j<=endCh;j++)
      if(expHist[i][j]!=0)
        {
          m_sum+=expHist[i][j]/((double)expHist[i][j]);
          s_sum+=simHist[i][j]/((double)expHist[i][j]);
          ss_sum+=simHist[i][j]*simHist[i][j]/((double)expHist[i][j]);
          ms_sum+=expHist[i][j]*simHist[i][j]/((double)expHist[i][j]);
          mi_sum+=expHist[i][j]*j/((double)expHist[i][j]);
          si_sum+=simHist[i][j]*j/((double)expHist[i][j]);
          i_sum+=j/((double)expHist[i][j]);
          ii_sum+=j*j/((double)expHist[i][j]);
          sum1+=1./((double)expHist[i][j]);
        }
      
  //calculate determinants
  detA=ss_sum*(sum1*ii_sum - i_sum*i_sum) - s_sum*(s_sum*ii_sum - i_sum*si_sum) + si_sum*(s_sum*i_sum - sum1*si_sum);
  detAi[0]=ms_sum*(sum1*ii_sum - i_sum*i_sum) - s_sum*(m_sum*ii_sum - i_sum*mi_sum) + si_sum*(m_sum*i_sum - sum1*mi_sum);
  detAi[1]=ss_sum*(m_sum*ii_sum - i_sum*mi_sum) - ms_sum*(s_sum*ii_sum - i_sum*si_sum) + si_sum*(s_sum*mi_sum - m_sum*si_sum);     
  detAi[2]=ss_sum*(sum1*mi_sum - m_sum*i_sum) - s_sum*(s_sum*mi_sum - m_sum*si_sum) + ms_sum*(s_sum*i_sum - sum1*si_sum);
  //get parameters (Cramer's rule)
  scaleFactor=detAi[0]/detA;
  bgA=detAi[1]/detA;
  bgB=detAi[2]/detA;
  printf("Fit linear background of form [A + B*channel], A = %0.3Lf, B = %0.3Lf\n",bgA,bgB);
  printf("Fit scaling factor: %f\n",scaleFactor);
}

//function handles plotting of data, using the gnuplot_i library
void plotSpectra(){

  //generate plot data
  double x[endSp-startSp+1][endCh-startCh+1];
  double y1[endSp-startSp+1][endCh-startCh+1];
  double y2[endSp-startSp+1][endCh-startCh+1];
  for (i=startSp;i<=endSp;i++)
    for (j=startCh;j<=endCh;j++)
      {
        x[i-startSp][j-startCh]=(double)j;
        y1[i-startSp][j-startCh]=(double)expHist[i][j];
        y2[i-startSp][j-startCh]=scaledSimHist[i][j];
      }
      
  handle=gnuplot_init();
  printf("\n\n");
  for(i=0;i<endSp-startSp+1;i++)
    {
      gnuplot_plot_xy(handle, x[i], y1[i], endCh-startCh+1, "Experiment");
      gnuplot_plot_xy(handle, x[i], y2[i], endCh-startCh+1, "Simulation");
      printf("Showing plot for spectrum %i, press any key to continue...", startSp+i);
      getc(stdin);
      gnuplot_resetplot(handle);
    }
  gnuplot_close(handle);
  
}
