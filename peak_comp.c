#include "peak_comp.h"
#include "read_config.c"

int main(int argc, char *argv[])
{

  FILE *expData,*simData;

  if(argc!=2)
    {
      printf("\npeak_comp parameter_file\n");
      printf("Compares the .mca spectra designated in the parameter file specified and generates cool statistics.\n\n");
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
  if(endSpectrum>=NSPECT)
    {
      printf("ERROR: A spectrum number specified in the parameter file is larger than the maximum value of %i.  Reduce it or increase NSPECT in peak_comp.h and recompile.\n",NSPECT);
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
  for (i=0;i<=endSpectrum;i++)
    if(fread(expHist[i],S32K*sizeof(int),1,expData)!=1)
      {
        printf("ERROR: Error reading file %s!\n",expDataName);
        printf("Verify that the format and number of spectra in the file are correct.\n");
        exit(-1);
      }
  for (i=0;i<=endSpectrum;i++)
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
      
      for (i=0;i<numSpectra;i++)
        {
          printf("Calculating integrals of data in spectrum %i...\n",spectrum[i]);
          expInt=0.;
          simInt=0.;
          for (j=startCh[i];j<=endCh[i];j++)
            {
              expInt+=(double)expHist[spectrum[i]][j];
              simInt+=(double)simHist[spectrum[i]][j];
            }
          scaleFactor[spectrum[i]]=expInt/simInt;
          printf("Experiment: %1.0f, Simulated: %1.0f\n",expInt,simInt);
          printf("Scaling simulated data by a factor of: %f\n",scaleFactor[spectrum[i]]);
        }
      
      
      for (i=0;i<numSpectra;i++)
        for (j=startCh[i];j<=endCh[i];j++)
          {
            scaledSimHist[spectrum[i]][j]=scaleFactor[spectrum[i]]*simHist[spectrum[i]][j];
          }
    
      compareSpectra();
    }
  else if(addBackground==1)//linear background addition
    {
    
      computeLinearBackground();//get background coefficients and scaling factor
      
      //scale simulated data
      for (i=0;i<numSpectra;i++)
        for (j=startCh[i];j<=endCh[i];j++)
          scaledSimHist[spectrum[i]][j]=scaleFactor[spectrum[i]]*simHist[spectrum[i]][j];
      //add background to simulated data
      for (i=0;i<numSpectra;i++)
        for (j=startCh[i];j<=endCh[i];j++)
          scaledSimHist[spectrum[i]][j]=scaledSimHist[spectrum[i]][j] + bgA[spectrum[i]] + bgB[spectrum[i]]*j;
      
      
      compareSpectra();
    }
      
  //print output
  printf("\nCOMPARISON DATA\n---------------\n");
  for (i=0;i<numSpectra;i++)
    printf("chisq (spectrum %i): %f\n",spectrum[i],spectChisq[i]);
  printf("chisq (total): %f\n",chisq);
  printf("number of bins: %i\n",numBinsUsed);
  printf("chisq (total) / number of bins: %f\n",chisq/(numBinsUsed));
  
  if(plotOutput==1)
    plotSpectra();

  return 0; //great success
}

//function compares spectra and gets chisq and other stats
void compareSpectra()
{

  //initialize values
  binsSkipped=0;
  numBinsUsed=0;
  
  //compute chisq for data in the spectra
  chisq=0;
  for (i=0;i<numSpectra;i++)
    {
      spectChisq[i]=0;
      for (j=startCh[i];j<=endCh[i];j++)
        {
          if(expHist[spectrum[i]][j]!=0)//avoid dividing by zero
            spectChisq[i]+=((expHist[spectrum[i]][j]-scaledSimHist[spectrum[i]][j])*(expHist[spectrum[i]][j]-scaledSimHist[spectrum[i]][j]))/((double)expHist[spectrum[i]][j]);
          else
            binsSkipped++;
        }
      chisq+=spectChisq[i];
    }
      
  //print warnings
  if(binsSkipped>0)
    {
      printf("Warning: some of the bins in the experiment data have values of zero.  These have been skipped when calculating chisq.\n");
      printf("Bins skipped: %i.\n",binsSkipped);
    }
  for (i=0;i<numSpectra;i++)
    numBinsUsed += (endCh[i]-startCh[i]+1);
  numBinsUsed -= binsSkipped;
}

//function computes background coefficients and scaling factor
//using an analytic solution to Cramer's rule minimizing chisq 
//for the expression:
//chisq=sum_i[(meas_i - scaleFactor*sim_i - A - B*i)^2 / meas_i]
void computeLinearBackground()
{
  printf("\n");  
  //get sums  
  for (i=0;i<numSpectra;i++)
    {
      m_sum=0.;
      s_sum=0.;
      ss_sum=0.;
      ms_sum=0.;
      mi_sum=0.;
      si_sum=0.;
      i_sum=0.;
      ii_sum=0.;
      sum1=0.;
      for (j=startCh[i];j<=endCh[i];j++)
        if(expHist[spectrum[i]][j]!=0)
          {
            m_sum+=expHist[spectrum[i]][j]/((double)expHist[spectrum[i]][j]);
            s_sum+=simHist[spectrum[i]][j]/((double)expHist[spectrum[i]][j]);
            ss_sum+=simHist[spectrum[i]][j]*simHist[spectrum[i]][j]/((double)expHist[spectrum[i]][j]);
            ms_sum+=expHist[spectrum[i]][j]*simHist[spectrum[i]][j]/((double)expHist[spectrum[i]][j]);
            mi_sum+=expHist[spectrum[i]][j]*j/((double)expHist[spectrum[i]][j]);
            si_sum+=simHist[spectrum[i]][j]*j/((double)expHist[spectrum[i]][j]);
            i_sum+=j/((double)expHist[spectrum[i]][j]);
            ii_sum+=j*j/((double)expHist[spectrum[i]][j]);
            sum1+=1./((double)expHist[spectrum[i]][j]);
          }
       
      //calculate determinants
      detA[spectrum[i]]=ss_sum*(sum1*ii_sum - i_sum*i_sum) - s_sum*(s_sum*ii_sum - i_sum*si_sum) + si_sum*(s_sum*i_sum - sum1*si_sum);
      detAi[spectrum[i]][0]=ms_sum*(sum1*ii_sum - i_sum*i_sum) - s_sum*(m_sum*ii_sum - i_sum*mi_sum) + si_sum*(m_sum*i_sum - sum1*mi_sum);
      detAi[spectrum[i]][1]=ss_sum*(m_sum*ii_sum - i_sum*mi_sum) - ms_sum*(s_sum*ii_sum - i_sum*si_sum) + si_sum*(s_sum*mi_sum - m_sum*si_sum);     
      detAi[spectrum[i]][2]=ss_sum*(sum1*mi_sum - m_sum*i_sum) - s_sum*(s_sum*mi_sum - m_sum*si_sum) + ms_sum*(s_sum*i_sum - sum1*si_sum);
      //get parameters (Cramer's rule)
      scaleFactor[spectrum[i]]=detAi[spectrum[i]][0]/detA[spectrum[i]];
      bgA[spectrum[i]]=detAi[spectrum[i]][1]/detA[spectrum[i]];
      bgB[spectrum[i]]=detAi[spectrum[i]][2]/detA[spectrum[i]];
      printf("Spectrum %i: fit linear background of form [A + B*channel], A = %0.3Lf, B = %0.3Lf\n",spectrum[i],bgA[spectrum[i]],bgB[spectrum[i]]);
      printf("Fit scaling factor: %f\n",scaleFactor[spectrum[i]]);   
          
       
     }
     
  
}

//function handles plotting of data, using the gnuplot_i library
void plotSpectra(){

  //generate plot data
  double x[numSpectra][maxNumCh];
  double y1[numSpectra][maxNumCh];
  double y2[numSpectra][maxNumCh];
  for (i=0;i<numSpectra;i++)
    for (j=startCh[i];j<=endCh[i];j++)
      {
        x[i][j-startCh[i]]=(double)j;
        y1[i][j-startCh[i]]=(double)expHist[spectrum[i]][j];
        y2[i][j-startCh[i]]=scaledSimHist[spectrum[i]][j];
      }
      
  handle=gnuplot_init();
  printf("\n");
  for(i=0;i<numSpectra;i++)
    {
      gnuplot_setstyle(handle,"steps");
      gnuplot_cmd(handle,"set xlabel 'Channel'");
      gnuplot_cmd(handle,"set ylabel 'Counts'");
      gnuplot_plot_xy(handle, x[i], y1[i], endCh[i]-startCh[i]+1, "Experiment");
      gnuplot_plot_xy(handle, x[i], y2[i], endCh[i]-startCh[i]+1, "Simulation");
      printf("Showing plot for spectrum %i, press [ENTER] to continue...", spectrum[i]);
      getc(stdin);
      gnuplot_resetplot(handle);
    }
  gnuplot_close(handle);
  
}
