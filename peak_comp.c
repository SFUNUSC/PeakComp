#include "peak_comp.h"
#include "read_config.c"


int main(int argc, char *argv[])
{

  FILE *expData,*simData[NSIMDATA];
  
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
  
  //initialize values
  addBackground=0;
  memset(expHist,0,sizeof(expHist));
  memset(fittedExpHist,0,sizeof(fittedExpHist));
  memset(simHist,0,sizeof(simHist));
  memset(fittedSimHist,0,sizeof(fittedSimHist));
  memset(scaledSimHist,0,sizeof(scaledSimHist));

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
  for (i=0;i<numSimData;i++)  
    if((simData[i]=fopen(simDataName[i],"r"))==NULL)
      {
        printf("ERROR: Cannot open the simulated data file %s!\n",simDataName[i]);
        exit(-1);
      }
  for (i=0;i<=endSpectrum;i++)
    {
      if(fread(expHist[i],S32K*sizeof(int),1,expData)!=1)
        {
          printf("ERROR: Error reading file %s!\n",expDataName);
          printf("Verify that the format and number of spectra in the file are correct.\n");
          exit(-1);
        }
      for (j=0;j<S32K;j++)
        fittedExpHist[i][j]=expHist[i][j];
    }
  int fi=0;//index for simulated data to be fitted
  for (i=0;i<numSimData;i++)
    {
      //read all simulated data in
      for (j=0;j<=endSpectrum;j++)
        if(fread(simHist[i][j],S32K*sizeof(int),1,simData[i])!=1)
          {
            printf("ERROR: Error reading file %s!\n",simDataName[i]);
            printf("Verify that the format and number of spectra in the file are correct.\n");
            exit(-1);
          }
      //determine whether simulated data is fitted and read into histograms for fitting as needed
      if(simDataFixedAmp[i]==0)
        {
          for (j=0;j<=endSpectrum;j++)
            for (k=0;k<S32K;k++)
              fittedSimHist[fi][j][k]=simHist[i][j][k];
          fi++;
        }
      else if(simDataFixedAmp[i]==2)//data scaling is fixed relative to the previous fitted dataset
        {
          if(fi>0)
            for (j=0;j<=endSpectrum;j++)
              for (k=0;k<S32K;k++)
                fittedSimHist[fi-1][j][k]+=simDataFixedAmpValue[i]*simHist[i][j][k];//add this data to the data that it is scaled relative to
        }
      else if(simDataFixedAmpValue[i]!=0.)//data scaling is fixed to a specified value (will not be fitted)
        for (j=0;j<=endSpectrum;j++)
          for (k=0;k<S32K;k++)
            fittedExpHist[j][k]-=simDataFixedAmpValue[i]*simHist[i][j][k];
    }
  fclose(expData);
  for (i=0;i<numSimData;i++) 
    fclose(simData[i]);
  printf("Spectra read in...\n");

  computeBackgroundandScaling(numFittedSimData,addBackground);//get background coefficients and scaling factors
      
  //scale simulated data
  for (i=0;i<numSimData;i++)
    for (j=0;j<numSpectra;j++)
      for (k=startCh[j];k<=endCh[j];k++)
        scaledSimHist[i][spectrum[j]][k]=scaleFactor[i][spectrum[j]]*simHist[i][spectrum[j]][k];
  //add background to simulated data
  if(addBackground==1)
    for (i=0;i<numSimData;i++)
      for (j=0;j<numSpectra;j++)
        for (k=startCh[j];k<=endCh[j];k++)
          scaledSimHist[i][spectrum[j]][k]=scaledSimHist[i][spectrum[j]][k] + bgA[spectrum[j]] + bgB[spectrum[j]]*k;
      
  compareSpectra();
  
  if(plotOutput>=1)
    plotSpectra();

  return 0; //great success
}


//function compares spectra and gets chisq and other stats
void compareSpectra()
{
  //initialize values
  double chisq=0;
  double spectChisq[NSPECT];
  int binsSkipped=0;
  int numBinsUsed=0;
  double sumSimValue=0;
  
  //compute chisq for data in the spectra
  for (i=0;i<numSpectra;i++)
    {
      spectChisq[i]=0;
      for (j=startCh[i];j<=endCh[i];j++)
        {  
          if(expHist[spectrum[i]][j]!=0)//avoid dividing by zero
            {
              //get the sum of all experimental data in the given bin 
              sumSimValue = bgA[spectrum[i]] + bgB[spectrum[i]]*j;
              for (k=0;k<numSimData;k++)
                sumSimValue+=scaledSimHist[k][spectrum[i]][j] - bgA[spectrum[i]] - bgB[spectrum[i]]*j;
              //increment the chisq value
              spectChisq[i]+=((expHist[spectrum[i]][j]-sumSimValue)*(expHist[spectrum[i]][j]-sumSimValue))/((double)expHist[spectrum[i]][j]);
            }
          else
            binsSkipped++;
        }
      chisq+=spectChisq[i];
    }
      
  //print warnings
  if(binsSkipped>0)
    {
      printf("\nWarning: some of the bins in the experiment data have values of zero.  These have been skipped when calculating chisq.\n");
      printf("Bins skipped: %i.\n",binsSkipped);
    }
  for (i=0;i<numSpectra;i++)
    numBinsUsed += (endCh[i]-startCh[i]+1);
  numBinsUsed -= binsSkipped;
  
  //print output
  printf("\nCOMPARISON DATA\n---------------\n");
  for (i=0;i<numSpectra;i++)
    printf("chisq (spectrum %i): %f\n",spectrum[i],spectChisq[i]);
  printf("chisq (total): %f\n",chisq);
  printf("number of bins: %i\n",numBinsUsed);
  printf("chisq (total) / number of bins: %f\n",chisq/(numBinsUsed));
  
}


//function computes background coefficients and scaling factors
//by analytically minimizing chisq for the expression:
//chisq=sum_i[(meas_i - A - B*i - scaleFactor_1*sim_1i - scaleFactor_2*sim_2i - ...)^2 / meas_i]
//addBG=0: no background addition
void computeBackgroundandScaling(int numData, int addBG)
{

  long double m_sum,s_sum[NSIMDATA],ss_sum[NSIMDATA][NSIMDATA],ms_sum[NSIMDATA],mi_sum,si_sum[NSIMDATA],i_sum,ii_sum,sum1; //sums needed to construct system of equations
  lin_eq_type linEq;
  printf("\n");  
  
  if(numData>0)
    {
      for (i=0;i<numSpectra;i++)
        {
          //initialize all sums to 0
          m_sum=0.;
          mi_sum=0.;
          i_sum=0.;
          ii_sum=0.;
          sum1=0.;
          for (j=0;j<numData;j++)
            {
              s_sum[j]=0.;
              ms_sum[j]=0.;
              si_sum[j]=0.;
              for (k=0;k<numData;k++)
                ss_sum[j][k]=0.;
            }
          
          //construct sums  
          for (j=startCh[i];j<=endCh[i];j++)
            if(expHist[spectrum[i]][j]!=0)
              {
                m_sum+=fittedExpHist[spectrum[i]][j]/((double)expHist[spectrum[i]][j]);
                mi_sum+=fittedExpHist[spectrum[i]][j]*j/((double)expHist[spectrum[i]][j]);
                i_sum+=j/((double)expHist[spectrum[i]][j]);
                ii_sum+=j*j/((double)expHist[spectrum[i]][j]);
                sum1+=1./((double)expHist[spectrum[i]][j]);
                for (k=0;k<numData;k++)
                  {
                    s_sum[k]+=fittedSimHist[k][spectrum[i]][j]/((double)expHist[spectrum[i]][j]);
                    ms_sum[k]+=fittedExpHist[spectrum[i]][j]*fittedSimHist[k][spectrum[i]][j]/((double)expHist[spectrum[i]][j]);
                    si_sum[k]+=fittedSimHist[k][spectrum[i]][j]*j/((double)expHist[spectrum[i]][j]);
                    for (l=0;l<numData;l++)
                      ss_sum[k][l]+=fittedSimHist[k][spectrum[i]][j]*fittedSimHist[l][spectrum[i]][j]/((double)expHist[spectrum[i]][j]);
                  }
              }
          
          //construct system of equations (matrix/vector entries) 
          if(addBG==0)
            {
              linEq.dim=numData;
              for (j=0;j<numData;j++)
                for (k=0;k<numData;k++)
                  linEq.matrix[j][k]=ss_sum[j][k];
              for (j=0;j<numData;j++)
                linEq.vector[j]=ms_sum[j];
            }
          else
            {
              linEq.dim=numData+2;
              
              //top-left 4 entries
              linEq.matrix[0][0]=sum1;
              linEq.matrix[0][1]=i_sum;
              linEq.matrix[1][0]=i_sum;
              linEq.matrix[1][1]=ii_sum;
              
              //regular simulated data entires (bottom-right)
              for (j=0;j<numData;j++)
                for (k=0;k<numData;k++)
                  linEq.matrix[j+2][k+2]=ss_sum[j][k];
              
              //remaining entires
              for (j=0;j<numData;j++)
                {     
                  linEq.matrix[0][2+j]=s_sum[j];
                  linEq.matrix[1][2+j]=si_sum[j];
                  linEq.matrix[2+j][0]=s_sum[j];
                  linEq.matrix[2+j][1]=si_sum[j];
                }
              
              linEq.vector[0]=m_sum;
              linEq.vector[1]=mi_sum;
              for (j=0;j<numData;j++)
                linEq.vector[j+2]=ms_sum[j];
            }
          
          //solve system of equations and assign values
          if(!(solve_lin_eq(&linEq)==1))
            {
              if(m_sum==0)
                printf("ERROR: Experiment data (spectrum %i) has no entires in the specified fitting region!\n",spectrum[i]);
              else
                printf("ERROR: Could not determine background and scaling parameters!\n");
              exit(-1);
            }
          
          if(addBG==0)
            {
              bgA[spectrum[i]]=0.;
              bgB[spectrum[i]]=0.;
              for (j=0;j<numData;j++)  
                fittedScaleFactor[j][spectrum[i]]=linEq.solution[j];
            }
          else
            {
              bgA[spectrum[i]]=linEq.solution[0];
              bgB[spectrum[i]]=linEq.solution[1];
              for (j=0;j<numData;j++)
                fittedScaleFactor[j][spectrum[i]]=linEq.solution[j+2];
            }

        }
    }
  else
    printf("NOTE: All parameters are fixed, no chisq minimzation was performed.\n");
    
  //generate scaling factors for all spectra, including those that weren't fitted
  int fd=0;//counter for number of datasets which have fit (not fixed amplitude)
  int ld=-1;//index of the last dataset for which a scaling factor was determined
  for (i=0;i<numSimData;i++)
    {
      if(simDataFixedAmp[i]==0)//data was fit
        {
          for (j=0;j<numSpectra;j++)
            scaleFactor[i][spectrum[j]]=fittedScaleFactor[fd][spectrum[j]];
          fd++;
          ld=i;
        }
      else if(simDataFixedAmp[i]==2)//data is scaled relative to the previous fit data
        {
          if(ld>=0)//has a previous dataset been fit?
            for (j=0;j<numSpectra;j++)
              scaleFactor[i][spectrum[j]]=simDataFixedAmpValue[i]*scaleFactor[ld][spectrum[j]];
          ld=i;
        }
      else//data wasn't fit
        {
          for (j=0;j<numSpectra;j++)
            scaleFactor[i][spectrum[j]]=simDataFixedAmpValue[i];
          ld=i;
        }
    }
  
  //print parameters  
  for (i=0;i<numSpectra;i++)
    if(addBG==0)
      {
        printf("Spectrum %i - ",spectrum[i]);
        for (j=0;j<numSimData;j++)
          if(simDataFixedAmp[j]==0)
            printf("Scaling factor for data from file %s: %f\n",simDataName[j],scaleFactor[j][spectrum[i]]);
          else
            printf("Scaling factor for data from file %s: %f [FIXED]\n",simDataName[j],scaleFactor[j][spectrum[i]]);
      }
    else
      {
        printf("Spectrum %i: fit linear background of form [A + B*channel], A = %0.3Lf, B = %0.3Lf\n",spectrum[i],bgA[spectrum[i]],bgB[spectrum[i]]);
        for (j=0;j<numSimData;j++)
          if(simDataFixedAmp[j]==0)
            printf("Scaling factor for data from file %s: %f\n",simDataName[j],scaleFactor[j][spectrum[i]]);
          else
            printf("Scaling factor for data from file %s: %f [FIXED]\n",simDataName[j],scaleFactor[j][spectrum[i]]);
      }
  
}


//function handles plotting of data, using the gnuplot_i library
void plotSpectra()
{
  char str[256];
  
  //generate plot data
  double x[numSpectra][maxNumCh];
  double yexp[numSpectra][maxNumCh];
  double ysim[NSIMDATA][numSpectra][maxNumCh];
  double ybackground[numSpectra][maxNumCh];
  double ysimsum[numSpectra][maxNumCh];
  
  for (i=0;i<numSpectra;i++)
    for (j=startCh[i];j<=endCh[i];j++)
      {
        x[i][j-startCh[i]]=(double)j;
        yexp[i][j-startCh[i]]=(double)expHist[spectrum[i]][j];
        if(addBackground==1)
          ybackground[i][j-startCh[i]]=bgA[spectrum[i]] + bgB[spectrum[i]]*j;
        ysimsum[i][j-startCh[i]]=ybackground[i][j-startCh[i]];
        for (k=0;k<numSimData;k++)
          {
            ysim[k][i][j-startCh[i]]=scaledSimHist[k][spectrum[i]][j] - ybackground[i][j-startCh[i]];
            ysimsum[i][j-startCh[i]]+=scaledSimHist[k][spectrum[i]][j] - ybackground[i][j-startCh[i]];
          }
      }

  plotOpen=1; 
  handle=gnuplot_init();
  printf("\n");
  for(i=0;i<numSpectra;i++)
    {
      gnuplot_setstyle(handle,"steps");
      gnuplot_cmd(handle,"set xlabel 'Channel'");
      gnuplot_cmd(handle,"set ylabel 'Counts'");
      gnuplot_plot_xy(handle, x[i], yexp[i], endCh[i]-startCh[i]+1, "Experiment");
      if(plotOutput>1)//detailed plot
        {
          if(addBackground==1)//plot background
            {
              gnuplot_setstyle(handle,"lines"); //plot background as a line
              gnuplot_plot_xy(handle, x[i], ybackground[i], endCh[i]-startCh[i]+1, "Background");
              gnuplot_setstyle(handle,"steps"); //set the plot style back
            }
          for (j=0;j<numSimData;j++)//plot individual sim data
            {
              sprintf(str,"Simulation (%s)",simDataName[j]);
              gnuplot_plot_xy(handle, x[i], ysim[j][i], endCh[i]-startCh[i]+1, str);
            }
          if((numSimData>1)||(addBackground==1))//plot sum
            {
              gnuplot_setcolor(handle, "black");
              gnuplot_plot_xy(handle, x[i], ysimsum[i], endCh[i]-startCh[i]+1, "Simulation and Background(sum)");
              gnuplot_unsetcolor(handle);
            }
        }
      else //simple plot
        gnuplot_plot_xy(handle, x[i], ysimsum[i], endCh[i]-startCh[i]+1, "Simulation and Background(sum)");
      if(i!=numSpectra-1)//check whether we're showing the last plot
        printf("Showing plot for spectrum %i, press [ENTER] to continue...", spectrum[i]);
      else
        printf("Showing plot for spectrum %i, press [ENTER] to exit.", spectrum[i]);
      getc(stdin);
      gnuplot_resetplot(handle);
    }
  gnuplot_close(handle);
  plotOpen=0;
  
}

//function run after CTRL-C, used to clean up temporary files generated
//by the plotting library
void sigint_cleanup()
{
  if(plotOpen==1)
    gnuplot_close(handle); //cleans up temporary files  
  exit(1); 
}
