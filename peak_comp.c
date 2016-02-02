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

  int i,j,k;
  pc_par parameters;
  
  //initialize values
  parameters.addBackground=0;
  memset(expHist,0,sizeof(expHist));
  memset(fittedExpHist,0,sizeof(fittedExpHist));
  memset(simHist,0,sizeof(simHist));
  memset(fittedSimHist,0,sizeof(fittedSimHist));
  memset(scaledSimHist,0,sizeof(scaledSimHist));

  readConfigFile(argv[1],&parameters); //grab data from the config file

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
      else if(parameters.simDataFixedAmpValue[i]!=0.)//data scaling is fixed to a specified value (will not be fitted)
        for (j=0;j<=parameters.endSpectrum;j++)
          for (k=0;k<S32K;k++)
            fittedExpHist[j][k]-=parameters.simDataFixedAmpValue[i]*simHist[i][j][k];
    }
  fclose(expData);
  for (i=0;i<parameters.numSimData;i++) 
    fclose(simData[i]);
  printf("Spectra read in...\n");

  computeBackgroundandScaling(&parameters);//get background coefficients and scaling factors
      
  //scale simulated data
  for (i=0;i<parameters.numSimData;i++)
    for (j=0;j<parameters.numSpectra;j++)
      for (k=0;k<S32K;k++)
        scaledSimHist[i][parameters.spectrum[j]][k]=scaleFactor[i][parameters.spectrum[j]]*simHist[i][parameters.spectrum[j]][k];
      
  compareSpectra(&parameters);
  
  if(parameters.plotOutput>=1)
    plotSpectra(&parameters);
  
  if(parameters.saveOutput==1)
    saveSpectra(&parameters);

  return 0; //great success
}


//function compares spectra and gets chisq and other stats
void compareSpectra(pc_par * par)
{
  //initialize values
  double chisq=0;
  double redChisq=0;
  double spectChisq[NSPECT];
  double spectRedChisq[NSPECT];
  int binsSkipped[NSPECT];
  int numBinsUsed[NSPECT];
  int sumBinsUsed=0;
  int sumBinsSkipped=0;
  double sumSimValue=0;
  int i,j,k;
  memset(spectRedChisq,0,sizeof(spectRedChisq));
  memset(spectChisq,0,sizeof(spectChisq));
  memset(binsSkipped,0,sizeof(binsSkipped));
  memset(numBinsUsed,0,sizeof(numBinsUsed));
  
  int numFittedParameters = par->numFittedSimData + par->addBackground*2;
  
  //compute chisq for data in the spectra
  for (i=0;i<par->numSpectra;i++)
      for (j=par->startCh[i];j<=par->endCh[i];j++)
          if(expHist[par->spectrum[i]][j]!=0)//avoid dividing by zero
            {
              //get the sum of all experimental data in the given bin 
              sumSimValue = bgA[par->spectrum[i]] + bgB[par->spectrum[i]]*j;
              for (k=0;k<par->numSimData;k++)
                sumSimValue+=scaledSimHist[k][par->spectrum[i]][j];
              //increment the chisq value
              spectChisq[i]+=((expHist[par->spectrum[i]][j]-sumSimValue)*(expHist[par->spectrum[i]][j]-sumSimValue))/((double)expHist[par->spectrum[i]][j]);
            }
          else
            binsSkipped[i]++;
      
  //print warnings
  for (i=0;i<par->numSpectra;i++)
    sumBinsSkipped+=binsSkipped[i];
  if(sumBinsSkipped>0)
    {
      printf("\nWarning: some of the bins in the experiment data have values of zero.  These have been skipped when calculating chisq.\n");
      printf("Bins skipped: %i.\n",sumBinsSkipped);
    }
    
  //compute total chisq and reduced total chisq
  for (i=0;i<par->numSpectra;i++)
    {
      numBinsUsed[i]=(par->endCh[i]-par->startCh[i]+1)-binsSkipped[i];
      sumBinsUsed+=numBinsUsed[i];
      chisq+=spectChisq[i]; 
      spectRedChisq[i]=spectChisq[i]/(numBinsUsed[i]-numFittedParameters-1);
    }
  redChisq=chisq/(sumBinsUsed-numFittedParameters-1);
  
  //print output
  printf("\nCOMPARISON DATA\n---------------\n");
  if(par->numSpectra>1)
    for (i=0;i<par->numSpectra;i++)
      printf("spectrum %i, channel %i to %i - chisq: %f, reduced chisq: %f\n",par->spectrum[i],par->startCh[i],par->endCh[i],spectChisq[i],spectRedChisq[i]);
  printf("chisq (total): %f\n",chisq);
  printf("number of bins (total): %i\n",sumBinsUsed);
  printf("number of fitted parameters: %i\n",numFittedParameters);
  printf("reduced chisq (total): %f\n",redChisq);
  
}


//function computes background coefficients and scaling factors
//by analytically minimizing chisq for the expression:
//chisq=sum_i[(meas_i - A - B*i - scaleFactor_1*sim_1i - scaleFactor_2*sim_2i - ...)^2 / meas_i]
//addBG=0: no background addition
void computeBackgroundandScaling(pc_par * par)
{

  long double m_sum,s_sum[NSIMDATA],ss_sum[NSIMDATA][NSIMDATA],ms_sum[NSIMDATA],mi_sum,si_sum[NSIMDATA],i_sum,ii_sum,sum1; //sums needed to construct system of equations
  int i,j,k,l;
  lin_eq_type linEq;
  printf("\n");
  
  if(par->numFittedSimData>0)
    {
      for (i=0;i<par->numSpectra;i++)
        {
          //initialize all sums to 0
          m_sum=0.;
          mi_sum=0.;
          i_sum=0.;
          ii_sum=0.;
          sum1=0.;
          memset(s_sum,0,sizeof(s_sum));
          memset(ms_sum,0,sizeof(ms_sum));
          memset(si_sum,0,sizeof(si_sum));
          memset(ss_sum,0,sizeof(ss_sum));
          
          //construct sums  
          for (j=par->startCh[i];j<=par->endCh[i];j++)
            if(expHist[par->spectrum[i]][j]!=0)
              {
                m_sum+=fittedExpHist[par->spectrum[i]][j]/((double)expHist[par->spectrum[i]][j]);
                mi_sum+=fittedExpHist[par->spectrum[i]][j]*j/((double)expHist[par->spectrum[i]][j]);
                i_sum+=j/((double)expHist[par->spectrum[i]][j]);
                ii_sum+=j*j/((double)expHist[par->spectrum[i]][j]);
                sum1+=1./((double)expHist[par->spectrum[i]][j]);
                for (k=0;k<par->numFittedSimData;k++)
                  {
                    s_sum[k]+=fittedSimHist[k][par->spectrum[i]][j]/((double)expHist[par->spectrum[i]][j]);
                    ms_sum[k]+=fittedExpHist[par->spectrum[i]][j]*fittedSimHist[k][par->spectrum[i]][j]/((double)expHist[par->spectrum[i]][j]);
                    si_sum[k]+=fittedSimHist[k][par->spectrum[i]][j]*j/((double)expHist[par->spectrum[i]][j]);
                    for (l=0;l<par->numFittedSimData;l++)
                      ss_sum[k][l]+=fittedSimHist[k][par->spectrum[i]][j]*fittedSimHist[l][par->spectrum[i]][j]/((double)expHist[par->spectrum[i]][j]);
                  }
              }
          
          //construct system of equations (matrix/vector entries) 
          if(par->addBackground==0)
            {
              linEq.dim=par->numFittedSimData;
              for (j=0;j<par->numFittedSimData;j++)
                for (k=0;k<par->numFittedSimData;k++)
                  linEq.matrix[j][k]=ss_sum[j][k];
              for (j=0;j<par->numFittedSimData;j++)
                linEq.vector[j]=ms_sum[j];
            }
          else
            {
              linEq.dim=par->numFittedSimData+2;
              
              //top-left 4 entries
              linEq.matrix[0][0]=sum1;
              linEq.matrix[0][1]=i_sum;
              linEq.matrix[1][0]=i_sum;
              linEq.matrix[1][1]=ii_sum;
              
              //regular simulated data entires (bottom-right)
              for (j=0;j<par->numFittedSimData;j++)
                for (k=0;k<par->numFittedSimData;k++)
                  linEq.matrix[j+2][k+2]=ss_sum[j][k];
              
              //remaining entires
              for (j=0;j<par->numFittedSimData;j++)
                {     
                  linEq.matrix[0][2+j]=s_sum[j];
                  linEq.matrix[1][2+j]=si_sum[j];
                  linEq.matrix[2+j][0]=s_sum[j];
                  linEq.matrix[2+j][1]=si_sum[j];
                }
              
              linEq.vector[0]=m_sum;
              linEq.vector[1]=mi_sum;
              for (j=0;j<par->numFittedSimData;j++)
                linEq.vector[j+2]=ms_sum[j];
            }
          
          //solve system of equations and assign values
          if(!(solve_lin_eq(&linEq)==1))
            {
              if(m_sum==0)
                printf("ERROR: Experiment data (spectrum %i) has no entires in the specified fitting region!\n",par->spectrum[i]);
              else
                printf("ERROR: Could not determine background and scaling parameters!\n");
              exit(-1);
            }
          
          if(par->addBackground==0)
            {
              bgA[par->spectrum[i]]=0.;
              bgB[par->spectrum[i]]=0.;
              for (j=0;j<par->numFittedSimData;j++)
                fittedScaleFactor[j][par->spectrum[i]]=linEq.solution[j];
            }
          else
            {
              bgA[par->spectrum[i]]=linEq.solution[0];
              bgB[par->spectrum[i]]=linEq.solution[1];
              for (j=0;j<par->numFittedSimData;j++)
                fittedScaleFactor[j][par->spectrum[i]]=linEq.solution[j+2];
            }

        }
    }
  else
    printf("NOTE: All parameters are fixed, no chisq minimzation was performed.\n");
    
  //generate scaling factors for all spectra, including those that weren't fitted
  int fd=0;//counter for number of datasets which have fit (not fixed amplitude)
  int ld=-1;//index of the last dataset for which a scaling factor was determined
  for (i=0;i<par->numSimData;i++)
    {
      if(par->simDataFixedAmp[i]==0)//data was fit
        {
          for (j=0;j<par->numSpectra;j++)
            scaleFactor[i][par->spectrum[j]]=fittedScaleFactor[fd][par->spectrum[j]];
          fd++;
          ld=i;
        }
      else if(par->simDataFixedAmp[i]==2)//data is scaled relative to the previous fit data
        {
          if(ld>=0)//has a previous dataset been fit?
            for (j=0;j<par->numSpectra;j++)
              scaleFactor[i][par->spectrum[j]]=par->simDataFixedAmpValue[i]*scaleFactor[ld][par->spectrum[j]];
          ld=i;
        }
      else//data wasn't fit
        {
          for (j=0;j<par->numSpectra;j++)
            scaleFactor[i][par->spectrum[j]]=par->simDataFixedAmpValue[i];
          ld=i;
        }
    }
  
  //print parameters  
  for (i=0;i<par->numSpectra;i++)
    if(par->addBackground==0)
      {
        printf("Spectrum %i, channel %i to %i - ",par->spectrum[i],par->startCh[i],par->endCh[i]);
        for (j=0;j<par->numSimData;j++)
          if(par->simDataFixedAmp[j]==0)
            printf("Scaling factor for data from file %s: %f\n",par->simDataName[j],scaleFactor[j][par->spectrum[i]]);
          else
            printf("Scaling factor for data from file %s: %f [FIXED]\n",par->simDataName[j],scaleFactor[j][par->spectrum[i]]);
      }
    else
      {
        printf("Spectrum %i, channel %i to %i - fit linear background of form [A + B*channel], A = %0.3Lf, B = %0.3Lf\n",par->spectrum[i],par->startCh[i],par->endCh[i],bgA[par->spectrum[i]],bgB[par->spectrum[i]]);
        for (j=0;j<par->numSimData;j++)
          if(par->simDataFixedAmp[j]==0)
            printf("Scaling factor for data from file %s: %f\n",par->simDataName[j],scaleFactor[j][par->spectrum[i]]);
          else
            printf("Scaling factor for data from file %s: %f [FIXED]\n",par->simDataName[j],scaleFactor[j][par->spectrum[i]]);
      }
  
}


//function handles plotting of data, using the gnuplot_i library
void plotSpectra(pc_par * par)
{
  char str[256];
  int i,j,k;
  
  //generate plot data
  double x[par->numSpectra][par->maxNumCh];
  double yexp[par->numSpectra][par->maxNumCh];
  double ysim[NSIMDATA][par->numSpectra][par->maxNumCh];
  double ybackground[par->numSpectra][par->maxNumCh];
  double ysimsum[par->numSpectra][par->maxNumCh];
  
  for (i=0;i<par->numSpectra;i++)
    for (j=par->startCh[i];j<=par->endCh[i];j++)
      {
        x[i][j-par->startCh[i]]=(double)j;
        yexp[i][j-par->startCh[i]]=(double)expHist[par->spectrum[i]][j];
        if(par->addBackground==1)
          ybackground[i][j-par->startCh[i]]=bgA[par->spectrum[i]] + bgB[par->spectrum[i]]*j;
        ysimsum[i][j-par->startCh[i]]=ybackground[i][j-par->startCh[i]];
        for (k=0;k<par->numSimData;k++)
          {
            ysim[k][i][j-par->startCh[i]]=scaledSimHist[k][par->spectrum[i]][j];
            ysimsum[i][j-par->startCh[i]]+=scaledSimHist[k][par->spectrum[i]][j];
          }
      }

  plotOpen=1; 
  handle=gnuplot_init();
  printf("\n");
  for(i=0;i<par->numSpectra;i++)
    {
      if(par->plotStyle==1)//log-lin plot
        gnuplot_cmd(handle,"set logscale y");
      gnuplot_setstyle(handle,"steps");
      gnuplot_cmd(handle,"set xlabel 'Channel'");
      gnuplot_cmd(handle,"set ylabel 'Counts'");
      gnuplot_plot_xy(handle, x[i], yexp[i], par->endCh[i]-par->startCh[i]+1, "Experiment");
      if(par->plotOutput>1)//detailed plot
        {
          if(par->addBackground==1)//plot background
            {
              gnuplot_setstyle(handle,"lines"); //plot background as a line
              gnuplot_plot_xy(handle, x[i], ybackground[i], par->endCh[i]-par->startCh[i]+1, "Background");
              gnuplot_setstyle(handle,"steps"); //set the plot style back
            }
          for (j=0;j<par->numSimData;j++)//plot individual sim data
            {
              sprintf(str,"Simulation (%s)",par->simDataName[j]);
              gnuplot_plot_xy(handle, x[i], ysim[j][i], par->endCh[i]-par->startCh[i]+1, str);
            }
          if((par->numSimData>1)||(par->addBackground==1))//plot sum
            {
              gnuplot_setcolor(handle, "black");
              gnuplot_plot_xy(handle, x[i], ysimsum[i], par->endCh[i]-par->startCh[i]+1, "Simulation and Background(sum)");
              gnuplot_unsetcolor(handle);
            }
        }
      else //simple plot
        gnuplot_plot_xy(handle, x[i], ysimsum[i], par->endCh[i]-par->startCh[i]+1, "Simulation and Background(sum)");
      if(i!=par->numSpectra-1)//check whether we're showing the last plot
        printf("Showing plot for spectrum %i, press [ENTER] to continue...", par->spectrum[i]);
      else
        printf("Showing plot for spectrum %i, press [ENTER] to exit.", par->spectrum[i]);
      getc(stdin);
      gnuplot_resetplot(handle);
    }
  gnuplot_close(handle);
  plotOpen=0;
  
}


//function handles saving of fitted data
void saveSpectra(pc_par * par)
{
  FILE *output;
  char str[256];
  int i,j,k;
  
  //allocate arrays
  int ***outHist=allocateArrayI3(par->numSimData,par->numSpectra,S32K);
  int **bgHist=allocateArrayI2(par->numSpectra,S32K);
  
  printf("Saving scaled simulation data to output file(s)...\n");
  
  //construct arrays
  for (i=0;i<par->numSpectra;i++)
    for (j=0;j<S32K;j++)
      {
        if(par->addBackground==1)
          bgHist[i][j]=(int)(bgA[par->spectrum[i]] + bgB[par->spectrum[i]]*j);
        for (k=0;k<par->numSimData;k++)
          outHist[k][i][j]=(int)(scaledSimHist[k][par->spectrum[i]][j]);
      }

  //save arrays to .mca files  
  if(par->addBackground==1)
    {
      if((output=fopen("fit_background.mca","w"))==NULL)
        {
          printf("ERROR: Cannot open the output file fit_background.mca!\n");
          exit(-1);
        }
      for (i=0;i<par->numSpectra;i++)
        fwrite(bgHist[i],S32K*sizeof(int),1,output);
      fclose(output);
    }
  for (i=0;i<par->numSimData;i++)
    {
      sprintf(str,"fit_sim%i.mca",i);
      if((output=fopen(str,"w"))==NULL)
        {
          printf("ERROR: Cannot open the output file %s!\n",str);
          exit(-1);
        }
      for (j=0;j<par->numSpectra;j++)
        fwrite(outHist[i][j],S32K*sizeof(int),1,output);
      fclose(output);
    }
  free(outHist);
  free(bgHist);
  
}


//function run after CTRL-C, used to clean up temporary files generated
//by the plotting library
void sigint_cleanup()
{
  if(plotOpen==1)
    gnuplot_close(handle); //cleans up temporary files  
  exit(1); 
}
