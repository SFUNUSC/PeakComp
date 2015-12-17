#include "peak_comp.h"
#include "read_config.c"

int main(int argc, char *argv[])
{

  FILE *expData,*simData[NSIMDATA];

  if(argc!=2)
    {
      printf("\npeak_comp parameter_file\n");
      printf("Compares the .mca spectra designated in the parameter file specified and generates cool statistics.\n\n");
      exit(-1);
    }
  printf("\n");
  
  //initialize values
  addBackground=0;
  for (i=0;i<NSIMDATA;i++)
    for (j=0;j<NSPECT;j++)
      for (k=0;k<S32K;k++)
        {
          expHist[j][k]=0;
          simHist[i][j][k]=0;
          scaledSimHist[i][j][k]=0.;
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
  for (i=0;i<numSimData;i++)  
    if((simData[i]=fopen(simDataName[i],"r"))==NULL)
      {
        printf("ERROR: Cannot open the simulated data file %s!\n",simDataName[i]);
        exit(-1);
      }
  for (i=0;i<=endSpectrum;i++)
    if(fread(expHist[i],S32K*sizeof(int),1,expData)!=1)
      {
        printf("ERROR: Error reading file %s!\n",expDataName);
        printf("Verify that the format and number of spectra in the file are correct.\n");
        exit(-1);
      }
  for (i=0;i<numSimData;i++)    
    for (j=0;j<=endSpectrum;j++)
      if(fread(simHist[i][j],S32K*sizeof(int),1,simData[i])!=1)
        {
          printf("ERROR: Error reading file %s!\n",simDataName[i]);
          printf("Verify that the format and number of spectra in the file are correct.\n");
          exit(-1);
        }
  fclose(expData);
  for (i=0;i<numSimData;i++) 
    fclose(simData[i]);
  printf("Spectra read in...\n");
  
  if(addBackground==0)//no background addition
    {
      //calculate integrals of data in each spectrum and scale the experiment and sim to each other
      
      for (i=0;i<numSimData;i++)
        {
          printf("Comparing experiment data to simulated data from file: %s\n",simDataName[i]);
          for (j=0;j<numSpectra;j++)
            {
              printf("Calculating integrals of data in spectrum %i...\n",spectrum[j]);
              expInt=0.;
              simInt[i]=0.;
            
              for (k=startCh[j];k<=endCh[j];k++)
                {
                  expInt+=(double)expHist[spectrum[j]][k];
                  simInt[i]+=(double)simHist[i][spectrum[j]][k];
                }
              scaleFactor[i][spectrum[j]]=expInt/simInt[i];
              printf("Experiment: %1.0f, Simulated: %1.0f\n",expInt,simInt[i]);
              printf("Scaling simulated data from file by a factor of: %f\n",scaleFactor[i][spectrum[j]]);
            }
        }
        
      for (i=0;i<numSimData;i++)  
        for (j=0;j<numSpectra;j++)
          for (k=startCh[j];k<=endCh[j];k++)
            {
              scaledSimHist[i][spectrum[j]][k]=scaleFactor[i][spectrum[j]]*simHist[i][spectrum[j]][k];
            }
    
      compareSpectra();
    }
  else if(addBackground==1)//linear background addition
    {
    
      computeLinearBackground(numSimData);//get background coefficients and scaling factor
      
      //scale simulated data
      for (i=0;i<numSimData;i++)
        for (j=0;j<numSpectra;j++)
          for (k=startCh[j];k<=endCh[j];k++)
            scaledSimHist[i][spectrum[j]][k]=scaleFactor[i][spectrum[j]]*simHist[i][spectrum[j]][k];
      //add background to simulated data
      for (i=0;i<numSimData;i++)
        for (j=0;j<numSpectra;j++)
          for (k=startCh[j];k<=endCh[j];k++)
            scaledSimHist[i][spectrum[j]][k]=scaledSimHist[i][spectrum[j]][k] + bgA[spectrum[j]] + bgB[spectrum[j]]*k;
      
      compareSpectra();
    }
      
  //print output
  printf("\nCOMPARISON DATA\n---------------\n");
  for (i=0;i<numSpectra;i++)
    printf("chisq (spectrum %i): %f\n",spectrum[i],spectChisq[i]);
  printf("chisq (total): %f\n",chisq);
  printf("number of bins: %i\n",numBinsUsed);
  printf("chisq (total) / number of bins: %f\n",chisq/(numBinsUsed));
  
  if(plotOutput>=1)
    plotSpectra();

  return 0; //great success
}

//function compares spectra and gets chisq and other stats
void compareSpectra()
{

  //initialize values
  int binsSkipped=0;
  numBinsUsed=0;
  double sumSimValue=0;
  
  //compute chisq for data in the spectra
  chisq=0;
  for (i=0;i<numSpectra;i++)
    {
      spectChisq[i]=0;
      for (j=startCh[i];j<=endCh[i];j++)
        {  
          if(expHist[spectrum[i]][j]!=0)//avoid dividing by zero
            {
              //get the sum of all experimental data in the given bin 
              sumSimValue=0;
              for (k=0;k<numSimData;k++)
                sumSimValue+=scaledSimHist[k][spectrum[i]][j];
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
      printf("Warning: some of the bins in the experiment data have values of zero.  These have been skipped when calculating chisq.\n");
      printf("Bins skipped: %i.\n",binsSkipped);
    }
  for (i=0;i<numSpectra;i++)
    numBinsUsed += (endCh[i]-startCh[i]+1);
  numBinsUsed -= binsSkipped;
}

//function computes background coefficients and scaling factors
//by analytically minimizing chisq for the expression:
//chisq=sum_i[(meas_i - A - B*i - scaleFactor_1*sim_1i - scaleFactor_2*sim_2i - ...)^2 / meas_i]
void computeLinearBackground(int numSimData)
{
  long double m_sum,s_sum[NSIMDATA],ss_sum[NSIMDATA][NSIMDATA],ms_sum[NSIMDATA],mi_sum,si_sum[NSIMDATA],i_sum,ii_sum,sum1; //sums needed to construct system of equations
  lin_eq_type linEq;
  printf("\n");  
  
  for (i=0;i<numSpectra;i++)
    {
      //initialize all sums to 0
      m_sum=0.;
      mi_sum=0.;
      i_sum=0.;
      ii_sum=0.;
      sum1=0.;
      for (j=0;j<numSimData;j++)
        {
          s_sum[j]=0.;
          ms_sum[j]=0.;
          si_sum[j]=0.;
          for (k=0;k<numSimData;k++)
            ss_sum[j][k]=0.;
        }
      
      //construct sums  
      for (j=startCh[i];j<=endCh[i];j++)
        if(expHist[spectrum[i]][j]!=0)
          {
            m_sum+=expHist[spectrum[i]][j]/((double)expHist[spectrum[i]][j]);
            mi_sum+=expHist[spectrum[i]][j]*j/((double)expHist[spectrum[i]][j]);
            i_sum+=j/((double)expHist[spectrum[i]][j]);
            ii_sum+=j*j/((double)expHist[spectrum[i]][j]);
            sum1+=1./((double)expHist[spectrum[i]][j]);
            for (k=0;k<numSimData;k++)
              {
                s_sum[k]+=simHist[k][spectrum[i]][j]/((double)expHist[spectrum[i]][j]);
                ms_sum[k]+=expHist[spectrum[i]][j]*simHist[k][spectrum[i]][j]/((double)expHist[spectrum[i]][j]);
                si_sum[k]+=simHist[k][spectrum[i]][j]*j/((double)expHist[spectrum[i]][j]);
                for (l=0;l<numSimData;l++)
                  ss_sum[k][l]+=simHist[k][spectrum[i]][j]*simHist[l][spectrum[i]][j]/((double)expHist[spectrum[i]][j]);
              }
          }
      
      //construct system of equations (matrix/vector entries)
      linEq.dim=numSimData+2;
      for (j=0;j<numSimData;j++)
        {
          linEq.matrix[j][j]=ss_sum[j][j];
          if((j+1)<numSimData)
            {
              linEq.matrix[j][j+1]=ss_sum[j][j+1];
              linEq.matrix[j+1][j]=ss_sum[j+1][j];
            }
          linEq.matrix[linEq.dim-2][j]=s_sum[j];
          linEq.matrix[linEq.dim-1][j]=si_sum[j];
          linEq.matrix[j][linEq.dim-2]=s_sum[j];
          linEq.matrix[j][linEq.dim-1]=si_sum[j];
          //printf("%i of %i\n",j,linEq.dim-1);
        }
      linEq.matrix[linEq.dim-2][linEq.dim-2]=sum1;
      linEq.matrix[linEq.dim-2][linEq.dim-1]=i_sum;
      linEq.matrix[linEq.dim-1][linEq.dim-2]=i_sum;
      linEq.matrix[linEq.dim-1][linEq.dim-1]=ii_sum;
      for (j=0;j<numSimData;j++)
        linEq.vector[j]=ms_sum[j];
      linEq.vector[linEq.dim-2]=m_sum;
      linEq.vector[linEq.dim-1]=mi_sum;
      
      //solve system of equations and assign values
      if(!(solve_lin_eq(&linEq)==1))
        {
          printf("ERROR: Could not determine background and scaling parameters!\n");
          exit(-1);
        }
      for (j=0;j<numSimData;j++)  
        scaleFactor[j][spectrum[i]]=linEq.solution[j];
      bgA[spectrum[i]]=linEq.solution[linEq.dim-2];
      bgB[spectrum[i]]=linEq.solution[linEq.dim-1];
      
      printf("Spectrum %i: fit linear background of form [A + B*channel], A = %0.3Lf, B = %0.3Lf\n",spectrum[i],bgA[spectrum[i]],bgB[spectrum[i]]);
      for (j=0;j<numSimData;j++)
        printf("Scaling factor for data from file %s: %f\n",simDataName[j],scaleFactor[j][spectrum[i]]);
          
       
    }
     
  
}

//function handles plotting of data, using the gnuplot_i library
void plotSpectra()
{

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
        ysimsum[i][j-startCh[i]]=0.;
        for (k=0;k<numSimData;k++)
          {
            ysim[k][i][j-startCh[i]]=scaledSimHist[k][spectrum[i]][j] - ybackground[i][j-startCh[i]];
            ysimsum[i][j-startCh[i]]+=scaledSimHist[k][spectrum[i]][j];
          }
      }

  gnuplot_ctrl *handle;    
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
            gnuplot_plot_xy(handle, x[i], ysimsum[i], endCh[i]-startCh[i]+1, "Simulation and Background(sum)");
        }
      else //simple plot
        gnuplot_plot_xy(handle, x[i], ysimsum[i], endCh[i]-startCh[i]+1, "Simulation and Background(sum)");
      printf("Showing plot for spectrum %i, press [ENTER] to continue...", spectrum[i]);
      getc(stdin);
      gnuplot_resetplot(handle);
    }
  gnuplot_close(handle);
  
}
