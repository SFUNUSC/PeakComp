//function computes background coefficients and scaling factors
//by analytically minimizing chisq for the expression:
//chisq=sum_i[(meas_i - A - B*i - scaleFactor_1*sim_1i - scaleFactor_2*sim_2i - ...)^2 / meas_i]
void computeBackgroundandScaling(const par * p, const data * d, fitpar * fp)
{

  long double m_sum,s_sum[NSIMDATA],ss_sum[NSIMDATA][NSIMDATA],ms_sum[NSIMDATA],
              mi_sum,mii_sum,si_sum[NSIMDATA],sii_sum[NSIMDATA],i_sum,ii_sum,
              iii_sum,iiii_sum,sum1; //sums needed to construct system of equations
  long double ind;
  double **fittedScaleFactor=allocateArrayD2(p->numFittedSimData,p->numSpectra);
  int i,j,k,l;
  lin_eq_type linEq;
  if(p->verbose>=0) printf("\n");
  
  if(p->numFittedSimData>0)
    {
      for (i=0;i<p->numSpectra;i++)
        {
          //initialize all sums to 0
          m_sum=0.;
          mi_sum=0.;
          mii_sum=0.;
          i_sum=0.;
          ii_sum=0.;
          iii_sum=0.;
          iiii_sum=0.;
          sum1=0.;
          memset(s_sum,0,sizeof(s_sum));
          memset(ms_sum,0,sizeof(ms_sum));
          memset(si_sum,0,sizeof(si_sum));
          memset(sii_sum,0,sizeof(sii_sum));
          memset(ss_sum,0,sizeof(ss_sum));
          
          
          //construct sums
          for (j=p->startCh[i];j<=p->endCh[i];j++)
            if(d->expHist[p->spectrum[i]][j]!=0)
              {
                m_sum+=d->fittedExpHist[p->spectrum[i]][j]/((double)d->expHist[p->spectrum[i]][j]);
                for (k=0;k<p->numFittedSimData;k++)
                  {
                    ms_sum[k]+=d->fittedExpHist[p->spectrum[i]][j]*(double)d->fittedSimHist[k][p->spectrum[i]][j]/((double)d->expHist[p->spectrum[i]][j]);//cast to double in numerator needed to prevent overflow
                    for (l=0;l<p->numFittedSimData;l++)
                      ss_sum[k][l]+=(double)d->fittedSimHist[k][p->spectrum[i]][j]*(double)d->fittedSimHist[l][p->spectrum[i]][j]/((double)d->expHist[p->spectrum[i]][j]);
                  }
              }
          if(p->fitAddBackground[i]>=1)
            for (j=p->startCh[i];j<=p->endCh[i];j++)
              if(d->expHist[p->spectrum[i]][j]!=0)
                {
                  ind=(long double)j;  
                  mi_sum+=d->fittedExpHist[p->spectrum[i]][j]*ind/((double)d->expHist[p->spectrum[i]][j]);
                  i_sum+=ind/((double)d->expHist[p->spectrum[i]][j]);
                  ii_sum+=ind*ind/((double)d->expHist[p->spectrum[i]][j]);
                  sum1+=1./((double)d->expHist[p->spectrum[i]][j]);
                  for (k=0;k<p->numFittedSimData;k++)
                    {
                      s_sum[k]+=d->fittedSimHist[k][p->spectrum[i]][j]/((double)d->expHist[p->spectrum[i]][j]);
                      si_sum[k]+=d->fittedSimHist[k][p->spectrum[i]][j]*ind/((double)d->expHist[p->spectrum[i]][j]);
                    }
                }
          if(p->fitAddBackground[i]>=2)
            for (j=p->startCh[i];j<=p->endCh[i];j++)
              if(d->expHist[p->spectrum[i]][j]!=0)
                {
                  ind=(long double)j;
                  mii_sum+=d->fittedExpHist[p->spectrum[i]][j]*ind*ind/((double)d->expHist[p->spectrum[i]][j]);
                  iii_sum+=ind*ind*ind/((double)d->expHist[p->spectrum[i]][j]);
                  iiii_sum+=ind*ind*ind*ind/((double)d->expHist[p->spectrum[i]][j]);
                  for (k=0;k<p->numFittedSimData;k++)
                    sii_sum[k]+=d->fittedSimHist[k][p->spectrum[i]][j]*ind*ind/((double)d->expHist[p->spectrum[i]][j]);
                }
          
          //construct system of equations (matrix/vector entries) 
          if(p->fitAddBackground[i]==0)
            {
              linEq.dim=p->numFittedSimData;
              for (j=0;j<p->numFittedSimData;j++)
                for (k=0;k<p->numFittedSimData;k++)
                  linEq.matrix[j][k]=ss_sum[j][k];
              for (j=0;j<p->numFittedSimData;j++)
                linEq.vector[j]=ms_sum[j];
            }
          else if(p->fitAddBackground[i]==1)
            {
              linEq.dim=p->numFittedSimData+2;
              
              //top-left 4 entries
              linEq.matrix[0][0]=sum1;
              linEq.matrix[0][1]=i_sum;
              linEq.matrix[1][0]=i_sum;
              linEq.matrix[1][1]=ii_sum;
              
              //regular simulated data entires (bottom-right)
              for (j=0;j<p->numFittedSimData;j++)
                for (k=0;k<p->numFittedSimData;k++)
                  linEq.matrix[2+j][2+k]=ss_sum[j][k];
              
              //remaining entires
              for (j=0;j<p->numFittedSimData;j++)
                {     
                  linEq.matrix[0][2+j]=s_sum[j];
                  linEq.matrix[1][2+j]=si_sum[j];
                  linEq.matrix[2+j][0]=s_sum[j];
                  linEq.matrix[2+j][1]=si_sum[j];
                }
              
              linEq.vector[0]=m_sum;
              linEq.vector[1]=mi_sum;
              for (j=0;j<p->numFittedSimData;j++)
                linEq.vector[j+2]=ms_sum[j];
            }
          else if(p->fitAddBackground[i]==2)
            {
              linEq.dim=p->numFittedSimData+3;
              
              //top-left 9 entries
              linEq.matrix[0][0]=sum1;
              linEq.matrix[0][1]=i_sum;
              linEq.matrix[0][2]=ii_sum;
              linEq.matrix[1][0]=i_sum;
              linEq.matrix[1][1]=ii_sum;
              linEq.matrix[1][2]=iii_sum;
              linEq.matrix[2][0]=ii_sum;
              linEq.matrix[2][1]=iii_sum;
              linEq.matrix[2][2]=iiii_sum;
              
              //regular simulated data entires (bottom-right)
              for (j=0;j<p->numFittedSimData;j++)
                for (k=0;k<p->numFittedSimData;k++)
                  linEq.matrix[3+j][3+k]=ss_sum[j][k];
              
              //remaining entires
              for (j=0;j<p->numFittedSimData;j++)
                {     
                  linEq.matrix[0][3+j]=s_sum[j];
                  linEq.matrix[1][3+j]=si_sum[j];
                  linEq.matrix[2][3+j]=sii_sum[j];
                  linEq.matrix[3+j][0]=s_sum[j];
                  linEq.matrix[3+j][1]=si_sum[j];
                  linEq.matrix[3+j][2]=sii_sum[j];
                }
              
              linEq.vector[0]=m_sum;
              linEq.vector[1]=mi_sum;
              linEq.vector[2]=mii_sum;
              for (j=0;j<p->numFittedSimData;j++)
                linEq.vector[j+3]=ms_sum[j];
            }
          
          //solve system of equations and assign values
          if(!(solve_lin_eq(&linEq)==1))
            {
              if(m_sum==0)
                printf("ERROR: Experiment data has no entries in the specified fitting region (spectrum %i, channel %i to %i).\n",p->spectrum[i],p->startCh[i],p->endCh[i]);
              else
                printf("ERROR: Could not determine background and scaling parameters for spectrum %i, channel %i to %i.\n",p->spectrum[i],p->startCh[i],p->endCh[i]);
              exit(-1);
            }
          
          if(p->fitAddBackground[i]==0)
            {
              fp->bgA[i]=0.;
              fp->bgB[i]=0.;
              fp->bgC[i]=0.;
              for (j=0;j<p->numFittedSimData;j++)
                fittedScaleFactor[j][i]=linEq.solution[j];
            }
          else if(p->fitAddBackground[i]==1)
            {
              fp->bgA[i]=linEq.solution[0];
              fp->bgB[i]=linEq.solution[1];
              fp->bgC[i]=0.;
              for (j=0;j<p->numFittedSimData;j++)
                fittedScaleFactor[j][i]=linEq.solution[j+2];
            }
          else if(p->fitAddBackground[i]==2)
            {
              fp->bgA[i]=linEq.solution[0];
              fp->bgB[i]=linEq.solution[1];
              fp->bgC[i]=linEq.solution[2];
              for (j=0;j<p->numFittedSimData;j++)
                fittedScaleFactor[j][i]=linEq.solution[j+3];
            }

        }
    }
  else if(p->verbose>=0)
    printf("NOTE: All scaling parameters are fixed, no chisq minimzation was performed.\n\n");
    
  //generate scaling factors for all spectra, including those that weren't fitted
  int fd=0;//counter for number of datasets which have fit (not fixed amplitude)
  int ld=-1;//index of the last dataset for which a scaling factor was determined
  for (i=0;i<p->numSimData;i++)
    {
      if(p->simDataFixedAmp[i]==0)//data was fit
        {
          for (j=0;j<p->numSpectra;j++)
            fp->scaleFactor[i][j]=fittedScaleFactor[fd][j];
          fd++;
          ld=i;
        }
      else if(p->simDataFixedAmp[i]==2)//data is scaled relative to the previous fit data
        {
          if(ld>=0)//has a previous dataset been fit?
            for (j=0;j<p->numSpectra;j++)
              fp->scaleFactor[i][j]=p->simDataFixedAmpValue[i]*fp->scaleFactor[ld][j];
          ld=i;
        }
      else//data wasn't fit
        {
          for (j=0;j<p->numSpectra;j++)
            fp->scaleFactor[i][j]=p->simDataFixedAmpValue[i];
          ld=i;
        }
    }
  //check for fixed background and generate background parameters if needed
  for (i=0;i<p->numSimData;i++)
    if(p->fixBG[i]==1)
      {
        if(p->addBackground>=1)
          {
            fp->bgA[i]=p->fixedBGPar[i][0];
            fp->bgB[i]=p->fixedBGPar[i][1];
          }
        if(p->addBackground>=2)
          {
            fp->bgC[i]=p->fixedBGPar[i][2];
          }
      }
  
  free(fittedScaleFactor);
  
  //print fit data
  if(p->verbose>=0)
    {
      printf("FIT DATA\n--------\n");
      for (i=0;i<p->numSpectra;i++)
        {
          if(p->addBackground==0)
            printf("Spectrum %i, channel %i to %i:\n",p->spectrum[i],p->startCh[i],p->endCh[i]);
          else if((p->addBackground==1)&&(p->fitAddBackground[i]==p->addBackground))
            printf("Spectrum %i, channel %i to %i:\nFit linear background of form [A + B*channel],\nA = %0.5LE, B = %0.5LE\n",p->spectrum[i],p->startCh[i],p->endCh[i],fp->bgA[i],fp->bgB[i]);
          else if(p->addBackground==1)
            printf("Spectrum %i, channel %i to %i:\nUsing linear background of form [A + B*channel],\nA = %0.5LE [FIXED], B = %0.5LE [FIXED]\n",p->spectrum[i],p->startCh[i],p->endCh[i],fp->bgA[i],fp->bgB[i]);
          else if((p->addBackground==2)&&(p->fitAddBackground[i]==p->addBackground))
            printf("Spectrum %i, channel %i to %i:\nFit quadratic background of form [A + B*channel + C*(channel^2)],\nA = %0.5LE, B = %0.5LE, C = %0.5LE\n",p->spectrum[i],p->startCh[i],p->endCh[i],fp->bgA[i],fp->bgB[i],fp->bgC[i]);
          else if(p->addBackground==2)
            printf("Spectrum %i, channel %i to %i:\nUsing quadratic background of form [A + B*channel + C*(channel^2)],\nA = %0.5LE [FIXED], B = %0.5LE [FIXED], C = %0.5LE [FIXED]\n",p->spectrum[i],p->startCh[i],p->endCh[i],fp->bgA[i],fp->bgB[i],fp->bgC[i]);
          for (j=0;j<p->numSimData;j++)
            {
              if(p->simDataFixedAmp[j]==0)
                printf("Scaling factor for data from file %s: %f\n",p->simDataName[j],fp->scaleFactor[j][i]);
              else
                printf("Scaling factor for data from file %s: %f [FIXED]\n",p->simDataName[j],fp->scaleFactor[j][i]);
            }
          printf("\n");
        }
    }

}

void applyBackgroundandScaling(const par * p, const fitpar * fp, const data * d, fitdata * fd)
{
  int i,j,k;
  //scale simulated data
  for (i=0;i<p->numSimData;i++)
    for (j=0;j<p->numSpectra;j++)
      for (k=0;k<S32K;k++)
        fd->scaledSimHist[i][j][k]=fp->scaleFactor[i][j]*d->simHist[i][p->spectrum[j]][k];//if the same spectrum in the simulated data 
                                                                                          //is compared over multiple windows, 
                                                                                          //generate different scaled data each time
  //generate background data
  for (i=0;i<p->numSpectra;i++)
    for (j=0;j<S32K;j++)
      if(p->addBackground>=1)
        fd->bgHist[i][j]=fp->bgA[i] + fp->bgB[i]*j + fp->bgC[i]*j*j;
}
