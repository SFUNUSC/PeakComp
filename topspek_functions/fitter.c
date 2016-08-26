//forward declarations
void generateSums(const par*, const data*, fitsum*, const int, const int);
void solveFitEq(const par*, const data*, const fitsum*, const int, const int, lin_eq_type*);

//function computes background coefficients and scaling factors
//by analytically minimizing chisq for the expression:
//chisq=sum_i[(meas_i - A - B*i - scaleFactor_1*sim_1i - scaleFactor_2*sim_2i - ...)^2 / meas_i]
void computeBackgroundandScaling(const par * p, const data * d, fitpar * fp)
{

  fitsum *fs=(fitsum*)malloc(sizeof(fitsum));
  double **fittedScaleFactor=allocateArrayD2(p->numSimData,p->numSpectra);//1st index will only be filled up to number of fitted data in each spectrum
  int i,j;
  lin_eq_type linEq;
  if(p->verbose>=0) printf("\n");
  
  if(p->commonScaling==1)
    {
      //generate sums over all spectra and then solve simultaneously
      for (i=0;i<p->numSpectra;i++)
        generateSums(p,d,fs,i,p->fitAddBackground[i]);
      solveFitEq(p,d,fs,0,p->fitAddBackground[0],&linEq);
      
      for (i=0;i<p->numSpectra;i++)
        {
          if(p->fitAddBackground[0]==0)
            {
              fp->bgA[i]=0.;
              fp->bgB[i]=0.;
              fp->bgC[i]=0.;
              for (j=0;j<p->numFittedSimData[i];j++)
                fittedScaleFactor[j][i]=linEq.solution[j];
            }
          else if(abs(p->fitAddBackground[0])==1)
            {
              fp->bgA[i]=linEq.solution[0];
              fp->bgB[i]=0.;
              fp->bgC[i]=0.;
              for (j=0;j<p->numFittedSimData[i];j++)
                fittedScaleFactor[j][i]=linEq.solution[j+1];
            }
          else if(p->fitAddBackground[0]==2)
            {
              fp->bgA[i]=linEq.solution[0];
              fp->bgB[i]=linEq.solution[1];
              fp->bgC[i]=0.;
              for (j=0;j<p->numFittedSimData[i];j++)
                fittedScaleFactor[j][i]=linEq.solution[j+2];
            }
          else if(p->fitAddBackground[0]==3)
            {
              fp->bgA[i]=linEq.solution[0];
              fp->bgB[i]=linEq.solution[1];
              fp->bgC[i]=linEq.solution[2];
              for (j=0;j<p->numFittedSimData[i];j++)
                fittedScaleFactor[j][i]=linEq.solution[j+3];
            }
        }
    }
  else
    {
      for (i=0;i<p->numSpectra;i++)
      	{
      		if(p->numFittedSimData[i]>0)
      	  	{        
	      	    //generate sums and solve spectrum-by-spectrum
    	      	generateSums(p,d,fs,i,p->fitAddBackground[i]);//generate sums for the spectrum
		          solveFitEq(p,d,fs,i,p->fitAddBackground[i],&linEq);//solve for the spectrum
          
		          if(p->fitAddBackground[i]==0)//no background fitting
	  	          {
  	  	          fp->bgA[i]=0.;
	        	      fp->bgB[i]=0.;
	          	    fp->bgC[i]=0.;
            	  	for (j=0;j<p->numFittedSimData[i];j++)
              	  	fittedScaleFactor[j][i]=linEq.solution[j];
		            }
		          else if(abs(p->fitAddBackground[i])==1)//constant background
		            {
    		          fp->bgA[i]=linEq.solution[0];
        		      fp->bgB[i]=0.;
            		  fp->bgC[i]=0.;
				          for (j=0;j<p->numFittedSimData[i];j++)
				            fittedScaleFactor[j][i]=linEq.solution[j+1];
				        }
				      else if(p->fitAddBackground[i]==2)//linear background
				        {
				          fp->bgA[i]=linEq.solution[0];
				          fp->bgB[i]=linEq.solution[1];
				          fp->bgC[i]=0.;
				          for (j=0;j<p->numFittedSimData[i];j++)
				            fittedScaleFactor[j][i]=linEq.solution[j+2];
				        }
				      else if(p->fitAddBackground[i]==3)//quadratic background
				        {
				          fp->bgA[i]=linEq.solution[0];
				          fp->bgB[i]=linEq.solution[1];
				          fp->bgC[i]=linEq.solution[2];
				          for (j=0;j<p->numFittedSimData[i];j++)
				            fittedScaleFactor[j][i]=linEq.solution[j+3];
				        }
		    		}
				  		else if(p->verbose>=0)
								printf("Spectrum %i: All scaling parameters are fixed, no chisq minimzation was performed.\n",i);
				}
  }
  
  free(fs);
    
  //generate scaling factors for all spectra, including those that weren't fitted
  int fd[NSPECT];//counter for number of datasets which have fit (not fixed amplitude)
  int ld[NSPECT];//index of the last dataset for which a scaling factor was determined
  for (j=0;j<p->numSpectra;j++)
  	ld[j]=-1;
  	
  for (i=0;i<p->numSimData;i++)
  	for (j=0;j<p->numSpectra;j++)
    	{
      	if(p->simDataFixedAmp[i][j]==0)//data was fit
        	{
            fp->scaleFactor[i][j]=fittedScaleFactor[fd[j]][j];
          	fd[j]++;
          	ld[j]=i;
        	}
      	else if(p->simDataFixedAmp[i][j]==1)//data wasn't fit (fixed amplitude)
        	{
            fp->scaleFactor[i][j]=p->simDataFixedAmpValue[i][j];
          	ld[j]=i;
        	}
				else if(p->simDataFixedAmp[i][j]==2)//data is scaled relative to the previous fit data
        	{
          	if(ld[j]>=0)//has a previous dataset been fit?
              fp->scaleFactor[i][j]=p->simDataFixedAmpValue[i][j]*fp->scaleFactor[ld[j]][j];
          	ld[j]=i;
        	}
    	}

  //check for fixed background and generate background parameters if needed
  for (i=0;i<p->numSimData;i++)
    if(p->fixBG[i]==1)
      {
        if(p->addBackground>=2)
          {
            fp->bgA[i]=p->fixedBGPar[i][0];
            fp->bgB[i]=p->fixedBGPar[i][1];
          }
        if(p->addBackground>=3)
          {
            fp->bgC[i]=p->fixedBGPar[i][2];
          }
      }
  
  free(fittedScaleFactor);
  
}

//checks for whether a re-fit should be performed
//return values:
//0 - don't re-fit
//1 - do re-fit 
int checkForRefit(par * p, fitpar * fp)
{
	int refit=0;
	int i,j,k;
	if(p->fitAddBackground[0]==-1)//force positive constant background
		{
			for (i=0;i<p->numSimData;i++)
				if(p->fixBG[i]==0)//make sure background isn't fixed
					if(fp->bgA[i]<0.)
						{
							p->fixBG[i]=0;//fix the background
							p->fixedBGPar[i][0]=0;//fix it to 0
							p->fixedBGPar[i][1]=0;
							p->fixedBGPar[i][2]=0;
							refit=1;
						}
		}
	if(p->forcePositiveS==1)//force positive scaling factors
  	{
  		if(refit==0)//check that refit hasn't previously been called (only refit one thing at a time)
				for (i=0;i<p->numSimData;i++)
					for (j=0;j<p->numSpectra;j++)
						if(p->simDataFixedAmp[i][j]==0)//if data was fit
							if(fp->scaleFactor[i][j]<0.)
								{
									//printf("fp->scaleFactor[%i][%i] = %f\n",i,j,fp->scaleFactor[i][j]);
									p->simDataFixedAmp[i][j]=1;//fix scaling
									p->simDataFixedAmpValue[i][j]=0.;//fix to 0
									//deal with data scaled relative to this data
									k=i+1;
									while(p->simDataFixedAmp[k][j]==2)
										{
											p->simDataFixedAmp[k][j]=1;//fix scaling
											p->simDataFixedAmpValue[k][j]=0.;//fix to 0
											k++;
										}
									p->numFittedSimData[j]-=1;
									refit=1;
								}
  	}
  if((p->verbose>=0)&&(refit==1))
  	printf("Some fit values were out of specified bounds.  Refitting...\n\n");
  return refit;
}


//functions prints parameter values from fit
void printFitData(const par * p, const fitpar * fp)
{
	int i,j;
	
  if(p->verbose>=0)
    {
      printf("FIT DATA\n--------\n");
      for (i=0;i<p->numSpectra;i++)
        {
          if(p->addBackground==0)
            printf("Spectrum %i, channel %i to %i:\n",p->spectrum[i],p->startCh[i],p->endCh[i]);
          else if((abs(p->addBackground)==1)&&(p->fixBG[i]==0))
            printf("Spectrum %i, channel %i to %i:\nFit constant background of amplitude A = %0.5LE\n",p->spectrum[i],p->startCh[i],p->endCh[i],fp->bgA[i]);
          else if(abs(p->addBackground)==1)
            printf("Spectrum %i, channel %i to %i:\nUsing constant background of amplitude A = %0.5LE [FIXED]\n",p->spectrum[i],p->startCh[i],p->endCh[i],fp->bgA[i]);
          else if((p->addBackground==2)&&(p->fixBG[i]==0))
            printf("Spectrum %i, channel %i to %i:\nFit linear background of form [A + B*channel],\nA = %0.5LE, B = %0.5LE\n",p->spectrum[i],p->startCh[i],p->endCh[i],fp->bgA[i],fp->bgB[i]);
          else if(p->addBackground==2)
            printf("Spectrum %i, channel %i to %i:\nUsing linear background of form [A + B*channel],\nA = %0.5LE [FIXED], B = %0.5LE [FIXED]\n",p->spectrum[i],p->startCh[i],p->endCh[i],fp->bgA[i],fp->bgB[i]);
          else if((p->addBackground==3)&&(p->fixBG[i]==0))
            printf("Spectrum %i, channel %i to %i:\nFit quadratic background of form [A + B*channel + C*(channel^2)],\nA = %0.5LE, B = %0.5LE, C = %0.5LE\n",p->spectrum[i],p->startCh[i],p->endCh[i],fp->bgA[i],fp->bgB[i],fp->bgC[i]);
          else if(p->addBackground==3)
            printf("Spectrum %i, channel %i to %i:\nUsing quadratic background of form [A + B*channel + C*(channel^2)],\nA = %0.5LE [FIXED], B = %0.5LE [FIXED], C = %0.5LE [FIXED]\n",p->spectrum[i],p->startCh[i],p->endCh[i],fp->bgA[i],fp->bgB[i],fp->bgC[i]);
          for (j=0;j<p->numSimData;j++)
            {
              if(p->simDataFixedAmp[j][i]==0)
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
      if(abs(p->addBackground)>=1)
        fd->bgHist[i][j]=fp->bgA[i] + fp->bgB[i]*j + fp->bgC[i]*j*j;
}


//generates sums for fitting routine for the given spectrum
//specNum is 0 indexed
void generateSums(const par * p, const data * d, fitsum * fs, const int specNum, const int bgMode)
{
  long double ind;
  int j,k,l;
  
  //initialize sums
  memset(fs,0,sizeof(fitsum));
  
  /*//debug
  printf("p->numFittedSimData[%i]: %i\n",specNum,p->numFittedSimData[specNum]);
  printf("p->startCh[%i]: %i, p->endCh[%i]: %i\n",specNum,p->startCh[specNum],specNum,p->endCh[specNum]);
  printf("d->fittedSimHist[0][%i][%i]: %f\n",specNum,p->startCh[specNum],d->fittedSimHist[0][specNum][p->startCh[specNum]]);*/
  
  //construct sums
  for (j=p->startCh[specNum];j<=p->endCh[specNum];j++)
    if(d->expHist[p->spectrum[specNum]][j]!=0)
      {
        fs->m_sum+=d->fittedExpHist[specNum][j]/((double)d->expHist[p->spectrum[specNum]][j]);
        for (k=0;k<p->numFittedSimData[specNum];k++)
          {
            fs->ms_sum[k]+=d->fittedExpHist[specNum][j]*(double)d->fittedSimHist[k][specNum][j]/((double)d->expHist[p->spectrum[specNum]][j]);//cast to double in numerator needed to prevent overflow
            for (l=0;l<p->numFittedSimData[specNum];l++)
              fs->ss_sum[k][l]+=(double)d->fittedSimHist[k][specNum][j]*(double)d->fittedSimHist[l][specNum][j]/((double)d->expHist[p->spectrum[specNum]][j]);
          }
      }
  if(abs(bgMode)>=1)
    for (j=p->startCh[specNum];j<=p->endCh[specNum];j++)
      if(d->expHist[p->spectrum[specNum]][j]!=0)
        {
          fs->sum1+=1./((double)d->expHist[p->spectrum[specNum]][j]);
          for (k=0;k<p->numFittedSimData[specNum];k++)
            fs->s_sum[k]+=d->fittedSimHist[k][specNum][j]/((double)d->expHist[p->spectrum[specNum]][j]);
        }
  if(abs(bgMode)>=2)
    for (j=p->startCh[specNum];j<=p->endCh[specNum];j++)
      if(d->expHist[p->spectrum[specNum]][j]!=0)
        {
          ind=(long double)j;  
          fs->mi_sum+=d->fittedExpHist[specNum][j]*ind/((double)d->expHist[p->spectrum[specNum]][j]);
          fs->i_sum+=ind/((double)d->expHist[p->spectrum[specNum]][j]);
          fs->ii_sum+=ind*ind/((double)d->expHist[p->spectrum[specNum]][j]);
          for (k=0;k<p->numFittedSimData[specNum];k++)
            {
              fs->si_sum[k]+=d->fittedSimHist[k][specNum][j]*ind/((double)d->expHist[p->spectrum[specNum]][j]);
            }
        }
  if(abs(bgMode)>=3)
    for (j=p->startCh[specNum];j<=p->endCh[specNum];j++)
      if(d->expHist[p->spectrum[specNum]][j]!=0)
        {
          ind=(long double)j;
          fs->mii_sum+=d->fittedExpHist[specNum][j]*ind*ind/((double)d->expHist[p->spectrum[specNum]][j]);
          fs->iii_sum+=ind*ind*ind/((double)d->expHist[p->spectrum[specNum]][j]);
          fs->iiii_sum+=ind*ind*ind*ind/((double)d->expHist[p->spectrum[specNum]][j]);
          for (k=0;k<p->numFittedSimData[specNum];k++)
            fs->sii_sum[k]+=d->fittedSimHist[k][specNum][j]*ind*ind/((double)d->expHist[p->spectrum[specNum]][j]);
        }
}

//given a set of sums, set up and solve the fit equation
void solveFitEq(const par * p, const data * d, const fitsum * fs, const int specNum, const int bgMode, lin_eq_type * linEq)
{
  int j,k;
  
  //construct system of equations (matrix/vector entries) 
  if(bgMode==0)//no background
    {
      linEq->dim=p->numFittedSimData[specNum];
      for (j=0;j<p->numFittedSimData[specNum];j++)
        for (k=0;k<p->numFittedSimData[specNum];k++)
          linEq->matrix[j][k]=fs->ss_sum[j][k];
      for (j=0;j<p->numFittedSimData[specNum];j++)
        linEq->vector[j]=fs->ms_sum[j];
    }
  else if(abs(bgMode)==1)//constant background
    {
      linEq->dim=p->numFittedSimData[specNum]+1;
      
      //top-left entry
      linEq->matrix[0][0]=fs->sum1;
      
      //regular simulated data entires (bottom-right)
      for (j=0;j<p->numFittedSimData[specNum];j++)
        for (k=0;k<p->numFittedSimData[specNum];k++)
          linEq->matrix[1+j][1+k]=fs->ss_sum[j][k];
      
      //remaining entires
      for (j=0;j<p->numFittedSimData[specNum];j++)
        {     
          linEq->matrix[0][1+j]=fs->s_sum[j];
          linEq->matrix[1+j][0]=fs->s_sum[j];
        }
      
      linEq->vector[0]=fs->m_sum;
      for (j=0;j<p->numFittedSimData[specNum];j++)
        linEq->vector[j+1]=fs->ms_sum[j];
    }
  else if(abs(bgMode)==2)//linear background
    {
      linEq->dim=p->numFittedSimData[specNum]+2;
      
      //top-left 4 entries
      linEq->matrix[0][0]=fs->sum1;
      linEq->matrix[0][1]=fs->i_sum;
      linEq->matrix[1][0]=fs->i_sum;
      linEq->matrix[1][1]=fs->ii_sum;
      
      //regular simulated data entires (bottom-right)
      for (j=0;j<p->numFittedSimData[specNum];j++)
        for (k=0;k<p->numFittedSimData[specNum];k++)
          linEq->matrix[2+j][2+k]=fs->ss_sum[j][k];
      
      //remaining entires
      for (j=0;j<p->numFittedSimData[specNum];j++)
        {     
          linEq->matrix[0][2+j]=fs->s_sum[j];
          linEq->matrix[1][2+j]=fs->si_sum[j];
          linEq->matrix[2+j][0]=fs->s_sum[j];
          linEq->matrix[2+j][1]=fs->si_sum[j];
        }
      
      linEq->vector[0]=fs->m_sum;
      linEq->vector[1]=fs->mi_sum;
      for (j=0;j<p->numFittedSimData[specNum];j++)
        linEq->vector[j+2]=fs->ms_sum[j];
    }
  else if(abs(bgMode)==3)//quadratic background
    {
      linEq->dim=p->numFittedSimData[specNum]+3;
      
      //top-left 9 entries
      linEq->matrix[0][0]=fs->sum1;
      linEq->matrix[0][1]=fs->i_sum;
      linEq->matrix[0][2]=fs->ii_sum;
      linEq->matrix[1][0]=fs->i_sum;
      linEq->matrix[1][1]=fs->ii_sum;
      linEq->matrix[1][2]=fs->iii_sum;
      linEq->matrix[2][0]=fs->ii_sum;
      linEq->matrix[2][1]=fs->iii_sum;
      linEq->matrix[2][2]=fs->iiii_sum;
      
      //regular simulated data entires (bottom-right)
      for (j=0;j<p->numFittedSimData[specNum];j++)
        for (k=0;k<p->numFittedSimData[specNum];k++)
          linEq->matrix[3+j][3+k]=fs->ss_sum[j][k];
      
      //remaining entires
      for (j=0;j<p->numFittedSimData[specNum];j++)
        {     
          linEq->matrix[0][3+j]=fs->s_sum[j];
          linEq->matrix[1][3+j]=fs->si_sum[j];
          linEq->matrix[2][3+j]=fs->sii_sum[j];
          linEq->matrix[3+j][0]=fs->s_sum[j];
          linEq->matrix[3+j][1]=fs->si_sum[j];
          linEq->matrix[3+j][2]=fs->sii_sum[j];
        }
      
      linEq->vector[0]=fs->m_sum;
      linEq->vector[1]=fs->mi_sum;
      linEq->vector[2]=fs->mii_sum;
      for (j=0;j<p->numFittedSimData[specNum];j++)
        linEq->vector[j+3]=fs->ms_sum[j];
    }
  
  //solve system of equations and assign values
  if(!(solve_lin_eq(linEq)==1))
    {
      if(p->commonScaling==1)
        {
          if(fs->m_sum==0)
            printf("ERROR: Experiment data has no entries in the spectra to be fit.\n");
          else
            printf("ERROR: Could not determine background and scaling parameters.\n");
        }
      else
        {  
          if(fs->m_sum==0)
            printf("ERROR: Experiment data has no entries in the specified fitting region (spectrum %i, channel %i to %i).\n",p->spectrum[specNum],p->startCh[specNum],p->endCh[specNum]);
          else
            printf("ERROR: Could not determine background and scaling parameters for spectrum %i, channel %i to %i.\n",p->spectrum[specNum],p->startCh[specNum],p->endCh[specNum]);
        }
      exit(-1);
    }

}
