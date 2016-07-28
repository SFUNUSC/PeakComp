//reads and processes all data from input files to prepare it for fitting
void readAndProcessData(par * p, data * d)
{
	int i,j,k;
	
	//check that the number of spectra being compared is fine
  if(p->endSpectrum>=NSPECT)
    {
      printf("ERROR: A spectrum number specified in the parameter file is larger than the maximum value of %i.  Reduce it or increase NSPECT in peak_comp.h and recompile.\n",NSPECT);
      exit(-1);
    }
	
  //read in the data files
  readDataFile(p->expDataName,p->endSpectrum+1,d->expHist);
  readDataFile(p->expDataName,p->endSpectrum+1,d->fittedExpHist);
  for (i=0;i<p->numSimData;i++) 
    {
      readDataFile(p->simDataName[i],p->endSpectrum+1,d->simHist[i]);
      readDataFile(p->simDataName[i],p->endSpectrum+1,d->fittedSimHist[i]);
    }
    
  //generate data for fitting
  for (j=0;j<=p->endSpectrum;j++)
    for (i=0;i<p->numSpectra;i++)
      if(p->spectrum[i]==j)
        {
          if(p->fixBG[i]!=0)
            for (k=0;k<S32K;k++)
              d->fittedExpHist[j][k]=d->expHist[j][k] - p->fixedBGPar[i][0] - p->fixedBGPar[i][1]*k - p->fixedBGPar[i][2]*k*k;
        }
    
  int fi=0;//index for simulated data to be fitted
  for (i=0;i<p->numSimData;i++)
    {
          
      //determine whether simulated data is fitted and read into histograms for fitting as needed
      if(p->simDataFixedAmp[i][0]==0)
        {
          for (j=0;j<=p->endSpectrum;j++)
            for (k=0;k<S32K;k++)
              d->fittedSimHist[fi][j][k]=d->simHist[i][j][k];
          fi++;
        }
      else if(p->simDataFixedAmp[i][0]==2)//data scaling is fixed relative to the previous fitted dataset
        {
          if(fi>0)
            for (j=0;j<=p->endSpectrum;j++)
              for (k=0;k<S32K;k++)
                d->fittedSimHist[fi-1][j][k]+=p->simDataFixedAmpValue[i][j]*d->simHist[i][j][k];//add this data to the data that it is scaled relative to
        }
      else if(p->simDataFixedAmp[i][0]==1)//data scaling is fixed to a specified value (will not be fitted)
        for (j=0;j<=p->endSpectrum;j++)
          for (k=0;k<S32K;k++)
            d->fittedExpHist[j][k]-=p->simDataFixedAmpValue[i][j]*d->simHist[i][j][k];
    }
	if(fi!=p->numFittedSimData)
		{
			printf("ERROR: not all fitted data was imported correctly!\n");
			exit(-1);
		}
    
  if(p->verbose>=0) printf("Spectra read in...\n");
  
  if(p->peakSearch==1)
    findFittingWindow(p,d);//find peaks and shift fitting windows (see peak_window.c)
}
