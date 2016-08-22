//reads and processes all data from input files to prepare it for fitting
void readAndProcessData(par * p, data * d)
{
	int i,j,k,l;
	
	//check that the number of spectra being compared is fine
  if(p->endSpectrum>=NSPECT)
    {
      printf("ERROR: A spectrum number specified in the parameter file is larger than the maximum value of %i.  Reduce it or increase NSPECT in peak_comp.h and recompile.\n",NSPECT);
      exit(-1);
    }
	
	//read in the data files
	memset(d->expHist,0,sizeof(d->expHist));
	memset(d->simHist,0,sizeof(d->simHist));
	readDataFile(p->expDataName,p->endSpectrum+1,d->expHist);
	for (i=0;i<p->numSimData;i++)
		readDataFile(p->simDataName[i],p->endSpectrum+1,d->simHist[i]);

	//Data read in from files has spectra indexed as in the original MCA file.
	//But the user can specify any range or order of spectra to analyze.
	//Need to convert to files with spectra indexed in the order of spectra to
	//be analyzed, starting with 0 (the first spectrum to be analyzed).
	
	//generate experiment data for fitting (0 indexed in spectra)
	memset(d->fittedExpHist,0,sizeof(d->fittedExpHist));
	for (i=0;i<p->numSpectra;i++)
  	for (j=0;j<=p->endSpectrum;j++)
      if(p->spectrum[i]==j)
        {
        	//printf("Spectrum: %i, index: %i\n",p->spectrum[i],i);
          if(p->fixBG[i]!=0)
          	{
		          for (k=0;k<S32K;k++)
		            d->fittedExpHist[i][k]=d->expHist[j][k] - p->fixedBGPar[i][0] - p->fixedBGPar[i][1]*k - p->fixedBGPar[i][2]*k*k;
            }
          else
          	{
		        	for (k=0;k<S32K;k++)
		        		d->fittedExpHist[i][k]=d->expHist[j][k];
          	}
        }
	
	//generate simulated data for fitting (0 indexed in spectra)
	int fi[NSPECT];//index for simulated data to be fitted
  memset(fi,0,sizeof(fi));
  memset(d->fittedSimHist,0,sizeof(d->fittedSimHist));
  for (i=0;i<p->numSimData;i++)
    {
    	for (j=0;j<p->numSpectra;j++)
  			for (k=0;k<=p->endSpectrum;k++)
  				if(p->spectrum[j]==k)
						{
							//determine whether simulated data is fitted and read into histograms for fitting as needed
							if(p->simDataFixedAmp[i][j]==0)
								{
								  for (l=0;l<S32K;l++)
								    d->fittedSimHist[fi[j]][j][l]=d->simHist[i][k][l];
								  fi[j]++;
								}
							else if(p->simDataFixedAmp[i][j]==2)//data scaling is fixed relative to the previous fit dataset
								{
								  if(fi[j]>0)
								    for (l=0;l<S32K;l++)
								      d->fittedSimHist[fi[j]-1][j][l]+=p->simDataFixedAmpValue[i][j]*d->simHist[i][k][l];//add this data to the data that it is scaled relative to
								}
							else if(p->simDataFixedAmp[i][j]==1)//data scaling is fixed to a specified value (will not be fitted)
								{
									for (l=0;l<S32K;l++)
										d->fittedExpHist[j][l]-=p->simDataFixedAmpValue[i][j]*d->simHist[i][k][l];
								}
						}
    }
  
  /*//debug  
  for(j=0;j<p->numSpectra;j++)
  	printf("fi[%i]: %i, numFittedSimData[%i]: %i\n",j,fi[j],j,p->numFittedSimData[j]);*/
  	
  for(j=0;j<p->numSpectra;j++)
		if(fi[j]!=p->numFittedSimData[j])
			{
				printf("ERROR: not all fitted data was imported correctly!\n");
				exit(-1);
			}
    
  if(p->verbose>=0) printf("Spectra read in...\n");
  
  if(p->peakSearch==1)
    findFittingWindow(p,d);//find peaks and shift fitting windows (see peak_window.c)
}
