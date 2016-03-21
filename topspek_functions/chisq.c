//function compares spectra and gets chisq and other stats
void compareSpectra(pc_par * par, histdata * data, fitteddata * fdata)
{
  //initialize values
  int numFittedParameters[NSPECT],sumFittedParameters;
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
  
  sumFittedParameters=0;
  for (i=0;i<par->numSpectra;i++)
    {
      numFittedParameters[i] = par->numFittedSimData;
      if(par->fixBG[i]==0)
        if(par->addBackground>=1)
          numFittedParameters[i] += par->addBackground + 1;
      sumFittedParameters+=numFittedParameters[i];
    }
  
  //compute chisq for data in the spectra
  for (i=0;i<par->numSpectra;i++)
    for (j=par->startCh[i];j<=par->endCh[i];j++)
      if(data->expHist[par->spectrum[i]][j]!=0)//avoid dividing by zero
        {
          //get the sum of all simulated data in the given bin 
          sumSimValue = fdata->bgHist[i][j];
          for (k=0;k<par->numSimData;k++)
            sumSimValue+=fdata->scaledSimHist[k][par->spectrum[i]][j];
          //increment the chisq value
          spectChisq[i]+=((data->expHist[par->spectrum[i]][j]-sumSimValue)*(data->expHist[par->spectrum[i]][j]-sumSimValue))/((double)data->expHist[par->spectrum[i]][j]);
        }
      else
        binsSkipped[i]++;  
  for (i=0;i<par->numSpectra;i++)
    sumBinsSkipped+=binsSkipped[i];
    
  //compute total chisq and reduced total chisq
  for (i=0;i<par->numSpectra;i++)
    {
      numBinsUsed[i]=(par->endCh[i]-par->startCh[i]+1)-binsSkipped[i];
      sumBinsUsed+=numBinsUsed[i];
      chisq+=spectChisq[i];
      spectRedChisq[i]=spectChisq[i]/(numBinsUsed[i]-numFittedParameters[i]-1);
    }
  redChisq=chisq/(sumBinsUsed-sumFittedParameters-1);
  
  //print output
  printf("COMPARISON DATA\n---------------\n");
  if(sumBinsSkipped>0)
    printf("Warning: some of the bins in the experiment data have values of zero.  These have been skipped when calculating chisq.  Bins skipped: %i\n\n",sumBinsSkipped);
  if(par->numSpectra>1)
    for (i=0;i<par->numSpectra;i++)
      printf("Spectrum %i, channel %i to %i - chisq: %f, reduced chisq: %f\n",par->spectrum[i],par->startCh[i],par->endCh[i],spectChisq[i],spectRedChisq[i]);
  printf("Chisq (total): %f\n",chisq);
  printf("Number of bins (total): %i\n",sumBinsUsed);
  printf("Number of fitted parameters (total): %i\n",sumFittedParameters);
  printf("Reduced chisq (total): %f\n",redChisq);
  
}
