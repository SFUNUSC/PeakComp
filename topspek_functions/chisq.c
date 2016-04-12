//function compares spectra and gets chisq and other stats
void compareSpectra(const par * p, const data * d, const fitdata * fd)
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
  for (i=0;i<p->numSpectra;i++)
    {
      numFittedParameters[i] = p->numFittedSimData;
      if(p->fixBG[i]==0)
        if(p->addBackground>=1)
          numFittedParameters[i] += p->addBackground + 1;
      sumFittedParameters+=numFittedParameters[i];
    }
  
  //compute chisq for data in the spectra
  for (i=0;i<p->numSpectra;i++)
    for (j=p->startCh[i];j<=p->endCh[i];j++)
      if(d->expHist[p->spectrum[i]][j]!=0)//avoid dividing by zero
        {
          //get the sum of all simulated data in the given bin 
          sumSimValue = fd->bgHist[i][j];
          for (k=0;k<p->numSimData;k++)
            sumSimValue+=fd->scaledSimHist[k][i][j];
          //increment the chisq value
          spectChisq[i]+=((d->expHist[p->spectrum[i]][j]-sumSimValue)*(d->expHist[p->spectrum[i]][j]-sumSimValue))/((double)d->expHist[p->spectrum[i]][j]);
        }
      else
        binsSkipped[i]++;  
  for (i=0;i<p->numSpectra;i++)
    sumBinsSkipped+=binsSkipped[i];
    
  //compute total chisq and reduced total chisq
  for (i=0;i<p->numSpectra;i++)
    {
      numBinsUsed[i]=(p->endCh[i]-p->startCh[i]+1)-binsSkipped[i];
      sumBinsUsed+=numBinsUsed[i];
      chisq+=spectChisq[i];
      spectRedChisq[i]=spectChisq[i]/(numBinsUsed[i]-numFittedParameters[i]-1);
    }
  redChisq=chisq/(sumBinsUsed-sumFittedParameters-1);
  
  //print output
  if(p->verbose>=0)
    {
      printf("COMPARISON DATA\n---------------\n");
      if(sumBinsSkipped>0)
        printf("Warning: some of the bins in the experiment data have values of zero.  These have been skipped when calculating chisq.  Bins skipped: %i\n\n",sumBinsSkipped);
      if(p->numSpectra>1)
        for (i=0;i<p->numSpectra;i++)
          printf("Spectrum %i, channel %i to %i - chisq: %f, reduced chisq: %f\n",p->spectrum[i],p->startCh[i],p->endCh[i],spectChisq[i],spectRedChisq[i]);
      printf("Chisq (total): %f\n",chisq);
      printf("Number of bins (total): %i\n",sumBinsUsed);
      printf("Number of fitted parameters (total): %i\n",sumFittedParameters);
      printf("Reduced chisq (total): %f\n",redChisq);
    }
  else if(p->verbose==-1)
    printf("%f\n",redChisq);//only print the reduced chisq
  
}
