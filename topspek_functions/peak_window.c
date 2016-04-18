void findFittingWindow(par * p, const data * d)
{
  int i;
  if(p->verbose>=0)
    printf("\nPEAK FINDER\n-----------\n");
  peak_search_par pspar;
  for (i=0;i<p->numSpectra;i++)
    {
      pspar.searchMin=p->startCh[i];
      pspar.searchMax=p->endCh[i];
      pspar.windowSize=(int)(pspar.searchMax-pspar.searchMin)/10;
      peak_fit_par pfpar = findPeak(&pspar, d->expHist[p->spectrum[i]], S32K);
      if(p->verbose>=0)
        printf("Spectrum %i: peak found with centroid at channel %i.\n",i,(int)pfpar.centroid);
      int range=abs(pspar.searchMax-pspar.searchMin);
      if(p->peakSearchWidth>0)
        range=p->peakSearchWidth;
      p->startCh[i]=pfpar.centroid-(int)(range/2);
      p->endCh[i]=pfpar.centroid+(int)(range/2);
    }
  if(p->verbose>=0)
    printf("Fit window(s) centered on peak(s)...\n");

}
