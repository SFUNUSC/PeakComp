void findFittingWindow(pc_par * par, const histdata * data)
{
  int i;
  if(par->verbose>=0)
    printf("\nPEAK FINDER\n-----------\n");
  peak_search_par pspar;
  for (i=0;i<par->numSpectra;i++)
    {
      pspar.searchMin=par->startCh[i];
      pspar.searchMax=par->endCh[i];
      pspar.windowSize=(int)(pspar.searchMax-pspar.searchMin)/10;
      peak_fit_par pfpar = findPeak(&pspar, data->expHist[par->spectrum[i]], S32K);
      if(par->verbose>=0)
        printf("Spectrum %i: peak found with centroid at channel %i.\n",i,(int)pfpar.centroid);
      int range=abs(pspar.searchMax-pspar.searchMin);
      par->startCh[i]=pfpar.centroid-(int)(range/2);
      par->endCh[i]=pfpar.centroid+(int)(range/2);
    }
  if(par->verbose>=0)
    printf("Fit window(s) centered on peak(s)...\n");

}
