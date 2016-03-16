//function handles plotting of data, using the gnuplot_i library
void plotSpectra(pc_par * par, histdata * data, fitteddata * fdata )
{
  char str[256];
  int i,j,k;
  
  //allocate arrays to hold plot data
  double** x=allocateArrayD2(par->numSpectra,par->maxNumCh);
  double** yexp=allocateArrayD2(par->numSpectra,par->maxNumCh);
  double*** ysim=allocateArrayD3(NSIMDATA,par->numSpectra,par->maxNumCh);
  double** ybackground=allocateArrayD2(par->numSpectra,par->maxNumCh);
  double** ysimsum=allocateArrayD2(par->numSpectra,par->maxNumCh);
  
  for (i=0;i<par->numSpectra;i++)
    for (j=par->startCh[i];j<=par->endCh[i];j++)
      {
        x[i][j-par->startCh[i]]=(double)j;
        yexp[i][j-par->startCh[i]]=(double)data->expHist[par->spectrum[i]][j];
        if(par->addBackground>=1)
          ybackground[i][j-par->startCh[i]]=fdata->bgHist[i][j];
        ysimsum[i][j-par->startCh[i]]=ybackground[i][j-par->startCh[i]];
        for (k=0;k<par->numSimData;k++)
          {
            ysim[k][i][j-par->startCh[i]]=fdata->scaledSimHist[k][par->spectrum[i]][j];
            ysimsum[i][j-par->startCh[i]]+=fdata->scaledSimHist[k][par->spectrum[i]][j];
          }
      }

  plotOpen=1; 
  handle=gnuplot_init();
  printf("\nDATA PLOTS\n----------\nUse 'l' in the plotting window to switch between linear and logarithmic scale.\n");
  for(i=0;i<par->numSpectra;i++)
    {
      gnuplot_setstyle(handle,"steps");
      gnuplot_cmd(handle,"set xlabel 'Channel'");
      gnuplot_cmd(handle,"set ylabel 'Counts'");
      gnuplot_plot_xy(handle, x[i], yexp[i], par->endCh[i]-par->startCh[i]+1, "Experiment");
      if(par->plotOutput>1)//detailed plot
        {
          if(par->addBackground>=1)//plot background
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
          if((par->numSimData>1)||(par->addBackground>=1))//plot sum
            {
              gnuplot_setcolor(handle, "black");
              gnuplot_plot_xy(handle, x[i], ysimsum[i], par->endCh[i]-par->startCh[i]+1, "Simulation and Background(sum)");
              gnuplot_unsetcolor(handle);
            }
        }
      else //simple plot
        gnuplot_plot_xy(handle, x[i], ysimsum[i], par->endCh[i]-par->startCh[i]+1, "Simulation and Background(sum)");
      if(i!=par->numSpectra-1)//check whether we're showing the last plot
        printf("Showing plot for spectrum %i, channel %i to %i.  Press [ENTER] to continue...", par->spectrum[i],par->startCh[i],par->endCh[i]);
      else
        printf("Showing plot for spectrum %i, channel %i to %i.  Press [ENTER] to exit.", par->spectrum[i],par->startCh[i],par->endCh[i]);
      getc(stdin);
      gnuplot_resetplot(handle);
    }
  gnuplot_close(handle);
  plotOpen=0;
  free(x);
  free(yexp);
  free(ysim);
  free(ybackground);
  free(ysimsum);
  
}

//function run after CTRL-C, used to clean up temporary files generated
//by the plotting library
void sigint_cleanup()
{
  if(plotOpen==1)
    gnuplot_close(handle); //cleans up temporary files  
  exit(1); 
}
