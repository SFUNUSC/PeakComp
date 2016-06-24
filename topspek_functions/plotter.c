//handles the gnuplot prompt
void plotPrompt(int cont)
{
	int c;
	char inp[256];
	if(cont==0)
 		printf("Enter 'g' for a gnuplot prompt or press [ENTER] to exit. ");
	else
		printf("Enter 'g' for a gnuplot prompt or press [ENTER] to continue. ");
	c=getc(stdin);
	if(c=='g')
		{
			printf("Enter 'exit' to return from the gnuplot prompt.\n");
			fgets(inp,256,stdin);
			while(strcmp(inp,"exit\n")!=0)
				{
					gnuplot_cmd(handle,inp);
					printf("gnuplot > ");
					fgets(inp,256,stdin);
				}
		}
	return;
}

//function handles plotting of data, using the gnuplot_i library
void plotSpectra(const par * p, const data * d, const fitdata * fd)
{
  char str[256];
  int i,j,k;
  
  //allocate arrays to hold plot data
  double** x=allocateArrayD2(p->numSpectra,p->maxNumCh);
  double** yexp=allocateArrayD2(p->numSpectra,p->maxNumCh);
  double*** ysim=allocateArrayD3(NSIMDATA,p->numSpectra,p->maxNumCh);
  double** ybackground=allocateArrayD2(p->numSpectra,p->maxNumCh);
  double** ysimsum=allocateArrayD2(p->numSpectra,p->maxNumCh);
  
  for (i=0;i<p->numSpectra;i++)
    for (j=p->startCh[i];j<=p->endCh[i];j++)
      {
        x[i][j-p->startCh[i]]=(double)j;
        yexp[i][j-p->startCh[i]]=(double)d->expHist[p->spectrum[i]][j];
        if(p->addBackground>=1)
          ybackground[i][j-p->startCh[i]]=fd->bgHist[i][j];
        ysimsum[i][j-p->startCh[i]]=ybackground[i][j-p->startCh[i]];
        for (k=0;k<p->numSimData;k++)
          {
            ysim[k][i][j-p->startCh[i]]=fd->scaledSimHist[k][i][j];
            ysimsum[i][j-p->startCh[i]]+=fd->scaledSimHist[k][i][j];
          }
      }

  plotOpen=1; 
  handle=gnuplot_init();
  printf("\nDATA PLOTS\n----------\nUse 'l' in the plotting window to switch between linear and logarithmic scale.\n");
  for(i=0;i<p->numSpectra;i++)
    {
      gnuplot_setstyle(handle,"steps");
      gnuplot_cmd(handle,"set xlabel 'Channel'");
      gnuplot_cmd(handle,"set ylabel 'Counts'");
      gnuplot_plot_xy(handle, x[i], yexp[i], p->endCh[i]-p->startCh[i]+1, "Experiment");
      if(p->plotOutput>1)//detailed plot
        {
          if(p->addBackground>=1)//plot background
            {
              gnuplot_setstyle(handle,"lines"); //plot background as a line
              gnuplot_plot_xy(handle, x[i], ybackground[i], p->endCh[i]-p->startCh[i]+1, "Background");
              gnuplot_setstyle(handle,"steps"); //set the plot style back
            }
          for (j=0;j<p->numSimData;j++)//plot individual sim data
            {
              sprintf(str,"Simulation (%s)",p->simDataName[j]);
              gnuplot_plot_xy(handle, x[i], ysim[j][i], p->endCh[i]-p->startCh[i]+1, str);
            }
          if((p->numSimData>1)||(p->addBackground>=1))//plot sum
            {
              gnuplot_setcolor(handle, "black");
              gnuplot_plot_xy(handle, x[i], ysimsum[i], p->endCh[i]-p->startCh[i]+1, "Simulation and Background(sum)");
              gnuplot_unsetcolor(handle);
            }
        }
      else //simple plot
        gnuplot_plot_xy(handle, x[i], ysimsum[i], p->endCh[i]-p->startCh[i]+1, "Simulation and Background(sum)");
      
      printf("Showing plot for spectrum %i, channel %i to %i.\n", p->spectrum[i],p->startCh[i],p->endCh[i]);
      if(i!=p->numSpectra-1)//check whether we're showing the last plot
      	plotPrompt(1);
      else
      	plotPrompt(0);
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
