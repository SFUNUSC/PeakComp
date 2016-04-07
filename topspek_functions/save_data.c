//function handles saving of fitted data
void saveSpectra(const pc_par * par, const fitteddata * fdata)
{
  FILE *output;
  char str[256];
  int i,j,k;
  
  int ***outHist=allocateArrayI3(par->numSimData,par->numSpectra,S32K); //allocate array
  
  if(par->verbose>=0)
    printf("Saving scaled simulation data to output file(s)...\n");
  
  //construct arrays
  for (i=0;i<par->numSpectra;i++)
    for (j=0;j<S32K;j++)
      {
        for (k=0;k<par->numSimData;k++)
          outHist[k][i][j]=(int)(fdata->scaledSimHist[k][par->spectrum[i]][j]);
      }

  //save arrays to .mca files  
  if(par->addBackground>=1)
    {
      if((output=fopen("fit_background.mca","w"))==NULL)
        {
          printf("ERROR: Cannot open the output file fit_background.mca!\n");
          exit(-1);
        }
      for (i=0;i<par->numSpectra;i++)
        fwrite(fdata->bgHist[i],S32K*sizeof(int),1,output);
      fclose(output);
    }
  for (i=0;i<par->numSimData;i++)
    {
      sprintf(str,"fit_sim%i.mca",i);
      if((output=fopen(str,"w"))==NULL)
        {
          printf("ERROR: Cannot open the output file %s!\n",str);
          exit(-1);
        }
      for (j=0;j<par->numSpectra;j++)
        fwrite(outHist[i][j],S32K*sizeof(int),1,output);
      fclose(output);
    }
  free(outHist);
  
}
