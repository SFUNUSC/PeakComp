//function handles saving of fitted data
void saveSpectra(const par * p, const fitdata * fd)
{
  FILE *output;
  char str[256];
  int i,j,k;
  
  //allocate output arrays
  int ***outHist=allocateArrayI3(p->numSimData,p->numSpectra,S32K);
  int **outBGHist=allocateArrayI2(p->numSpectra,S32K); 
  
  //construct arrays
  for (i=0;i<p->numSimData;i++)
    for (j=0;j<p->numSpectra;j++)
      for (k=0;k<S32K;k++)
        outHist[i][j][k]=(int)(fd->scaledSimHist[i][j][k]);
  for (i=0;i<p->numSpectra;i++)
    for (j=0;j<S32K;j++)
      outBGHist[i][j]=(int)(fd->bgHist[i][j]);
        
  if(p->verbose>=0)
    printf("Saving scaled simulation data to output file(s)...\n");

  //save arrays to .mca files  
  if(p->addBackground>=1)
    {
      if((output=fopen("fit_background.mca","w"))==NULL)
        {
          printf("ERROR: Cannot open the output file fit_background.mca!\n");
          exit(-1);
        }
      for (i=0;i<p->numSpectra;i++)
        fwrite(outBGHist[i],S32K*sizeof(int),1,output);
      fclose(output);
    }
  for (i=0;i<p->numSimData;i++)
    {
      sprintf(str,"fit_sim%i.mca",i);
      if((output=fopen(str,"w"))==NULL)
        {
          printf("ERROR: Cannot open the output file %s!\n",str);
          exit(-1);
        }
      for (j=0;j<p->numSpectra;j++)
        fwrite(outHist[i][j],S32K*sizeof(int),1,output);
      fclose(output);
    }
  free(outHist);
  free(outBGHist);
  
}
