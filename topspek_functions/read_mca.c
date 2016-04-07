//function reads an .mca file into an integer array and returns the array
void readMCA(const char * filename, const int numSpec, int outHist[NSPECT][S32K])
{
  FILE *inp;
  int i;
  
  if((inp=fopen(filename,"r"))==NULL)
    {
      printf("ERROR: Cannot open the .mca file: %s\n",filename);
      exit(-1);
    }

  for (i=0;i<numSpec;i++)
    if(fread(outHist[i],S32K*sizeof(int),1,inp)!=1)
      {
        printf("ERROR: Error reading spectrum %i from the .mca file: %s\n",i,filename);
        printf("Verify that the format and number of spectra in the file are correct.\n");
        exit(-1);
      }

  fclose(inp);
  
}
