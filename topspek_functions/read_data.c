//function reads an .mca file into an integer array and returns the array
void readMCA(FILE * inp, const char * filename, const int numSpec, int outHist[NSPECT][S32K])
{
  int i;

  for (i=0;i<numSpec;i++)
    if(fread(outHist[i],S32K*sizeof(int),1,inp)!=1)
      {
        printf("ERROR: Cannot read spectrum %i from the .mca file: %s\n",i,filename);
        printf("Verify that the format and number of spectra in the file are correct.\n");
        exit(-1);
      }
  
}

//function reads an .spe file into an integer array and returns the array
void readSPE(FILE * inp, const char * filename, const int numSpec, int outHist[NSPECT][S32K])
{
	int i;
  char header[36];
  float inpHist[S32K];
  memset(outHist,0,sizeof(*outHist));

	if(fread(header,36,1,inp)!=1)
    {
      printf("ERROR: Cannot read header from the .spe file: %s\n",filename);
      printf("Verify that the format of the file is correct.\n");
      exit(-1);
    }
  if(fread(inpHist,4096*sizeof(float),1,inp)!=1)
    {
      printf("ERROR: Cannot read spectrum from the .spe file: %s\n",filename);
      printf("Verify that the format of the file is correct.\n");
      exit(-1);
    }
  
  //convert input data to integer
  for(i=0;i<S32K;i++)
  	outHist[0][i]=(int)inpHist[i];
  
  //copy data to each spectrum
  for(i=1;i<numSpec;i++)
  	memcpy(outHist[i],outHist[0],S32K*sizeof(int));
  
}

void readDataFile(const char * filename, const int numSpec, int outHist[NSPECT][S32K])
{
	FILE *inp;
	if((inp=fopen(filename,"r"))==NULL)
    {
      printf("ERROR: Cannot open the input file: %s\n",filename);
      printf("Check that the file exists.\n");
      exit(-1);
    }
	
	const char *dot = strrchr(filename, '.');//get the file extension
	if(strcmp(dot + 1,"mca")==0)
		readMCA(inp, filename, numSpec, outHist);
	else if(strcmp(dot + 1,"spe")==0)
		readSPE(inp, filename, numSpec, outHist);
	else
		{
			printf("ERROR: Improper type of input file: %s\n",filename);
      printf("Integer array (.mca) and radware (.spe) files are supported.\n");
      exit(-1);
		}
	
	fclose(inp);
}
