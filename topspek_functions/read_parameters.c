//function reads parameter files for the topspek code
void readParFile(const char * fileName, par * p) 
{
  FILE *config;
  char str[256],str1[256],str2[256];
  int index=0;
  p->numSpectra=0;
  p->endSpectrum=0;
  p->maxNumCh=0;
  p->numSimData=0;
  p->numFittedSimData=0;
  memset(p->fixBG,0,sizeof(p->fixBG));
  if((config=fopen(fileName,"r"))==NULL)
    {
      printf("ERROR: Cannot open the config file %s!\n",fileName);
      exit(-1);
    }
  while(!(feof(config)))//go until the end of file is reached
    {
      if(fgets(str,256,config)!=NULL)
        {
        
          if(p->numSimData<NSIMDATA)
            if(sscanf(str,"%i %i %i",&p->spectrum[index],&p->startCh[index],&p->endCh[index])!=3) //no spectrum and channel data
              if(sscanf(str,"%s %s %lf",p->simDataName[p->numSimData],str1,&p->simDataFixedAmpValue[p->numSimData])==3) //simulated dataset info
                {
                  if(strcmp(str1,"yes")==0)
                    p->simDataFixedAmp[p->numSimData]=1;
                  else if(strcmp(str1,"rel")==0)
                    p->simDataFixedAmp[p->numSimData]=2;
                  else
                    {
                      p->simDataFixedAmp[p->numSimData]=0;
                      strcpy(p->fittedSimDataName[p->numFittedSimData],p->simDataName[p->numSimData]);
                      p->numFittedSimData++;
                    }
                  p->numSimData++;
                }
              
          if(index<NSPECT)
            {
              if(sscanf(str,"%i %i %i %s",&p->spectrum[index],&p->startCh[index],&p->endCh[index],str1)==3) //spectrum and channel data
                {
                  if(p->spectrum[index]>p->endSpectrum)
                    p->endSpectrum=p->spectrum[index];
                  if((p->endCh[index]-p->startCh[index]+1)>p->maxNumCh)
                    p->maxNumCh=p->endCh[index]-p->startCh[index]+1;
                  index++;
                  p->numSpectra++;
                }
              if(sscanf(str,"%i %i %i %s %lf %lf %lf",&p->spectrum[index],&p->startCh[index],&p->endCh[index],str1,&p->fixedBGPar[index][0],&p->fixedBGPar[index][1],&p->fixedBGPar[index][2])>=5) //spectrum, channel, and background data
                {
                  if(p->spectrum[index]>p->endSpectrum)
                    p->endSpectrum=p->spectrum[index];
                  if((p->endCh[index]-p->startCh[index]+1)>p->maxNumCh)
                    p->maxNumCh=p->endCh[index]-p->startCh[index]+1;
                  if(strcmp(str1,"yes")==0)
                    p->fixBG[index]=1;
                  index++;
                  p->numSpectra++;
                }
            }
              
          if(sscanf(str,"%s %s",str1,str2)==2) //single parameter data
            {
              if(strcmp(str1,"EXPERIMENT_DATA")==0)
                strcpy(p->expDataName,str2);
              if(strcmp(str1,"ADD_BACKGROUND")==0)
                {
                  if(strcmp(str2,"quad")==0)
                    p->addBackground=3;
                  else if((strcmp(str2,"yes")==0)||(strcmp(str2,"lin")==0))
                    p->addBackground=2;
                  else if(strcmp(str2,"const")==0)
                    p->addBackground=1;
                  else
                    p->addBackground=0;
                }
              if(strcmp(str1,"PEAK_SEARCH")==0)
                {
                  if(strcmp(str2,"yes")==0)
                    p->peakSearch=1;
                  else
                    p->peakSearch=0;
                }
              if(strcmp(str1,"PEAK_SEARCH_SET_WINDOW")==0)
                {
                  sscanf(str2,"%i",&p->peakSearchWidth);
                }
              if(strcmp(str1,"COMMON_SCALING")==0)
                {
                  if(strcmp(str2,"yes")==0)
                    p->commonScaling=1;
                  else
                    p->commonScaling=0;
                }
              if(strcmp(str1,"PLOT_OUTPUT")==0)
                {
                  if(strcmp(str2,"yes")==0)
                    p->plotOutput=1;
                  else if(strcmp(str2,"detailed")==0)
                    p->plotOutput=2;
                  else
                    p->plotOutput=0;
                }
              if(strcmp(str1,"SAVE_OUTPUT")==0)
                {
                  if(strcmp(str2,"yes")==0)
                    p->saveOutput=1;
                  else
                    p->saveOutput=0;
                }
              if(strcmp(str1,"VERBOSITY")==0)
                {
                  if(strcmp(str2,"chisq")==0)
                    p->verbose=-1;
                  else
                    p->verbose=0;
                }
            }
          
          if(sscanf(str,"%s %s",str1,str2)==1) //listing of simulated data
            {
              if(strcmp(str1,"<---END_OF_PARAMETERS--->")==0)
                break;
              else if(strcmp(str1,"SIMULATED_DATA")!=0)
                if(p->numSimData<NSIMDATA)
                  {
                    strcpy(p->simDataName[p->numSimData],str1);
                    p->simDataFixedAmp[p->numSimData]=0;
                    p->simDataFixedAmpValue[p->numSimData]=1.;
                    p->numFittedSimData++;
                    p->numSimData++;
                  }
            }
        }
    }
  fclose(config);
  
  //correct parameters
  for(index=0;index<p->numSpectra;index++)
    if(p->fixBG[index]==0)
      p->fitAddBackground[index]=p->addBackground;
    else if(p->fixBG[index]==1)
      p->fitAddBackground[index]=0;
  
  if(p->addBackground==1)
    for(index=0;index<p->numSpectra;index++)
      {
        p->fixedBGPar[index][1]=0;
        p->fixedBGPar[index][2]=0;
      }
  if(p->addBackground==2)
    for(index=0;index<p->numSpectra;index++)
      p->fixedBGPar[index][2]=0;
  
  //print parameters read from the file
  if(p->verbose>=0)
    {
      
      if(strcmp(p->expDataName,"")==0)
        {
          printf("ERROR: No experiment data file specified in the parameter file!\n");
          exit(-1);
        }
      else  
        printf("\nTaking experiment data from file: %s\n",p->expDataName);
      for(index=0;index<p->numSimData;index++)
        {
          printf("Taking simulated data from file (%i of %i): %s\n",index+1,p->numSimData,p->simDataName[index]);
          if(p->simDataFixedAmp[index]==1)
            printf("Fixing scaling factor for this data to %lf\n",p->simDataFixedAmpValue[index]);
          if(p->simDataFixedAmp[index]==2)
            printf("Fixing scaling factor for this data to a factor of %lf relative to the last fitted data.\n",p->simDataFixedAmpValue[index]);
        }
      if(p->peakSearch==0)
        for(index=0;index<p->numSpectra;index++)
          printf("Will compare spectrum %i from channels %i to %i.\n",p->spectrum[index],p->startCh[index],p->endCh[index]);
      else
        {
          for(index=0;index<p->numSpectra;index++)
            printf("Will search for a peak in spectrum %i from channels %i to %i.\n",p->spectrum[index],p->startCh[index],p->endCh[index]);
          if(p->peakSearchWidth>0)
            printf("Will set fitting window width to %i channels around each peak found.\n",p->peakSearchWidth);
        }
      if(p->commonScaling==1)
        printf("Will use common scaling and background for all spectra in each data file.\n");
      if(p->addBackground==0)
        printf("Will not add background to simulated data.\n");
      if(p->addBackground==1)
        printf("Will add a constant background to simulated data.\n");
      if(p->addBackground==2)
        printf("Will add a linear background to simulated data.\n");
      if(p->addBackground==3)
        printf("Will add a quadratic background to simulated data.\n");
      if(p->addBackground==1)
        for(index=0;index<p->numSpectra;index++)
          if(p->fixBG[index]==1)
            printf("Fixing background amplitude to %lf for spectrum %i, channels %i to %i.\n",p->fixedBGPar[index][0],p->spectrum[index],p->startCh[index],p->endCh[index]);
      if(p->addBackground==2)
        for(index=0;index<p->numSpectra;index++)
          if(p->fixBG[index]==1)
            printf("Fixing background parameters to A = %lf, B = %lf for spectrum %i, channels %i to %i.\n",p->fixedBGPar[index][0],p->fixedBGPar[index][1],p->spectrum[index],p->startCh[index],p->endCh[index]);
      if(p->addBackground==3)
        for(index=0;index<p->numSpectra;index++)
          if(p->fixBG[index]==1)
            printf("Fixing background parameters to A = %lf, B = %lf, C = %lf for spectrum %i, channels %i to %i.\n",p->fixedBGPar[index][0],p->fixedBGPar[index][1],p->fixedBGPar[index][2],p->spectrum[index],p->startCh[index],p->endCh[index]);
      if(p->plotOutput==0)
        printf("Will not plot output data.\n");
      if(p->plotOutput==1)
        printf("Will plot output data.\n");
      if(p->plotOutput==2)
        printf("Will plot detailed output data.\n");
      if(p->saveOutput==0)
        printf("Will not save fitted simulation data.\n");
      if(p->saveOutput==1)
        {
          if(p->addBackground==0)
            {
              if(p->numSimData>1)
                printf("Will save fitted simulation data to files fit_sim0.mca through fit_sim%i.mca.\n",p->numSimData-1);
              else
                printf("Will save fitted simulation data to file fit_sim0.mca.\n");
            }
          else if(p->addBackground==1)
            {
              if(p->numSimData>1)
                printf("Will save fitted simulation data to files fit_background.mca and fit_sim0.mca through fit_sim%i.mca.\n",p->numSimData-1);
              else
                printf("Will save fitted simulation data to files fit_background.mca and fit_sim0.mca.\n");
            }
        }
      
      printf("Finished reading parameter file...\n");
    }
  
}
