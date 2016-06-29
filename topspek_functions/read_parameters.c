#include "read_parameters.h"

//parses a string into an option
par_option * readOption(char * str)
{
	par_option *opt=(par_option*)calloc(1,sizeof(par_option));
	char *tok;
	tok=strtok (str,"(,)");
	strcpy(opt->name,tok);
	int i=0;
	while (tok != NULL)
	{
		tok = strtok (NULL, "(,)");
		if(tok!=NULL)
			{
				strcpy(opt->par[i],tok);
				i++;
				if(i>=MAX_OPT_PAR)
					{
						printf("WARNING: too many parameters specified in parameter file option: %s\n",opt->name);
						break;
					}
			}
	}
	opt->numPar=i-1;
	return (par_option*)opt;
}

//function reads parameter files for the topspek code
void readParFile(const char * fileName, par * p) 
{
  FILE *config;
  char str[256],str1[256];
  int i,j;
  par_option *opt[MAX_NUM_OPT];
  for(i=0;i<MAX_NUM_OPT;i++)
  	opt[i]=(par_option*)malloc(sizeof(par_option));
  int numOpt=0;
  
  //default values
  p->numSpectra=0;
  p->endSpectrum=0;
  p->maxNumCh=0;
  p->numSimData=0;
  p->numFittedSimData=0;
  p->channelScaling=1.;
  memset(p->fixBG,0,sizeof(p->fixBG));
  
  //open the file and read all parameters
  if((config=fopen(fileName,"r"))==NULL)
    {
      printf("ERROR: Cannot open the parameter file %s!\n",fileName);
      exit(-1);
    }
  while(!(feof(config)))//go until the end of file is reached
    {
			if(fgets(str,256,config)!=NULL)
				{
					sscanf(str,"%s",str1);
					if(strcmp(str1,"<---END_OF_PARAMETERS--->")==0)
          	break;
        	opt[numOpt]=readOption(str);
        	numOpt++;
        	if(numOpt>=MAX_NUM_OPT)
        		{
        			printf("ERROR: The maximum number of options (%i) in the parameter file %s has been exceeded!\nIncrease the value of MAX_NUM_OPT in read_parameters.h and recompile.\n",MAX_NUM_OPT,fileName);
							exit(-1);
        		}
				}
		}
	fclose(config);
	
	//parse the parameters which were just read in
	
	//preprocess and deal with the spectra first (other parameters depend on this) 
	for(i=0;i<numOpt;i++)
		{
			if(strcmp(opt[i]->name,"SP")==0)
				{
					if(opt[i]->numPar<3)
						{
							printf("ERROR: Not enough information was specified for one or more spectra.\n");
							exit(-1);
						}
					p->spectrum[p->numSpectra]=atoi(opt[i]->par[0]);//spectrum number
					p->startCh[p->numSpectra]=atoi(opt[i]->par[1]);//start channel
					p->endCh[p->numSpectra]=atoi(opt[i]->par[2]);//end channel
					if(strcmp(opt[i]->par[3],"yes")==0)
          	p->fixBG[i]=1;//fix background
          p->fixedBGPar[i][0]=atof(opt[i]->par[4]);
          p->fixedBGPar[i][1]=atof(opt[i]->par[5]);
          p->fixedBGPar[i][2]=atof(opt[i]->par[6]);
					if(p->spectrum[p->numSpectra]>p->endSpectrum)
						p->endSpectrum=p->spectrum[p->numSpectra];
					if((p->endCh[p->numSpectra]-p->startCh[p->numSpectra]+1)>p->maxNumCh)
						p->maxNumCh=p->endCh[p->numSpectra]-p->startCh[p->numSpectra]+1;
					p->numSpectra++;
				}
		}
	
	//deal with all other parameters
	for(i=0;i<numOpt;i++)
		{
			if(strcmp(opt[i]->name,"EXPERIMENT_DATA")==0)
				{
					strcpy(p->expDataName,opt[i]->par[0]);
				}
			else if(strcmp(opt[i]->name,"ADD_BACKGROUND")==0)
        {
          if(strcmp(opt[i]->par[0],"quad")==0)
            p->addBackground=3;
          else if((strcmp(opt[i]->par[0],"yes")==0)||(strcmp(opt[i]->par[0],"lin")==0))
            p->addBackground=2;
          else if(strcmp(opt[i]->par[0],"const")==0)
            p->addBackground=1;
          else
            p->addBackground=0;
        }
      else if(strcmp(opt[i]->name,"PEAK_SEARCH")==0)
        {
          if(strcmp(opt[i]->par[0],"yes")==0)
            p->peakSearch=1;
          else
            p->peakSearch=0;
        }
      else if(strcmp(opt[i]->name,"PEAK_SEARCH_SET_WINDOW")==0)
        {
          p->peakSearchWidth=atoi(opt[i]->par[0]);
        }
      else if(strcmp(opt[i]->name,"COMMON_SCALING")==0)
        {
          if(strcmp(opt[i]->par[0],"yes")==0)
            p->commonScaling=1;
          else
            p->commonScaling=0;
        }
      else if(strcmp(opt[i]->name,"PLOT_OUTPUT")==0)
        {
          if(strcmp(opt[i]->par[0],"yes")==0)
            p->plotOutput=1;
          else if(strcmp(opt[i]->par[0],"detailed")==0)
            p->plotOutput=2;
          else
            p->plotOutput=0;
        }
      else if(strcmp(opt[i]->name,"SAVE_OUTPUT")==0)
        {
          if(strcmp(opt[i]->par[0],"yes")==0)
            p->saveOutput=1;
          else
            p->saveOutput=0;
        }
      else if(strcmp(opt[i]->name,"VERBOSITY")==0)
        {
          if(strcmp(opt[i]->par[0],"chisq")==0)
            p->verbose=-1;
          else
            p->verbose=0;
        }
      else if(strcmp(opt[i]->name,"CHANNEL_SCALING")==0)
        {
        	p->channelScaling=atof(opt[i]->par[0]);
        	if(p->channelScaling<=0.)
        		p->channelScaling=1.;
        }
      else if(strcmp(opt[i]->name,"FORCE_POSITIVE_SCALING")==0)
        {
        	if(strcmp(opt[i]->par[0],"yes")==0)
            p->forcePositiveS=1;
          else
            p->forcePositiveS=0;
        }
      else if(strcmp(opt[i]->name,"INDEPENDENT_SP")==0)
        {
        	if(strcmp(opt[i]->par[0],"yes")==0)
            p->indSpectra=1;
          else
            p->indSpectra=0;
        }
      else if(strcmp(opt[i]->name,"DATA")==0)
        {
        	p->simDataCommonScaling[p->numSimData]=1;
        	strcpy(p->simDataName[p->numSimData],opt[i]->par[0]);//filename
					
					//set up scaling
    			if((strcmp(opt[i]->par[1],"rel_scaling")==0)||(strcmp(opt[i]->par[1],"rel")==0))
    				p->simDataFixedAmp[p->numSimData]=2;
    			else if((strcmp(opt[i]->par[1],"abs_scaling")==0)||(strcmp(opt[i]->par[1],"abs")==0))
    				p->simDataFixedAmp[p->numSimData]=1;
    			else
    				{
							p->simDataFixedAmp[p->numSimData]=0;
							strcpy(p->fittedSimDataName[p->numFittedSimData],p->simDataName[p->numSimData]);
							p->numFittedSimData++;
    				}		
        	if(opt[i]->numPar==3)//scaling is the same for each spectrum
        		{
        			for(j=0;j<p->numSpectra;j++)
        				p->simDataFixedAmpValue[p->numSimData][j]=atof(opt[i]->par[2]);
        		}
        	else if(opt[i]->numPar>3)//scaling is different for each spectrum
        		{
        			for(j=0;j<p->numSpectra;j++)
        				if(j+2<MAX_OPT_PAR)
        					p->simDataFixedAmpValue[p->numSimData][j]=atof(opt[i]->par[j+2]);
        			p->simDataCommonScaling[p->numSimData]=0;
        		}
        	p->numSimData++;
        }
		}
  
  //correct parameters
  if(p->channelScaling!=1.)
  	for(i=0;i<p->numSpectra;i++)
  		{
  			//rescale channel ranges
  			p->startCh[i]=(int)(p->startCh[i]*p->channelScaling);
  			p->endCh[i]=(int)(p->endCh[i]*p->channelScaling);
  			if((p->endCh[i]-p->startCh[i]+1)>p->maxNumCh)
        	p->maxNumCh=p->endCh[i]-p->startCh[i]+1;
  		}
  for(i=0;i<p->numSpectra;i++)
    if(p->fixBG[i]==0)
      p->fitAddBackground[i]=p->addBackground;
    else if(p->fixBG[i]==1)
      p->fitAddBackground[i]=0;
  
  if(p->addBackground==1)
    for(i=0;i<p->numSpectra;i++)
      {
        p->fixedBGPar[i][1]=0;
        p->fixedBGPar[i][2]=0;
      }
  if(p->addBackground==2)
    for(i=0;i<p->numSpectra;i++)
      p->fixedBGPar[i][2]=0;
  
  //print parameters read from the file
  if(p->verbose>=0)
    {
      if(strcmp(p->expDataName,"")==0)
        {
          printf("ERROR: No experiment data file specified in the parameter file!\n");
          exit(-1);
        }
      else  
        printf("%i line(s) read from the parameter file: %s\n",numOpt,fileName);
      for(i=0;i<p->numSimData;i++)
        {
          printf("Taking simulated data from file (%i of %i): %s\n",i+1,p->numSimData,p->simDataName[i]);
          if(p->simDataFixedAmp[i]==1)
            printf("Fixing scaling factor for this data to %lf\n",p->simDataFixedAmpValue[i][0]);
          if(p->simDataFixedAmp[i]==2)
          	{
		        	if(p->simDataCommonScaling[i]==1)
		          	printf("Fixing scaling factor for this data to a factor of %lf relative to the last fitted data.\n",p->simDataFixedAmpValue[i][0]);
		          else
		          	{
		          		for(j=0;j<p->numSpectra;j++)
		          			printf("Fixing scaling factor for this data to a factor of %lf relative to the last fitted data for spectrum %i.\n",p->simDataFixedAmpValue[i][j],p->spectrum[j]);
		          	}
            }
        }
      if(p->channelScaling!=1.)
				printf("Channel scaling factor of %lf will be used.\n",p->channelScaling);
      if(p->peakSearch==0)
        for(i=0;i<p->numSpectra;i++)
          printf("Will compare spectrum %i from channels %i to %i.\n",p->spectrum[i],p->startCh[i],p->endCh[i]);
      else
        {
          for(i=0;i<p->numSpectra;i++)
            printf("Will search for a peak in spectrum %i from channels %i to %i.\n",p->spectrum[i],p->startCh[i],p->endCh[i]);
          if(p->peakSearchWidth>0)
            printf("Will set fitting window width to %i channels around each peak found.\n",p->peakSearchWidth);
        }
      if(p->commonScaling==1)
        printf("Will use common scaling and background for all spectra in each data file.\n");
      if(p->forcePositiveS==1)
      	printf("Will force scaling factors to positive values when fiting.\n");
      if(p->addBackground==0)
        printf("Will not add background to simulated data.\n");
      if(p->addBackground==1)
        printf("Will add a constant background to simulated data.\n");
      if(p->addBackground==2)
        printf("Will add a linear background to simulated data.\n");
      if(p->addBackground==3)
        printf("Will add a quadratic background to simulated data.\n");
      if(p->addBackground==1)
        for(i=0;i<p->numSpectra;i++)
          if(p->fixBG[i]==1)
            printf("Fixing background amplitude to %lf for spectrum %i, channels %i to %i.\n",p->fixedBGPar[i][0],p->spectrum[i],p->startCh[i],p->endCh[i]);
      if(p->addBackground==2)
        for(i=0;i<p->numSpectra;i++)
          if(p->fixBG[i]==1)
            printf("Fixing background parameters to A = %lf, B = %lf for spectrum %i, channels %i to %i.\n",p->fixedBGPar[i][0],p->fixedBGPar[i][1],p->spectrum[i],p->startCh[i],p->endCh[i]);
      if(p->addBackground==3)
        for(i=0;i<p->numSpectra;i++)
          if(p->fixBG[i]==1)
            printf("Fixing background parameters to A = %lf, B = %lf, C = %lf for spectrum %i, channels %i to %i.\n",p->fixedBGPar[i][0],p->fixedBGPar[i][1],p->fixedBGPar[i][2],p->spectrum[i],p->startCh[i],p->endCh[i]);
      if(p->indSpectra==1)
      	printf("Will treat spectra as independent measurements for reduced chisq calculation.\n");
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
