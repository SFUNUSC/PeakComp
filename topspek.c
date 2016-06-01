//definitions
#include "topspek.h"
//functions and logic
#include "chisq.c"
#include "fitter.c"
#include "plotter.c"
#include "peak_window.c"
#include "read_parameters.c"
#include "read_data.c"
#include "read_and_process_data.c"
#include "save_data.c"

int main(int argc, char *argv[])
{
 
  //set up handler to take action upon SIGINT (CTRL-C command)
  struct sigaction sigIntHandler;
  sigIntHandler.sa_handler = sigint_cleanup;
  sigaction(SIGINT, &sigIntHandler, NULL);

  if(argc!=2)
    {
      printf("\ntopspek parameter_file\n----------------------\n\n");
      printf("Compares the spectra designated in the parameter file specified and generates cool statistics.\n\n");
      exit(-1);
    }
  
  //allocate structures
  par *p=(par*)malloc(sizeof(par));
  fitpar *fp=(fitpar*)malloc(sizeof(fitpar));
  data *d=(data*)malloc(sizeof(data));
  fitdata *fd=(fitdata*)malloc(sizeof(fitdata));

  readParFile(argv[1],p); //grab data from the parameter file (see read_config.c)

  readAndProcessData(p,d);//read and prepare data for fitting (see read_and_process_data.c)

  computeBackgroundandScaling(p,d,fp);//get background coefficients and scaling factors (see fitter.c)
   
  applyBackgroundandScaling(p,fp,d,fd);//apply fit parameters to data (see fitter.c)

  compareSpectra(p,d,fd);//generate chisq stats (see chisq.c)
  
  if((p->plotOutput>=1)&&(p->verbose>=0))
    plotSpectra(p,d,fd);//see plotter.c
  
  if(p->saveOutput==1)
    saveSpectra(p,fd);//see save_data.c
  
  //free structures
  free(p);
  free(fp);
  free(d);
  free(fd);

  return 0; //great success
}
