#ifndef RP_H
#define RP_H

#define MAX_OPT_PAR      100
#define MAX_NUM_OPT      200

//structure defining optional parameters to be read in from the parameter file
typedef struct
{
  char name[256];//name of the option
  char par[MAX_OPT_PAR][256];//parameters in the option
  int numPar;//number of parameters
}par_option; //fit parameters

#endif
