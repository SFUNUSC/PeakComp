#include "peak_find.h"

//function to find a peak, given an approximate value for its centroid
peak_fit_par findPeak(peak_search_par * spar, int * data, int dataLength)
{ 
  //set default return parameters
  peak_fit_par par;
  par.centroid=(spar->searchMin + spar->searchMax)/2;
  par.width=spar->searchMax-spar->searchMin;
  
  double win1val,win2val,tfval;//holds average values for trapezoidal filter window(s)
  double maxval,minval;//holds maximum and minimum values from the filter
  int maxch,minch;//holds indices for maximum and minimum values from the filter
  int i,j,maxFlag,minFlag;
  
  
  //set bounds for trapezoidal filter
  int startCh=spar->searchMin;
  int endCh=spar->searchMax;

  maxval=0;
  minval=BIG_NUMBER;
  minFlag=0;
  maxFlag=0;
  for(i=startCh;i<endCh;i++)
    if((i+2*spar->windowSize)<dataLength)//check that we won't go out of bounds
      {
        win1val=0;
        win2val=0;
        for(j=0;j<spar->windowSize;j++)
          {
            win1val+=data[i+j];
            win2val+=data[i+j+spar->windowSize];
          }
        win1val/=spar->windowSize;//get averages
        win2val/=spar->windowSize;
        tfval=win2val-win1val;
        if(tfval>maxval)
          {
            maxval=tfval;
            maxch=i+spar->windowSize;
            maxFlag=1;
          }
        if(tfval<minval)
          {
            minval=tfval;
            minch=i+spar->windowSize;
            minFlag=1;
          }
      }
  
  if((maxFlag==1)&&(minFlag==1))
    {
      par.width=abs(minch-maxch);
      par.centroid=(minch+maxch)/2;
    } 
  
  return par;
}
