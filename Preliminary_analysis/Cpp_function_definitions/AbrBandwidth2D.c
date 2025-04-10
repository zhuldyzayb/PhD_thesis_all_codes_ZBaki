#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define PI 3.1415926535
/// The functions are checked using the small sample size in R and calculation results are correct.

 double lambda(int nl, double *sv1, double *sv2, double sp1, double sp2, double h){
     /* sv1 & sv2=UTM coordinates of y
      sp1, sp2=x_0 coordinates (from sample data)
      h=bandwidth (scalar)
     */
     double val=0;
     int i;

    for(i=0; i<nl; i++){
        val+=exp(-0.5*( pow(sv1[i]-sp1,2) + pow(sv2[i]-sp2,2) )/pow(h,2));
    }
    return(val/(2*PI*pow(h,2)));
}


 void firstOpt(int *nl, double *sv1, double *sv2, double *h, double *lw, double *ans){
    int nLocs=*nl;
    double band=*h;
    double Measure=*lw;
    int i;
    double val=0;

    for(i=0; i<nLocs; i++){
        val+=1/lambda(nLocs, sv1, sv2, sv1[i], sv2[i], band);
    }
   ans[0]=fabs(val-Measure);
 }