#include <stdio.h>
#include <math.h>
#define PI 3.1415926535

 void weightFuncNew(int *np, double *gp1, double *gp2, double *p1, double *p2, double *h, double *c, double *gridSize, double *ans){
     /* gp1 & gp2=UTM coordinates of grid x (vector)
      p1, p2=UTM coordinates of observed y (point)
      h= adaptive bandwidth
      c=\hat{c}(y, h_g)
     */
    int nPoints=*np;
    double y1=*p1;
    double y2=*p2;
    double hA=*h;
    double gs=*gridSize;
    double val=0;
    int i;
     
    for(i=0; i<nPoints; i++){
        val+=gs*exp( -0.5*( pow(gp1[i]-y1,2) + pow(gp2[i]-y2,2) )/pow(c[i]*hA,2) )/(2*PI*pow(c[i]*hA,2));
    }
    
    ans[0]=val;
       
} // checked with R: correct

