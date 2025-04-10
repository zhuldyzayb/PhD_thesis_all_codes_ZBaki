#include <stdio.h>
#include <math.h>
#define PI 3.1415926535

double lambda(int nl, double *sv1, double *sv2, double sp1, double sp2, double h){
     /* sv1 & sv2=sample vector (observations)
      sp1, sp2= sample point 
      h=bandwidth (scalar)
     */
     double val=0;
     int i;

    for(i=0; i<nl; i++){
        val+=exp(-0.5*( pow(sv1[i]-sp1,2) + pow(sv2[i]-sp2,2) )/pow(h,2));
    }
    return(val/(2*PI*pow(h,2)));
}

 double weightFunc(int np, double *gp1, double *gp2, double p1, double p2, double h, double *gridSize){
     /* gp1 & gp2=UTM coordinates of grid x (vector)
      p1, p2=UTM coordinates of observed y (point)
      h= global bandwidth (already selected in the previous stage)
     */
     double val=0;
     int i;
     double gg=gridSize[0]*gridSize[1];
     
    for(i=0; i<np; i++){
        val+=gg*exp(-0.5*( pow(gp1[i]-p1,2) + pow(gp2[i]-p2,2) )/pow(h,2));
    }
    return(val/(2*PI*pow(h,2)));
       
} // checked with R: correct

double lambdaCorr(int nl, int np, double *sp1, double *sp2, double y1, double y2, double *gp1, double *gp2, double h, double *gridSize){
     /* sp1 & sp2=sample vector (observations)
      gp1, gp2= grid points 
      h=global bandwidth
     */

    double val=0;
    int i;

    for (i=0; i<nl; i++){
        val+=exp(-0.5*( pow(sp1[i]-y1,2) + pow(sp2[i]-y2,2) )/pow(h,2))/weightFunc(np, gp1, gp2, sp1[i], sp2[i], h, gridSize);
        }
    return(val/(2*PI*pow(h,2)));

}  /// checked with R: correct

// Retrieve the product of edge-corrected weights separately


double myCum(int nl, int np, double *sp1, double *sp2, double *gp1, double *gp2, double h, double *gridSize, double power){
    double Cum=1;
    int i;

    for(i=0; i<nl; i++){
        Cum*=lambdaCorr(nl, np, sp1, sp2, sp1[i], sp2[i], gp1, gp2, h, gridSize);
    }
return(pow(Cum, power));
} /// checked with R: correct


void retrieveC(int *nl, int *np, double *sp1, double *sp2, double *gp1, double *gp2, 
               double *gridSize, double *hGlob, double *power, double *ans){
    int nLocs=*nl;
    int nPoints=*np;
    double hG=*hGlob;
    //double gs=*gridSize;
    double P=*power;
    double CC=0;
    int i;

    CC=myCum(nLocs, nPoints, sp1, sp2, gp1, gp2, hG, gridSize, P);

    for(i=0; i<nLocs; i++){
        ans[i]=pow(lambdaCorr(nLocs,nPoints, sp1, sp2, sp1[i], sp2[i], gp1, gp2, hG, gridSize)/CC, -0.5);
    }   

}







 
