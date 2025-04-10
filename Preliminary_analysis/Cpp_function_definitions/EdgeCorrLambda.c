#include <math.h>
#include <stdio.h>
#include <stdlib.h>


double weightsFunc(double spt1, double spt2, double tpt, int np, int nt, 
                double *gpt1, double *gpt2, double *gt, double *h)
{
/*
This functio calculates thw weights in edge corrected lambda for every given obs point y (26 locations by 6209 days)

spt1 & spt2 = UTM coordinates of a point y
tpt = time coordinate of a point y
np = number of points in space
nt = number of points in time
gpt1 & gpt2 = gridded location UTM coordinates z(s,) (1018)
gt = gridded time coordinates z( ,t)
h = global bandwidth [h_S1, h_S2, h_T]

*/

double sdiff=0;
double tdiff=0;
double temp=0;
int i,j;

for(i=0; i<np; i++){
    sdiff=(gpt1[i]-spt1)*(gpt1[i]-spt1)/(h[0]*h[0])+(gpt2[i]-spt2)*(gpt2[i]-spt2)/(h[1]*h[1]);
    for(j=0; j<nt; j++){
        tdiff=(gt[j]-tpt)*(gt[j]-tpt)/(h[2]*h[2]);
        temp+=exp(-0.5*(sdiff+tdiff));
    }
}

return temp/(h[0]*h[1]*h[2]);

}

double correctedLambda(int nTime, int nLocs, int nPoints, double sPoint1, 
             double sPoint2, double tPoint, double *osp1, double *osp2, 
             double *ot, double *h, double *gpt1, double *gpt2, double *gt)  
{
/*  
This function calculates the edge corrected lambda_c( ) for any given x (1018 locs by 6209 days)
 
 sPoint1 & sPoint2 = UTM coordinates of x [1*1]
 tPoint = time coordinate of x [1*1]
 osp1 & osp2 = vectors with observed UTM locations [26*1] each
 ot = vector with observed time pts [6209*1]
 h = global bandwidth 
 gpt1 & gpt2 = gridded location UTM coordinates z(s,) (1018)
 gt = gridded time coordinates z( ,t)

*/

   int i, j;
   double td=0;
   double sd=0;
   double w=0;
   
/* Initialise */

   double temp=0;
   
/* Fill */

    for( i=0; i<nLocs; i++ ) {
         sd=(osp1[i]-sPoint1)*(osp1[i]-sPoint1)/(h[0]*h[0])+(osp2[i]-sPoint2)*(osp2[i]-sPoint2)/(h[1]*h[1]);
       for(j=0; j<nTime; j++ ) {
           td=(tPoint-ot[j])*(tPoint-ot[j])/(h[2]*h[2]);
           w=weightsFunc(osp1[i], osp2[i], ot[j], nPoints, nTime, gpt1, gpt2, gt, h);
           temp+=exp(-0.5*(sd+td))/w;
       } 
     }
    return(temp/(h[0]*h[1]*h[2]));
}


double cAdaptive(int nTime, int nLocs, int nPoints, double sPoint1, 
             double sPoint2, double tPoint,  int N, double *osp1, double *osp2, 
             double *ot, double *h, double *gpt1, double *gpt2, double *gt)
{
/*
This function calculates c(s,t) for from Abramson estimator
 
 sPoints1 & sPoints2 = UTM coordinates of choosen y (26 by 6209)
 tPoint = time coordinate of chosen y (26 by 6209)
 osp1, osp2, ot = UTM and tme values of all observed data y
 gpt1, gpt2, gt = UTM and time coordinates of gridded data x
 h = global bandwidth (selected from AbrGlobalBand)
 nTime = 6209
 nLocs = 26
 nPoints = 1018 
 N = N(\Psi \cap W_S\times W_T)

*/

int i,j;
double c=0;

double prod=1; 
for(i=0; i<nLocs;i++){
    for(j=0; j<nTime; j++){
        prod*=correctedLambda(nTime, nLocs, nPoints, osp1[i], osp2[i], ot[j], osp1, osp2, ot, h, gpt1, gpt2, gt);  
    }
}

c=pow(correctedLambda(nTime, nLocs, nPoints, sPoint1, sPoint2, tPoint, osp1, osp2, ot, h, gpt1, gpt2, gt)/(pow(prod, 1/N)), -0.5);
 return c;   
}

double AbramsonLambda(int nTime, int nLocs, int nPoints, double sPoint1, 
             double sPoint2, double tPoint, double *osp1, double *osp2, 
             double *ot, double *h, double *gpt1, double *gpt2, double *gt, int N)
{
/*
This function calculates Abramson estimator lambda_A
 
 sPoints1 & sPoints2 = UTM coordinates of choosen x_0 
 tPoint = time coordinate of chosen x_0 
 osp1, osp2, ot = UTM and tme values of all observed data y
 gpt1, gpt2, gt = UTM and time coordinates of gridded data x
 h = global bandwidth (selected from AbrGlobalBand)
 nTime = 6209
 nLocs = 26
 nPoints = 1018 
 N = N(\Psi \cap W_S\times W_T)

*/

int i, j;
double val=0; 
double sd=0;
double td=0;
double c=0; 

for(i=0; i<nLocs;i++){
    sd=(osp1[i]-sPoint1)*(osp1[i]-sPoint1)/(h[0]*h[0])+(osp2[i]-sPoint2)*(osp2[i]-sPoint2)/(h[1]*h[1]);
    for(j=0; j<nTime; j++){
        td=(ot[j]-tPoint)*(ot[j]-tPoint)/(h[2]*h[2]);
        c=cAdaptive(nTime, nLocs, nPoints, osp1[i],osp2[i], ot[j], N, osp1, osp2, ot, h,  gpt1,  gpt2, gt);
        val+= exp(-0.5*(td+sd)/(c*c));
        }
    }

return val/(pow(c,3)*h[0]*h[1]*h[2]);
}            



void newOptimObj(int *nt, int *nl, int *np, double *osp1, double *osp2, 
             double *ot, double *h, double *gpt1, double *gpt2, double *gt, int *num, double *lw1, double *lw2, double *ans)
{
/*
This is an objective function for minimization in adaptive bandwidth selection

 osp1, osp2, ot = UTM and tme values of all observed data y
 gpt1, gpt2, gt = UTM and time coordinates of gridded data x
 h = global bandwidth (selected from AbrGlobalBand)
 nTime = 6209
 nLocs = 26
 nPoints = 1018 
 num = N(\Psi \cap W_S\times W_T)

*/

int nTime=*nt;
int nLocs=*nl;
int nPoints=*np;
int N=*num; 
double LW1=*lw1;
double LW2=*lw2;

int i,j;
double value=0;

for(i=0; i<nPoints;i++){
    for (j=0; j<nTime; j++){
        value+=1/AbramsonLambda(nTime, nLocs, nPoints, gpt1[i],gpt2[i], gt[j], osp1, osp2, ot, h, gpt1, gpt2, gt, N);
    }
}

for(i=0; i<1; i++){ans[i]=value-LW1*LW2;}
}