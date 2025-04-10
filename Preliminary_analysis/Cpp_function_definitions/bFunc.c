#include <stdio.h>
#include <math.h>

void BfuncIntegral(int *nl, int *nt, double *gp1, double *gp2, double *gt, double *alphas, double *Asigma, double *cp, double *ans){
    /* nl = number of spatial locations
       gp1 & gp2 = vectors with UTM coordinates of grided locations
       nt = number of points in time
       gt = time points vector
       alphas = parameter estimates from m(s,t)
       Asigma = parameter
       cp = coefficient c from replacement S(s,t)=c*p(s,t)
    */
    //int nLocs=*nl;
    int nTime=*nt;
    double asigma=*Asigma;
    double c=*cp;
    int i;
    double temp=0;
     
    for(i=0; i<nTime; i++){
        temp=c*gt[i]*gt[i]/asigma;
        //for(j=0; j<nLocs; j++){
        //    temp+=0.5*exp(c*(alphas[0]+alphas[1]*gt[i]+alphas[2]*gp1[j]+alphas[3]*gp2[j]+alphas[4]*gt[i]*gt[i]+alphas[5]*gp1[j]*gp1[j]+alphas[6]*gp2[j]*gp2[j]+alphas[7]*gt[i]*gp1[j]+alphas[8]*gt[i]*gp2[j]+alphas[9]*gp1[j]*gp2[j])/asigma);
        //}
    }
    
    ans[0]=temp;
       
} 
// 