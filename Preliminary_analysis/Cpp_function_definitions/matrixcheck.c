#include <stdio.h>
#include <math.h>
//calculates the sum R(s,t)\delta s for fixed values of t

double myGamma(double p1, double p2, double time, double *a, double c, double asigma, double Bval, double c2, double c3){
    double temp=0;
    temp=(Bval+c3)*exp((-1)*c*(a[0]+a[1]*time+a[2]*p1+a[3]*p2+a[4]*time*time+a[5]*p1*p1+a[6]*p2*p2+a[7]*time*p1+a[8]*time*p2+a[9]*p1*p2)/asigma)/c2;
    return(temp);  
}
//checked with R: correct

double myRiskFunc(double p1, double p2, double time, double *a, double *theta, double g1, double g2, double Bval){
    //theta=[r, tau, Asigma,c, c2,c3, b0, b1, b2]
    double finval=0;
    finval=exp(theta[6]+theta[7]*g1+theta[8]*g2)*theta[0]/(theta[1]*myGamma(p1,p2,time,a,theta[3], theta[2], Bval, theta[4], theta[5]));
    return(finval);
}
//checked with R: correct
 
void myRiskInt(int *nl, double *s1, double *s2, double *t, double *a, double *theta, double *gas1, double *gas2, double *B, double *ans){
    //theta=[r, tau, Asigma,c, c2, c3, b0, b1, b2]
    //s1 & s2 =vectors of locations
    //t is time point value
    //gas1 is a vector of spatial values of annual gas production at time t in all locations s
    //gas2 is a single cumulative value
    int nLocs=*nl;
    double time=*t;
    double Bval=*B;
    double gcum=*gas2;
    double temp=0;
    int i;
    for(i=0; i<nLocs; i++){
        temp+=0.25*myRiskFunc(s1[i], s2[i],time, a, theta, gas1[i], gcum, Bval);
    }
    ans[0]=temp;
}
//checked with R: correct