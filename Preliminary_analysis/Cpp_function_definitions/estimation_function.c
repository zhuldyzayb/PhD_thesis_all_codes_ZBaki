//  Created by Zhuldyzay Baki on 01/02/2023.
//  This file is created in order to speed up the calculations of the parameter estimates
//  for theta=[\gamma_0, alpha, eta, theta_1, theta_2] using Waagepetersen estimation function

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// FUNCTIONS
// M(s,t)
void findMst(double *a, double *b, double *c, double *beta, double *ans){
    double x=*a;
    double y=*b;
    double t=*c;
    double res;
    x=x-750.0;
    y=y-5900.0;
    res=beta[0]+beta[1]*t+beta[2]*pow(t,2)+beta[3]*x+beta[4]*y+beta[5]*pow(x,2)+beta[6]*x*y+beta[7]*pow(y,2)+beta[8]*pow(x,3)+beta[9]*pow(x,2)*y+beta[10]*x*pow(y,2)+beta[11]*pow(y,3)+beta[12]*pow(x,4)+beta[13]*pow(x,3)*y+beta[14]*pow(x,2)*pow(y,2)+beta[15]*x*pow(y,3)+beta[16]*pow(y,4)+beta[17]*t*x+beta[18]*t*y+beta[19]*t*pow(x,2)+beta[20]*t*x*y+beta[21]*t*pow(y,2)+beta[22]*t*pow(x,3)+beta[23]*t*pow(x,2)*y+beta[22]*t*x*pow(y,2)+beta[25]*t*pow(y,3);
    ans[0]=res;
}
