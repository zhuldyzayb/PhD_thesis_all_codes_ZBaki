#include <stdio.h>    
#include <math.h>
#include <complex.h>
# define PI 3.1415926535


double Covariance(double *point1, double *point2, double *theta){
// theta=[w, C_0, C, a] for spherical variogram
// point1 and point2 are 1*2 vectors of coordinates
double h;
double gamma;
h=sqrt(pow(point1[0]-point2[0],2)+pow(point1[1]-point2[1],2));
if(h==0){
    gamma=0;
}
else if (h<=theta[3]){
    gamma=theta[1]+theta[2]*(1.5*h/theta[3]-pow(h,3)/(2*pow(theta[3],3)));
}
else {gamma=theta[1]+theta[2];}
return(theta[0]-gamma/2);
}//checked:correct


int main()
{ 
    double vec[4][2]={{3.5,0},{3,0.5}, {3,2}, {2,0.5}}; 
    double mytheta[4]={2,2.3,4.1,4};
    int n=3;
    double Sigma[n][n];
    int i,j,k;
    
    
    // define upper part of Sigma
    //for(i=0; i<n; i++){
      //  for(j=i;j<n; j++){
        //    double *p1=vec[i];
        //    double *p2=vec[j];
        //    Sigma[i][j]=Covariance(p1, p2,mytheta);
        //}
    //}

    Sigma[0][0]=1;
    Sigma[0][1]=2;
    Sigma[0][2]=3;
    Sigma[1][1]=1;
    Sigma[1][2]=2;
    Sigma[2][2]=1;
    //mirror the lower part
    for(i=1; i<n; i++){
        for(j=0;j<i; j++){
            Sigma[i][j]=Sigma[j][i];
        }
    }

    // define Fourier matrix F
    double complex F[n][n];
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            if (i==0 || j==0){F[i][j]=1/pow(n,0.5);}
            else{
                k=i*j;
                F[i][j]=cos(2*PI*k/n)/pow(n,0.5)-sin(2*PI*k/n)/pow(n,0.5)*I;
            }
        }
    }
    
    // define Hermitian of Fourier matrix F^H
    double complex Fh[n][n];
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            Fh[i][j]=creal(F[j][i])-cimag(F[j][i])*I;
        }
    }

    //define Lambda matrix of eigenvalues L:
    double complex first[n][n];
    for(i=0;i<n;i++){    
        for(j=0;j<n;j++){    
            first[i][j]=0;    
            for(k=0;k<n;k++){    
                first[i][j]+=Sigma[i][k]*F[k][j];  
                }    
            }    
        }    
    double complex L[n][n];
    for(i=0;i<n;i++){    
        for(j=0;j<n;j++){    
            L[i][j]=0;    
            for(k=0;k<n;k++){    
                L[i][j]+=Fh[i][k]*first[k][j];  
                }    
            }    
        }
    
    
    printf("My sigma matrix: \n");
    for (int i = 0; i < n; ++i){
    for (int j = 0; j < n; ++j)
    {
      printf("%.2f\t", Sigma[i][j]);
    }
        printf("\n");
    }
    printf("My Fourier matrix: \n");
    for (int i = 0; i < n; ++i){
    for (int j = 0; j < n; ++j)
    {
      printf("%.2f+%.2f*i\t", creal(F[i][j]), cimag(F[i][j]));
    }
        printf("\n");
    }
    printf("My Hermitian of Fourier matrix: \n");
    for (int i = 0; i < n; ++i){
    for (int j = 0; j < n; ++j)
    {
      printf("%.2f+%.2f*i\t", creal(Fh[i][j]), cimag(Fh[i][j]));
    }
        printf("\n");
    }
    
    printf("My Lambda matrix: \n");
    for (int i = 0; i < n; ++i){
    for (int j = 0; j < n; ++j)
    {
      printf("%.2f+%.2f*i\t", creal(L[i][j]), cimag(L[i][j]));
    }
        printf("\n");
    }
    


    return 0;
}