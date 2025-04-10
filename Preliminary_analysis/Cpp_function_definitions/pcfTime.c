#include <math.h>
#include <stdio.h>
#include <stdlib.h>


void pcfTimeSmooth( int *nlags, double *lagmax, 
          int *npoints, double *pt, double *maxT, 
          double *pcoef, double *ans, double *eps )
{

/* 
 
 Point pattern of length npoints, coordinates (px, py) and 
 intensity estimators pcoef.
 Pcf density estimator with translation edge correction 
 and Gaussian kernel standard error epsilon

 Time interval assumed to be [0, T];
 
*/
   int nPoints = *npoints;
   int nLags = *nlags;
   double maxLag = *lagmax;
   double maxTime = *maxT;
   double stepSize = maxLag / nLags;
   double epsilon = *eps;

   int i, j, k;
   double dist, exponent;
   double tk;

/* Initialise */

   for( k=0; k <= nLags; k++ ) ans[k] = 0.0;

/* Fill */

   for( k=0; k <= nLags; k++ ) { 
     tk =  k * stepSize;
     for( i=0; i<nPoints; i++ ) {
       for(j=0; j<nPoints; j++ ) {
         if( i != j ) { 
           dist = pt[j] - pt[i];
           exponent =  ( tk - ( pt[j] - pt[i] ) ) / epsilon;
           ans[k] += exp( - exponent * exponent / 2.0 ) / ( ( pcoef[i] * pcoef[j] ) * ( maxTime-fabs(dist) ) ); }
       }
     }
   }  

   for( k=0; k<=nLags; k++ ) ans[k] = ans[k] / ( sqrt( 2 * M_PI ) * epsilon );
}

void pcfTime( int *nlags, double *lagmax, 
          int *npoints, double *pt, double *maxT, 
          double *pcoef, double *ans, double *eps )
{

/* 
 
 Point pattern of length npoints, coordinates (px, py) and 
 intensity estimators pcoef.
 Pcf density estimator with translation edge correction 
 and box kernel radius epsilon

 Time interval assumed to be [0, T];
 
*/
   int nPoints = *npoints;
   int nLags = *nlags;
   double maxLag = *lagmax;
   double maxTime = *maxT;
   double stepSize = maxLag / nLags;
   double epsilon = *eps;

   int i, j, k, kmin, kmax;
   double dist;

/* Initialise */

   for( k=0; k <= nLags; k++ ) ans[k] = 0.0;

/* Fill */

   for( i=0; i<nPoints; i++ ) {
      for(j=0; j<nPoints; j++ ) {
         if( i != j ) { 
           dist =  pt[j] - pt[i];

           kmin = ceil( ( dist - epsilon ) / stepSize );
           if( kmin < 0 ) kmin = 0;
           kmax = floor( ( dist + epsilon ) / stepSize );
           if( kmax > nLags ) kmax = nLags;

/* fprintf(stderr, " dist = %f, kmin = %d, kmax = %d\n", dist, kmin, kmax ); */

           for( k = kmin; k <= kmax; k++ ) { 
              ans[k] += 1.0 / ( ( pcoef[i] * pcoef[j] ) * ( maxTime - fabs(dist) ) ); }
         }
      }
   }  


   for( k=0; k<=nLags; k++ ) ans[k] = ans[k] / ( 2 * epsilon );
}

void IntensityTime( int *npoints, double *pt, double *maxT, 
                    double *ans, double *h )
{

/* 
 
 Point pattern of length npoints, coordinates (px, py) 
 Gaussian kernel standard deviation h 
 Returns classic kernel density estimator without
 edge correction and using all the points (leave-one-out is FALSE)
 
*/
   int nPoints = *npoints;
   double h2 = *h;
   h2 = h2 * h2;

   int i, j, k;
   double dist2, kappa;


/* Initialise */

   for( k=0; k < nPoints; k++ ) ans[k] = 0.0;

/* Fill */

   for( i=0; i<nPoints; i++ ) {
      for(j=0; j<nPoints; j++ ) {
            dist2 = ( pt[i] - pt[j] ) * (pt[i] - pt[j] );
            kappa = exp( - dist2  / ( 2 * h2 ) ); 
            ans[i] += kappa; 
            if( isnan( ans[i] ) ) { 
              fprintf(stderr, " i = %d, j = %d, kappa = %f\n", i, j, kappa ); 
              exit(1); 
            }  
       }
   }

   for( k=0; k<nPoints; k++ ) ans[k] = ans[k] / sqrt( h2 * 2.0 * M_PI );

}


void IntensityBins( int *ntimes, int *npoints, double *pt, double *wts, 
                    double *maxT, double *ans, double *h )
{

/*

 Point pattern of length npoints, coordinates (px, py)
 Gaussian kernel standard deviation h
 Returns classic weighted kernel density estimator without
 edge correction at ntimes pixels  between 0 and maxT

*/
   int nPoints = *npoints;
   int nTimes = *ntimes;
   double maxTime = *maxT;
   double h2 = *h;
   h2 = h2 * h2;

   int i, j, k;
   double dist2, kappa;
   double ti;

/* Initialise */

   for( k=0; k<nTimes; k++ ) ans[k] = 0.0;

/* Fill */
   for( i=0; i<nTimes; i++ ) {
     ti =  (i+0.5) * (maxTime/nTimes);
         for(j=0; j<nPoints; j++ ) {
            dist2 = ( ti - pt[j] ) * ( ti - pt[j] );
            kappa = exp( - dist2  / ( 2 * h2 ) );
            ans[i] += kappa * wts[j];
            if( isnan( ans[i] ) ) {
              fprintf(stderr, " i = %d, j = %d, kappa = %f\n", i, j, kappa );
              exit(1);
            }
       }
   }

   for( k=0; k<nTimes; k++ ) ans[k] = ans[k] / sqrt( h2 * 2 * M_PI );
}

