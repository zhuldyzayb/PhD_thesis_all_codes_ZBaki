#include <R.h>
#include <math.h>

void AbramsonPixels( int *nx, double *xmax, int *ny, double *ymax, 
                     int *nt, double *tmax, int *window,
                     int *npoints, double *px, double *py, double *pt,
                     double *pcoef, double *ans, double *hs, double *ht,
                     double *alpha )
{

/* 
   Voxel image of nx by ny pixels representing [0, xmax] x [0, ymax] 
   with time in [0, tmax]
   Spatial window indicator window

   It is written to a matrix, nrow=ny, ncol=nx, by row. Hence
   for i in 1, ..., nx and j = 1, ..., ny
      the (j-1)*nx + i th member of ans corresponds to matrix[j,i]
      or in Euclidean terms to
         ( ( i - 0.5 ) * xmax / nx, ( j - 0.5 ) * ymax / ny )
   In C, the first element is labelled 0

   pcoef contains density estimators for the points

   Gaussian kernels in space and time with local edge correction

*/

   int nxCoord = *nx;
   int nyCoord = *ny;
   int ntCoord = *nt;

   int nIm = nxCoord * nyCoord;
   double maxx = *xmax;
   double maxy = *ymax;
   double maxt = *tmax;

   double xstep = maxx / nxCoord;
   double ystep = maxy / nyCoord;
   double tstep = maxt / ntCoord;

/*
  Point pattern of length npoints, coordinates (px, py, pt) 
  pcoef contains pilot values of intensity
*/
  
   int nPoints = *npoints;
   double h2 = *hs; 
   h2 = h2 * h2;
   double h1 = *ht;
   double abram = *alpha;

   int k, i, j, t;
   double xi, yj, tt;
   double dist2, kappa2; 
   double dist1, kappa1;
   double geomMean, bw2;
   double edge[nPoints];
   for( k=0; k < nPoints; k++ ) edge[k] = 0.0;
   double c = exp( log(2 * M_PI) * 3 / 2 );
   double cAbram;

/* Calculate normalisation constant */

   geomMean = 0.0;
   for( k=0; k<nPoints; k++ ) geomMean += log( pcoef[k] );
   geomMean = exp( geomMean / nPoints);

/* Fill edge correction weight vector */

   for( k=0; k<nPoints; k++ ) {
      bw2 = exp( log(pcoef[k] / geomMean) * 2 * abram );

      for( t=0; t<ntCoord; t++ ) {
         tt = (t+0.5) * tstep;

         for( i=0; i<nxCoord; i++ ) { 
            xi = (i+0.5) * xstep; 
 
            for( j=0; j<nyCoord; j++ ) {
               yj = (j+0.5) * ystep;

               if( window[j*nxCoord + i] == 1 ) {

               dist2 = ( px[k] - xi ) * ( px[k] - xi ) + ( py[k] - yj ) * ( py[k] - yj ); 
               kappa2 = exp( - dist2 / ( 2 * h2 * bw2 ) );

               dist1 = ( pt[k] - tt ) * ( pt[k] - tt );
               kappa1 = exp( - dist1 / ( 2 * h1 * h1 * bw2 ) );

               edge[k] += kappa1 * kappa2;
               if( isnan( edge[k] ) ) { 
                 fprintf(stderr, " k = %d  bw2 = %f\n", k, bw2 ); 
                 exit(1); 
               }
               }
            }
         }
      }
      cAbram = c * exp( log( bw2 ) * 3 / 2 );
      edge[k] = ( ( xstep * ystep * tstep ) * edge[k] ) / ( h1 * h2 * cAbram );
      //fprintf(stderr, "edge k = %f\n", edge[k] ); 
   }

/* Initialise */

   for( k=0; k < nIm*ntCoord; k++ ) ans[k] = 0.0;

/* Fill */

   for( t=0; t<ntCoord; t++ ) {
      tt = (t+0.5) * tstep;

      for( i=0; i<nxCoord; i++ ) {
         xi = (i+0.5) * xstep;

         for(j=0; j<nyCoord; j++ ) {
            yj = (j+0.5) * ystep;

            if( window[j*nxCoord + i] == 1 ) {
            for( k=0; k<nPoints; k++ ) {
               bw2 = exp(  log(pcoef[k] / geomMean) * 2 * abram );
               dist2 = ( px[k] - xi ) * ( px[k] - xi ) + ( py[k] - yj ) * ( py[k] - yj );
               kappa2 = exp( - dist2  / ( 2 * h2 * bw2 ) ); 

               dist1 = ( pt[k] - tt ) * ( pt[k] - tt );
               kappa1 = exp( - dist1  / ( 2 * h1 * h1 * bw2 ) ); 

               ans[t*nxCoord*nyCoord + j*nxCoord + i] += 
                  ( kappa1 * kappa2 ) * exp( - 3 * log(bw2) / 2 ) / edge[k];

               if( isnan( ans[t*nxCoord*nyCoord + j*nxCoord + i] ) ) { 
                 fprintf(stderr, " i = %d, j = %d, bw2 = %f\n", i, j, bw2 ); 
                 exit(1);
               }
            }
            }
         }
      }
   }  


   for( k=0; k<nIm*ntCoord; k++ ) ans[k] = ans[k] / ( h1 * h2  * c );

}


void AbramsonPoints( int *npoints, double *px, double *py, double *pt,
                     double *pcoef, double *ans, double *hs, double *ht, 
                     double *alpha )
{

/*
  Space time point pattern of length npoints, coordinates (px, py, pt) 
  hs the spatial bandwidth, ht the temporal one; 
  Gaussian kernels in space and time, no edge correction;
  pcoef contains pilot values of intensity; alpha = - 0.5 for classic case
*/
  
   int nPoints = *npoints;
   double h2 = *hs; 
   h2 = h2 * h2;
   double h1 = *ht;
   double abram = *alpha;

   int k, i;
   double xi, yi, ti;
   double dist2, kappa2;
   double dist1, kappa1;
   double geomMean, bw2;

/* Initialise */

   for( k=0; k < nPoints; k++ ) ans[k] = 0.0;

/* Calculate normalisation constant */

   geomMean = 0.0;
   for( k=0; k<nPoints; k++ ) geomMean += log( pcoef[k] );
   geomMean = exp( geomMean / nPoints);

/* Fill */

   for( i=0; i<nPoints; i++ ) {
      xi = px[i];
      yi = py[i];
      ti = pt[i];

         for( k=0; k<nPoints; k++ ) {
            bw2 = exp(  log(pcoef[k] / geomMean) * 2 * abram );

            dist2 = ( px[k] - xi ) * (px[k] - xi ) + ( py[k] - yi ) * ( py[k] - yi );
            kappa2 = exp( - dist2  / ( 2 * h2 * bw2 ) ); 

            dist1 = ( pt[k] - ti ) * ( pt[k] - ti );
            kappa1 = exp( - dist1 / ( 2* h1 * h1 * bw2 ) );

            ans[i] += ( kappa1 * kappa2 ) * exp( - 3 * log(bw2) / 2 );
            if( isnan( ans[i] ) ) { 
              fprintf(stderr, " i = %d  bw2 = %f\n", i, bw2 ); 
              exit(1); 
            }
         }
   }  

   double c = exp( log( 2 * M_PI ) * 3 / 2 );
   for( k=0; k<nPoints; k++ ) ans[k] = ans[k] / ( h1 * h2 * c );
}


void IntensityPixels( int *nx, double *xmax, int *ny, double *ymax, 
                      int *nt, double *tmax, int *window,
                      int *npoints, double *px, double *py, double *pt,
                      double *ans, double *hs, double *ht )
{

/* 
   Voxel image of nx by ny pixels representing [0, xmax] x [0, ymax] 
   with time coordinate in [0, tmax]

   Spatial window indicator window
   Writtin to a matrix, nrow=ny, ncol=nx, by row. Hence, 
   for i in 1, ..., nx and j = 1, ..., ny, 
      the (j-1)*nx +i the member of ans corresponds to matrix[j,i]
      or in Euclidean terms to 
      ( ( i - 0.5 ) * xmax / nx, ( j - 0.5 ) * ymax / ny )
   In C, the first element is labelled 0

   Gaussian kernels in space and time, local edge correction factor

*/

   int nxCoord = *nx;
   int nyCoord = *ny;
   int ntCoord = *nt;

   int nIm = nxCoord * nyCoord;
   double maxx = *xmax;
   double maxy = *ymax;
   double maxt = *tmax;

   double xstep = maxx / nxCoord;
   double ystep = maxy / nyCoord;
   double tstep = maxt / ntCoord;

/*
  Point pattern of length npoints, coordinates (px, py, pt) 
*/
  
   int nPoints = *npoints;
   double h2 = *hs;
   h2 = h2 * h2;
   double h1 = *ht;

   int k, i, j, t;
   double xi, yj, tt;
   double dist2, kappa2;
   double dist1, kappa1;

   double c = exp( log(2 * M_PI) * 3 / 2 );
   double edge[nPoints];
   for( k=0; k < nPoints; k++ ) edge[k] = 0.0;

/* Fill edge correction weight vector */

   for( k=0; k<nPoints; k++ ) {

      for( t=0; t<ntCoord; t++ ) {
         tt = (t+0.5) * tstep;

         for( i=0; i<nxCoord; i++ ) {
            xi = (i+0.5) * xstep;

            for( j=0; j<nyCoord; j++ ) {
               yj = (j+0.5) * ystep;

            if( window[j*nxCoord+i] == 1 ) {
               dist2 = ( px[k] - xi ) * ( px[k] - xi ) + ( py[k] - yj ) * ( py[k] - yj );
               kappa2 = exp( - dist2 / ( 2 * h2 ) );

               dist1 = ( pt[k] - tt ) * ( pt[k] - tt );
               kappa1 = exp( - dist1 / ( 2 * h1 * h1 ) );

               edge[k] += kappa1 * kappa2;
               if( isnan( edge[k] ) ) {
                 fprintf(stderr, " k = %d  edge[k] = %f\n", k, edge[k] );
                 exit(1);
               }
            }
            }
         }
      }

      edge[k] = ( ( xstep * ystep * tstep ) * edge[k] ) / ( h1 * h2 * c );
      //fprintf(stderr, "edge k = %f; ", edge[k] );
   }


/* Initialise */

   for( k=0; k < nIm*ntCoord; k++ ) ans[k] = 0.0;

/* Fill */

   for( t=0; t<ntCoord; t++ ) {
      tt = (t+0.5) * tstep;

      for( i=0; i<nxCoord; i++ ) {
         xi = (i+0.5) * xstep;

         for( j=0; j<nyCoord; j++ ) {
            yj = (j+0.5) * ystep;
      
            if( window[j*nxCoord+i] == 1 ) {

            for( k=0; k<nPoints; k++ ) {
               dist2 = ( px[k] - xi ) * (px[k] - xi ) + ( py[k] - yj ) * ( py[k] - yj );
               kappa2 = exp( - dist2  / ( 2 * h2 ) ); 

               dist1 = ( pt[k] - tt ) * (pt[k] - tt );
               kappa1 = exp( - dist1  / ( 2 * h1 * h1 ) ); 

               ans[t*nyCoord*nxCoord + j*nxCoord + i] += 
                  ( kappa1 * kappa2 / edge[k] ); 
               if( isnan( ans[t*nyCoord*nxCoord + j*nxCoord + i] ) ) { 
                 fprintf(stderr, " i = %d, j = %d, t = %d\n", i, j, t ); 
                 // fprintf(stderr, " k = %d, edge[k] = %f\n", k, edge[k] ); 
                 // exit(1); 
               }
            }
         }
         }
      }
   }  

   for( k=0; k<nIm*ntCoord; k++ ) ans[k] = ans[k] / ( h1 * h2 * c );
}


void IntensityPoints( int *npoints, double *px, double *py, double *pt,
                      double *ans, double *hs, double *ht )
{

/*
  Space time point pattern of length npoints, coordinates (px, py, pt) 
  hs the spatial bandwidth, ht the temporal ones
  Gaussian kernels in space and time; no edge correction
*/
  
   int nPoints = *npoints;
   double h2 = *hs;
   h2 = h2 * h2;
   double h1 = *ht;

   int k, i;
   double xi, yi, ti;
   double dist2, kappa2;
   double dist1, kappa1;

/* Initialise */

   for( k=0; k<nPoints; k++ ) ans[k] = 0.0;


/* Fill */

   for( i=0; i<nPoints; i++ ) {
      xi = px[i];
      yi = py[i];
      ti = pt[i];

         for( k=0; k<nPoints; k++ ) {
            dist2 = ( px[k] - xi ) * ( px[k] - xi ) + ( py[k] - yi ) * ( py[k] - yi );
            kappa2 = exp( - dist2  / ( 2 * h2 ) ); 

            dist1 = ( pt[k] - ti ) * ( pt[k] - ti );
            kappa1 = exp( - dist1 / ( 2 * h1 * h1 ) );

            ans[i] +=  kappa1 * kappa2;
            if( isnan( ans[i] ) ) { 
              fprintf(stderr, " i = %d  kappa1 = %f, kappa2 = %f\n", i,
                 kappa1, kappa2 ); 
              exit(1); 
            }
         }
   }  

   double c = exp( log(2 * M_PI) * 3 / 2 );
   for( k=0; k<nPoints; k++ ) ans[k] = ans[k] / ( h1 * h2 * c );

}

