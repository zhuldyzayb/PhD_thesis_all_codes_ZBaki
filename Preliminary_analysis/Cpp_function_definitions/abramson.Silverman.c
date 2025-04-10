#include <R.h>
#include <math.h>

void AbramsonPixels( int *nx, double *xmax, int *ny, double *ymax, 
                     int *npoints, double *px, double *py, double *pcoef,
                     double *ans, double *h, double *alpha )
{

/* 
   Image of nx by ny pixels representing [0, xmax] x [0, ymax] 
   It is written to a matrix, nrow=ny, ncol=nx, by row. Hence
   for i in 1, ..., nx and j = 1, ... ny
      the (j-1)*nx + i th member of ans corresponds to matrix[j,i]
      or in Euclidean terms to
         ( ( i - 0.5 ) * xmax / nx, ( j - 0.5 ) * ymax / ny )
   In C, the first element is labelled 0

   pcoef contains density estimators for the points

*/

   int nxCoord = *nx;
   int nyCoord = *ny;
   int nIm = nxCoord * nyCoord;
   double maxx = *xmax;
   double maxy = *ymax;
   double xstep = maxx / nxCoord;
   double ystep = maxy / nyCoord;

/*
  Point pattern of length npoints, coordinates (px, py) 
  pcoef contains pilot values of intensity
*/
  
   int nPoints = *npoints;
   double h2 = *h; 
   h2 = h2 * h2;
   double abram = *alpha;

   int k, i, j;
   double xi, yj;
   double dist2, kappa, bw2;

   double geomMean = 0.0;
   for( k=0; k<nPoints; k++ ) geomMean += log( pcoef[k] );
   geomMean = exp( geomMean / nPoints);

   double edge[nPoints];
   for( k=0; k<nPoints; k++ ) edge[k] = 0.0;

   for( k=0; k<nPoints; k++ ) {
      bw2 = exp( log(pcoef[k] /geomMean) * 2 * abram );

      for( i=0; i<nxCoord; i++ ) { 
         xi = (i+0.5) * xstep; 
 
         for( j=0; j<nyCoord; j++ ) {
            yj = (j+0.5) * ystep;

            dist2 = ( px[k] - xi ) * ( px[k] - xi ) + ( py[k] - yj ) * ( py[k] - yj ); 
            kappa = exp( - dist2 / ( 2 * h2 * bw2 ) );

            edge[k] += kappa;
            if( isnan( edge[k] ) ) { 
              fprintf(stderr, " k = %d  bw2 = %f\n", k, bw2 ); 
              exit(1); 
            }
         }
      }

      edge[k] = ( ( xstep * ystep ) * edge[k] ) / ( bw2 * h2 * 2 * M_PI );
   }


/* Initialise */

   for( k=0; k < nIm; k++ ) ans[k] = 0.0;

/* Fill */

   for( i=0; i<nxCoord; i++ ) {
      xi = (i+0.5) * xstep;

      for(j=0; j<nyCoord; j++ ) {
         yj = (j+0.5) * ystep;

         for( k=0; k<nPoints; k++ ) {
            dist2 = ( px[k] - xi ) * (px[k] - xi ) + ( py[k] - yj ) * ( py[k] - yj );
            bw2 = exp(  log(pcoef[k]/geomMean) * 2 * abram );
            kappa = exp( - dist2  / ( 2 * h2 * bw2 ) ); 

            ans[j*nxCoord + i] += ( kappa / ( bw2 * edge[k] ) );
            if( isnan( ans[j*nxCoord + i] ) ) { 
              fprintf(stderr, " i = %d, j = %d, bw2 = %f\n", i, j, bw2 ); 
              exit(1);
            }
         }
      }
   }  


   for( k=0; k<nIm; k++ ) ans[k] = ans[k] / ( h2  * 2 * M_PI );

}


void AbramsonPoints( int *npoints, double *px, double *py, double *pcoef,
                     double *ans, double *h, double *alpha )
{

/*
  Point pattern of length npoints, coordinates (px, py) 
  pcoef contains pilot values of intensity
*/
  
   int nPoints = *npoints;
   double h2 = *h; 
   h2 = h2 * h2;
   double abram = *alpha;

   int k, i;
   double xi, yi;
   double dist2, kappa, bw2;

   double geomMean = 0.0;
   for( k=0; k<nPoints; k++ ) geomMean += log( pcoef[k] );
   geomMean = exp( geomMean / nPoints);

/* Initialise */

   for( k=0; k < nPoints; k++ ) ans[k] = 0.0;

/* Fill */

   for( i=0; i<nPoints; i++ ) {
      xi = px[i];
      yi = py[i];

         for( k=0; k<nPoints; k++ ) {
            dist2 = ( px[k] - xi ) * (px[k] - xi ) + ( py[k] - yi ) * ( py[k] - yi );
            bw2 = exp(  log(pcoef[k] / geomMean ) * 2 * abram );
            kappa = exp( - dist2  / ( 2 * h2 * bw2 ) ); 

            ans[i] += ( kappa / bw2 );
            if( isnan( ans[i] ) ) { 
              fprintf(stderr, " i = %d  bw2 = %f\n", i, bw2 ); 
              exit(1); 
            }
         }
   }  


   for( k=0; k<nPoints; k++ ) ans[k] = ans[k] / ( h2 * 2 * M_PI );

}



