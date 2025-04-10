
is.inside <- function( z = NULL, box = NULL ) {
  
  if( missing(z) || !is.numeric(z) ) stop("time not given or not numeric")
  
  inside <- FALSE
  if( box$zrange[1] <= z & z<box$zrange[2] )  inside <- TRUE
  
  return( list( z=z, box=box, inside=inside ) )
}

#box is list( zrange ) with zrange = c(a,b)


#Construct boxes: it returns all dummy points and all boxes around them as a list
# list_centers: dummy points
# list_box: box around dummy points
# indexd: 1

dummy <- function( box=NULL, tfrac=NULL )
{
  if( missing(tfrac) || !is.numeric(tfrac) ) { stop("tfrac missing or non-numeric") }
  
  c <- box$zrange[2] - box$zrange[1]
  tfrac <- round(tfrac)
  
  high <- data.frame()
  for( k in 1:tfrac ) high <- rbind( high, c( box$zrange[1] + (c/tfrac)*k ) )
  
  low <- data.frame()
  for( k in 1:tfrac ) low <- rbind( low, c( box$zrange[1] + (c/tfrac)*(k-1) ) )
  
  list_box <- list()
  for( i in 1:dim(high)[1]) list_box[[i]] <- list( zrange=c( low[i,1], high[i,1] ) )
  
  list_centers <- list()
  for( i in 1:length(list_box) ) 
    list_centers[[i]] <- c( ( list_box[[i]]$zrange[1] + list_box[[i]]$zrange[2] ) / 2 )
  
  dummy_list <- list() 
  for( i in 1:length(list_box) ) 
    dummy_list[[i]] <- list( dum=list_centers[[i]], box=list_box[[i]], indexd=1 )
  
  return( dummy_points=dummy_list )
  
}


volume.box <- function( box ) 
{
  v <- (box$zrange[2]-box$zrange[1] )
  return( list( volume=v, box=box ) )
}


pseudo_likelihood <- function( X, boxTime, tfrac ) {
  # X a list of times on (boxTime[1], boxTime[2]), 
  # tfrac governs dummy points
  
  if( is.null(X) ) { stop("X not given: needs to be a data frame") }
  if( !is.data.frame(X) ) { stop("X needs to be a data frame, time the first component") }
  if( missing(tfrac) || !is.numeric(tfrac) ) { 
    stop("Fraction tfrac not given or not numeric") }
  
  
  dummy_list <- dummy( box=boxTime, tfrac = tfrac )
  data_list <- list()
  for( i in 1:dim(X)[1] ) data_list[[i]] <- list( point= as.numeric(X[i,1]), indexp=1 )
  
  # Box counts for dummy and data points
  
  for( i in 1:length(dummy_list) ) {
    for( j in 1:length(data_list) ) { 
      if( is.inside( data_list[[j]]$point, box=dummy_list[[i]]$box )$inside==TRUE ) {
        dummy_list[[i]]$indexd <- dummy_list[[i]]$indexd + 1
      }
    }
  }
  
  # Update index for data
  
  for( j in 1:length(data_list) ) {
    for( i in 1:length(dummy_list) ){ 
      if( is.inside( data_list[[j]]$point, box=dummy_list[[i]]$box )$inside==TRUE ) {
        data_list[[j]]$indexp <- dummy_list[[i]]$indexd
      }
    }
  }
  
  if( sum( sapply( X=dummy_list, FUN=function(x) { sum(x$indexd) } ) ) 
      != ( length(data_list) + length(dummy_list) ) ) {
    stop("The sum of indexes is not equal to dummy+data") }
  
  # combine dummy with data points: two columns, (nr dummy + nr data) rows 
  
  u_points <- data.frame( cbind( 
    c( sapply( dummy_list, function(x) { x$dum } ),
       sapply( data_list, function(x) { x$point } ) ) ) , 
    c( sapply( dummy_list, function(x) { x$indexd } ), 
       sapply( data_list, function(x) { x$indexp } ) ) )
  names( u_points ) <- c("time", "index")
  
  #quadrature weights
  area_box <- volume.box( dummy_list[[1]]$box)$volume 
  
  u_points <- cbind(u_points, w = area_box / u_points$index )
  
  #indicators
  u_points <- cbind( u_points, zz=c(rep(0, length(dummy_list)), 
                                    rep(1, length(data_list)) ) )
  
  u_points <- cbind( u_points, yy = u_points$zz / u_points$w )
  
  return( list( u_points = u_points, n_dummy = length(dummy_list), 
                n_data = length(data_list), volume=area_box ) )
}
