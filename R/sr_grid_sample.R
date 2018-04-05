#' @title Regular Grid Sample Placement

#' @description Creates a regular grid sampling design

#' @import raster
#' @import gstat
#' @import rgeos
#' @import sp
#' @import gtools
#' @import rgdal

#' @importFrom utils globalVariables

#' @param m a list of all objects to be used.
#' @param plot.results logical. shall results be plotted. default=TRUE.

#fix no visible binding for global variables
if(getRversion() >= "2.15.1")  globalVariables(c('strat2.sp', 'b3.sp'))

#define function
grid_sample<-function(m, plot.results){
  
#unlist objects
list.a<-m
for (h in 1:length(list.a)){
  n<-names(list.a)[h]
  g<-unlist(list.a[h], use.names=F, recursive = T)
  if(is.list(g))g<-g[[1]]
  if(!is.null(n))assign(x=n, value= g)
}

#set graphical parameters
par(mfrow=c(1,1))

#find centroids of strata
centroids.sp<-gCentroid(strat2.sp, byid=T) ##create sentroids to be used for grid sampling
crs(centroids.sp)<-crs(b3.sp)
centroids.sp<-crop(centroids.sp, b3.sp)
crs(centroids.sp)<-crs(b3.sp)

#convert to data.frame
p.df<-data.frame(centroids.sp@coords)
p.df$no<-1:nrow(p.df)

#convert to spatialpointsdata.frame
p.sp<-p.df
coordinates(p.sp)<-~x+y

#plot results
if(plot.results==T){
  plot(b.sp)
  plot(strat.sp, add=T)
  plot(p.sp, pch=16, add=T)
}

#give feedback
print('-I have completed your sampling!')

#return objects
rm('list.a')
return(list(mget(c('p.sp', 'strat.sp')))[[1]])
}
