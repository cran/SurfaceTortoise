#' @title Creates Stratification Grid

#' @description Creates stratification grid according to user specifications

#' @import raster
#' @import gstat
#' @import rgeos
#' @import sp
#' @import gtools
#' @import rgdal

#' @importFrom utils globalVariables

#' @param m a list of all objects to be used
#' @param s cell side (m) of a square staratification grid.
#' @param method sampling alogorithm: 'dir' = directed (sr algorithm),
#' 'stratdir' = stratified directed (sr algorithm), 'grid' = regular grid and
#' 'stratrand' = random stratified.

#fix no visible binding for global variables
if(getRversion() >= "2.15.1")  globalVariables(c('as', 'method', 'b2.sp', 'stop.n'))

#define function
stratify<-function(m, s, method, stop.n){

#unlist objects
list.a<-m
for (h in 1:length(list.a)){
  n<-names(list.a)[h]
  g<-unlist(list.a[h], use.names=F, recursive = T)
  if(is.list(g))g<-g[[1]]
  if(!is.null(n))assign(x=n, value= g)
}

#create stratification polygons
ee<-raster::extent(b3.sp)+c(-s, s, -s, s) ##buffer extent
strat.r<-raster(ext=ee,resolution=s) ##create raster
values(strat.r)<-1:ncell(strat.r)
strat.sp<-as(strat.r,'SpatialPolygonsDataFrame')
crs(strat.sp)<-crs(b3.sp)
strat2.sp<-strat.sp ##copy to be used for grid and stratrand sampling
strat.sp<-crop(strat.sp, b3.sp)

if(method=='dir'){
  strat.sp<-b2.sp #overwrites stratification with one polygon,if area is not to be stratified
  stop.dens2<-stop.n+1
}

if(method=='dir'|method=='stratdir'){
#convert stratification polygons to stratification raster
strat.r<-rasterize(strat.sp, r)

#specify the number of remaining samples to take in each stratum
no.r<-r
values(no.r)<-1
strat.sp$no<-as.vector(extract(x=no.r, y=strat.sp, fun=sum))
strat.sp$no<-strat.sp$no/max(strat.sp$no, na.rm=T)

#create diagniostic objects
refills<-0
refill.needed<-0
}

#give feedback
if(method!='dir'& method!='grid')print('-I have created a stratification grid!')
print('-Now I will place the samples, hang on!')

#return objects'
rm('list.a')
if(method=='dir'|method=='stratdir')return(list(mget(c('b.sp', 'b2.sp','b3.sp','r','strat.r', 'strat.sp', 'strat2.sp', 'area.r', 'refills', 'refill.needed')))[[1]])
if(method=='grid'|method=='stratrand')return(list(mget(c('b.sp', 'b2.sp','b3.sp','strat.sp','strat2.sp')))[[1]])
}

