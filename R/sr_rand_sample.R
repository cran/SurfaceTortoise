#' @title Random Stratified Sample Placement

#' @description Creates a random stratified sampling design

#' @import raster
#' @import gstat
#' @import rgeos
#' @import sp
#' @import gtools
#' @import rgdal

#' @importFrom utils globalVariables

#' @param m a list of all objects to be used.
#' @param mindist minimum distance allowed between samples. default=0 m.
#' @param stop.dens2 a stopping criterium. No more samples will be added when
#' this number of samples per stratum has been reached. default = 100. valid
#' for 'stratdir' and 'stratrand' methods.
#' @param plot.results logical. shall results be plotted. default=TRUE.

#fix no visible binding for global variables
if(getRversion() >= "2.15.1") globalVariables(c('strat2.sp', 'stop.dens2', 'b3.sp'))

#define function
rand_sample<-function(m, mindist, plot.results, stop.dens2){

#unlist objects
list.a<-m
for (h in 1:length(list.a)){
    n<-names(list.a)[h]
    g<-unlist(list.a[h], use.names=F, recursive=T)
    if(is.list(g))g<-g[[1]]
    if(!is.null(n))assign(x=n, value= g)
}


#set graphical parameters
par(mfrow=c(1,1))

#select strata to sample(i.e. strata with >50 percent of the area within the buffered area boundary)
strat3.sp<-crop(strat2.sp, b3.sp)
strat3.sp$no<-1:nrow(strat3.sp)
strat3.sp$nsample<-round(stop.dens2*area(strat3.sp)/max(area(strat3.sp)))

#create empty data.frame to fill out
p.df<-as.data.frame(matrix(nrow=sum(strat3.sp$nsample), ncol=4, data=NA))
colnames(p.df)<-c('stratum', 'x','y','no')
for (k in unique((strat3.sp$nsample))){
  nnk<-rep(as.vector(strat3.sp$no[strat3.sp$nsample==k]), k)
  if(exists('nn')) nn<-c(nn, nnk) else  nn<- nnk
}
p.df$stratum<-nn
p.df$no<-1:nrow(p.df)

#add one random data point remove a buffer zone around that point and add the next random point
for (i in 1:nrow(p.df)){
  #identify row idex for the stratum
  a<-p.df[i,'stratum']
  ##crop selected polygon by buffered samples in order to avoid sampling too close to other samples.
  if(i==1)cc<-strat3.sp[a,]
  if(i>1)cc<-erase(strat3.sp[a,], buff.sp)
  ##make random sampling of x points
  pi.sp<-spsample(cc, n=1, type='random', byid=T) ##create centroids to be used for grid sampling
  ##add coordinatesto regristry df
  p.df[i,c('x','y')]<-pi.sp@coords
  #convert to spatialpointsdata.frame
  p.sp<-p.df[complete.cases(p.df),]
  coordinates(p.sp)<-~x+y
  crs(p.sp)<-crs(strat2.sp)
  #buffer sampling points
  buff.sp<-buffer(x=p.sp, width=mindist ,dissolve=T)
  ##plot
  if(plot.results==T){
    plot(b.sp)
    plot(strat3.sp, add=T)
    plot(p.sp, pch=16,add=T)
  }

}
if(plot.results==T){
  p.sp<-crop(p.sp, b3.sp)
  plot(strat2.sp, border='white')
  plot(b.sp, add=T, col='white')
  plot(strat.sp, add=T, col='white')
  plot(p.sp, pch=16, add=T)
}

#give feedback
print('-I have completed your sampling!')

#return objects'
rm('list.a')
return(list(mget(c('p.sp', 'strat.sp')))[[1]])
}
