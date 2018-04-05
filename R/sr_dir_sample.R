#' @title Directed Sample Placement

#' @description Creates a spaital sampling design by the Surface Reconstructor algorithm

#' @import raster
#' @import gstat
#' @import rgeos
#' @import sp
#' @import gtools
#' @import rgdal

#' @importFrom utils globalVariables

#' @param m a list of all objects to be used.
#' @param p.idw the power exponent used for idw-interpolation (sr algorithm)
#' @param nmax.idw the number of neighbouring samples used for idw-interpolation
#'  (sr algorithm)
#' @param mindist minimum distance allowed between samples. default=0 m.
#' @param stop.n a stopping criterium. No more samples will be added when this
#' number of samples has been reached. default = 100. valid for the 'dir',
#' 'stratdir' and 'stratrand' methods.
#' @param stop.dens1 a stopping criterium. No more samples will be added when
#' this number of samples per hectare has been reached. default = 100. valid
#' for all methods.
#' @param stop.dens2 a stopping criterium. No more samples will be added when
#' this number of samples per stratum has been reached. default = 100. valid
#' for 'stratdir' and 'stratrand' methods.
#' @param stop.change a stopping criterium. No more samples will be added when
#' the improvement in mean absolute error obtained by the last added sample is
#' less than this percentage of the hitherto improvement since sample no 3.
#' default = 0.0001. valid for the 'dir' and the 'stratdir'.
#' @param plot.results logical. shall results be plotted. default=TRUE.

#fix no visible binding for global variables
if(getRversion() >= "2.15.1") globalVariables(c('r', 'area.r', 'b.sp'))

#define function
dir_sample<-function(m, p.idw, nmax.idw, mindist, stop.n,stop.dens1, stop.dens2,
  stop.change, plot.results)
{
  
#unlist objects
list.a<-m
for (h in 1:length(list.a)){
  n<-names(list.a)[h]
  g<-unlist(list.a[h], use.names=F, recursive = T)
  if(is.list(g))g<-g[[1]]
  if(!is.null(n))assign(x=n, value= g)
}

#set graphical parameters
par(mfrow=c(2,2))

#create empty data.frame to fill out
p.df<-as.data.frame(matrix(nrow=stop.n,ncol=7,data=NA))
names(p.df)<-c('x','y','r.value','mae','no','density','change')

#run sampling loop
for (k in 1:stop.n){

  #if needed, refill number of samples left to take in strata
  if(k>1)if(refill.needed==1){
    strat.sp$no<-strat.sp$no+1
    refills<-refills+1
    refill.needed<-0
  }

  ##calculate prediction raster (idw.r)
  if(k==1) {
    r.mean<-mean(area.r@data@values, na.rm=T)
    idw.r <- r
    idw.r[is.na(idw.r)==F]<-r.mean
  }

  ##calculate aboslute error raster (ae.r)
  ae.r<-abs(idw.r-r)

  #restrict sampling area
  m1.sp<-strat.sp[round(strat.sp$no)>=1,] ## mask away strata with no samples left to take
  if(nrow(m1.sp)>0)ae2.r<-mask(ae.r, m1.sp)
    if(k==1)ae3.r<-ae2.r ## mask away areas too close to other samples
  if(k>1)ae3.r<-mask(ae2.r, erase(b.sp,buffer(p.sp, mindist)))

  #convert masked raster to a data.frame
  ae.df<-data.frame(coordinates(ae3.r),values(area.r), values(ae3.r))
  names(ae.df)<-c('x','y','r.value','ae')

  #identify cell with maximum absolute error (ae), or for the first ssample the minimum ae
  if(k==1)pi.df<-ae.df[which.min(ae.df$ae),]
  if(k!=1)pi.df<-ae.df[which.max(ae.df$ae),]
  pi.df<-pi.df[round(runif(1, 1, nrow(pi.df))),] ## in case several cells has the same ae, on of these cells are randomly chosen.

  #enter data into the register data.frame (p.df)
  p.df[k,1:4]<-pi.df[1,1:4]
  p.df[k,'no']<-k
  p.df[k,'density']<-k*10000/sum(area(b.sp))

  #convert samples to SpatialPointsDataFrame
  p.sp<-p.df[complete.cases(p.df[c('x','y')]),]
  coordinates(p.sp)<-~x+y
  crs(p.sp)<-crs(area.r)

  #adjust number of samples left to target in strata
  pi.sp<-p.sp[nrow(p.sp),]
  ff<-over(pi.sp, strat.sp)
  strat.sp[strat.sp$layer==ff$layer,'no']<-max(0,ff$no-1)

  #again, mask away strata with no samples left to take and sample buffers, and check if refill is needed
  m1.sp<-strat.sp[round(strat.sp$no)>=1,]
  if(nrow(m1.sp)>0){
  ae4.r<-mask(ae3.r,m1.sp)
  m2.sp<-erase(b.sp,buffer(p.sp, mindist))
  ae5.r<-mask(ae4.r, m2.sp)
  if(sum(values(ae5.r), na.rm=T)==0) refill.needed<-1
  } else refill.needed<-1

  #compute new map
  if(k>1){
  idw.g <- gstat(id = 'r.value', formula = r.value~1, locations = ~x+y, data=p.df[complete.cases(p.df[,c('x','y','r.value')]),], nmax=nmax.idw, set=list(idp=p.idw))
  idw.r <- interpolate(r, idw.g, debug.level = 0)
  idw.r <- mask(idw.r, r)
  }

  #compute mean absolute error (mae) of new map
  mae<-mean(values(abs(idw.r-r)), na.rm=T)
  p.df[k,'mae']<-mae
  p.sp[k,'mae']<-mae

  #compute mae inprovement of new map
  current.improvement<-(p.df[k,'mae']-p.df[3,'mae'])
  previous.improvement<-(p.df[(k-1),'mae']-p.df[3,'mae'])
  change.improvement<-abs((current.improvement-previous.improvement)/previous.improvement)
  if(k>3) p.df[k,'change']<-100*change.improvement

  #set plot parameters
  if(plot.results==T & k==1){
    p1<-ceiling(mae)
    r.brk<-pretty(values(r),na.rm=T)
    ae.brk<-pretty(values(ae.r),na.rm=T)
  }

  #plot results
  if(plot.results==T){
    plot(r, col=rev(terrain.colors(length(r.brk)-1, alpha = 1)), breaks = r.brk, legend=F, main='Surface (green=high, white=low)', xaxt="n", yaxt="n")
    plot(strat.sp,add=T)
    plot(p.sp, add=T,pch = 16)
    plot(idw.r, col=rev(terrain.colors(length(r.brk)-1, alpha = 1)), breaks = r.brk, legend=F, main='Reconstructed surface', xaxt="n", yaxt="n")
    plot(ae.r, col=rev(terrain.colors(length(ae.brk)-1, alpha = 1)), breaks = ae.brk, legend=F, main='Mean absolute error (MAE)', xaxt="n", yaxt="n")
    plot(x=p.df$no,y=p.df$mae, ylim=c(0,p1), xlab = NA, ylab=NA, main='MAE vs no of samples', xaxt="n", yaxt="n")
    }

  #break loop when any of the specified stopping criteria is met
  if (k>1){
    if (p.df[k,'no']>=stop.n){print('The specified number of samples has been reached!'); break}
    if (p.df[k,'density']>=stop.dens1){print('-The specified number of samples per hectare has been reached!'); break}
    if (k>3)if(abs(p.df[k,'change'])<=stop.change){print('-The improvement by adding more samples is now smaller than the specified limit'); break}
    if (refills==(stop.dens2)){print('-The specified number of samples per stratum has been reached!'); break}
  }
}

#give feedback
print('-I have completed your sampling!')

#return objects
rm('list.a')
return(list(mget(c('p.sp', 'strat.sp')))[[1]])
}

