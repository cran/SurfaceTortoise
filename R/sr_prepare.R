#' @title Prepares Data

#' @description Prepares data to be used for sampling

#' @import raster
#' @import gstat
#' @import rgeos
#' @import sp
#' @import gtools
#' @import rgdal
 
#' @param m a list of all objects to be used.
#' @param mean.filter optional. side of the square window (number of
#' raster cells) used for mean filtering of the raster(s). default = 3
#' @param buff buffer zone (m) inside the sample area border, where sampling is
#' prohibited. default = 0.
#' @param method sampling alogorithm: 'dir' = directed (sr algorithm),
#' 'stratdir' = stratified directed (sr algorithm), 'grid' = regular grid and
#' 'stratrand' = random stratified.

prepare<-function(m, method, mean.filter, buff){

#unlist objects
list.a<-m
for (h in 1:length(list.a)){
  n<-names(list.a)[h]
  g<-unlist(list.a[h], use.names=F, recursive = T)
  if(is.list(g))g<-g[[1]]
  if(!is.null(n))assign(x=n, value=g)
}

#create boundary for area to sample
##buffer imported boundary polygon
b2.sp<-buffer(b.sp,-buff) #b2.sp is the buffered bundary polygon.Samples will be targeted within the intersect of this polygon and the imported raster(s).
b3.sp<-aggregate(b2.sp,dissolve=T)

#mean filter raster data
if(method=='dir'|method=='stratdir'){
aa<-is.na(r[[1]])
if(even(mean.filter))mean.filter<-mean.filter+1
if(!is.na(mean.filter)){
  for(j in 1:nlayers(r)){
    rj<-r[[j]]
    filter.m<-matrix(1,mean.filter,mean.filter)
    rj<-raster::focal(rj,w=filter.m,fun=mean,na.rm=T,pad=T)
    if(nlayers(r)>1)r[[j]]<-rj
    if(nlayers(r)==1)r<-rj
  }
}
r[aa]<-NA

#Derive first principal component raster(pc1 raster) from raster stack
if(nlayers(r)>1){
  pca_data<-as.data.frame(r)
  pca_data<-pca_data[complete.cases(pca_data),]
  bb<-sample(x=1:nrow(pca_data),size=min(5000,nrow(pca_data)))
  pca_data<-pca_data[bb,]
  frml<-as.formula(paste0('~',paste0(names(r),collapse='+')))
  mod<-prcomp(formula=frml,data=pca_data,subset=1:ncell(r),na.action='na.omit',scale=T,center=T)
  r$pc1<-predict(mod,as.data.frame(r))[,1]
  r<-r$pc1
}
##crop raster
area.r<-mask(r,b3.sp)
}

#give feedback
print('-and now I have prepared your data!')

#return objects
rm('list.a')
if(method=='dir'|method=='stratdir')return(list(mget(c('r', 'b.sp', 'b2.sp', 'b3.sp', 'area.r')))[[1]])
if(method=='grid'|method=='stratrand')return(list(mget(c('b.sp', 'b2.sp', 'b3.sp')))[[1]])
}
