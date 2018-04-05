#' @title Imports Data

#' @description imports specified datasets

#' @import raster
#' @import gstat
#' @import rgeos
#' @import sp
#' @import gtools
#' @import rgdal

#' @importFrom utils globalVariables
#' @importFrom grDevices is.raster

#' @param m a list of all objects to be used.
#' @param r_name1 a raster layer. additionally four raster layers (x2, x3, x4 and x5)
#' can be used. These will be transferred to the same extent and resulution as x1.
#' All rasters will be centered and scaled. The first principal component of the
#' rasters will be used as the surface for wich to optimize the sampling.
#' @param r_name2 a raster layer, see x1.
#' @param r_name3 a raster layer, see x1.
#' @param r_name4 a raster layer, see x1.
#' @param r_name5 a raster layer, see x1.
#' @param b_name a SpatialPolygonsDataframe delineating the area to be sampled.
#' The coordinate system shall be the same as for the rasters. If not completely
#' overlapping, the itersection of x1 and y will be sampled.
#' @param cc the epsg code (numeric) for the spatial reference system onto
#' which x1 and y are projected, see e.g. http://spatialreference.org/ref/epsg/

#fix no visible binding for global variables
if(getRversion() >= "2.15.1") globalVariables(c('list.a'))

#define function
import<-function(m, r_name1, r_name2, r_name3, r_name4, r_name5,b_name,cc){

#import shapefile
  if(is.character(b_name)==T)b.sp<-shapefile(b_name)
  if(is.character(b_name)==F)b.sp<-b_name

#import and stack rasters
if(is.character(r_name1)==T){
r_names<-c(r_name1,r_name2,r_name3,r_name4,r_name5)
for(i in 1:length(r_names)){
  if(!is.na(r_names[i])){
    ri<-raster(r_names[i])
    if(exists('r')){
      ri<-resample(x=ri,y=r)
      ri[is.na(r[[1]])]<-NA
      r[is.na(ri)]<-NA
      r<-stack(r,ri)
      }
    else r<-ri
    }
  }
}
if(is.character(r_name1)==F){
 r<-r_name1
 if(is.raster(r_name2)){r2<-resample(x=r_name2,y=r);  r2[is.na(r)]<-NA; r[is.na(r2)]<-NA; r<-stack(r,r2)}
 if(is.raster(r_name3)){r3<-resample(x=r_name3,y=r);  r3[is.na(r)]<-NA; r[is.na(r3)]<-NA; r<-stack(r,r3)}
 if(is.raster(r_name4)){r4<-resample(x=r_name4,y=r);  r4[is.na(r)]<-NA; r[is.na(r4)]<-NA; r<-stack(r,r4)}
 if(is.raster(r_name5)){r5<-resample(x=r_name5,y=r);  r5[is.na(r)]<-NA; r[is.na(r5)]<-NA; r<-stack(r,r5)}
}
print(1)
#define projection
if(exists('cc')&is.raster(r)){
crs(r)<-NA;proj4string(r)<-CRS(paste0("+init=epsg:",cc))
crs(b.sp)<-NA;proj4string(b.sp)<-CRS(paste0("+init=epsg:",cc))
}

#give feedback
print('Tortoise says:')
print('-I have imported your data!')

#return objects'
return(list(mget(c('r', 'b.sp')))[[1]])
}

