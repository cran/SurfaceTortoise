#' @title Exports Data

#' @description Exports data to the specified output folder

#' @import raster
#' @import gstat
#' @import rgeos
#' @import sp
#' @import gtools
#' @import rgdal

#' @importFrom utils globalVariables
 
#' @param m a list of all objects to be used.
#' @param output.folder the output folder (path) to which the output files shall be
#' exported.
#' @param output.prefix a prefix (character) for the output filenames.
#' @param method sampling alogorithm: 'dir' = directed (sr algorithm),
#' 'stratdir' = stratified directed (sr algorithm), 'grid' = regular grid and
#' 'stratrand' = random stratified.

#fix no visible binding for global variables
if(getRversion() >= "2.15.1") globalVariables(c('p.sp', 'strat.sp', 'p.df'))

#define function
export<-function(m, method, output.folder, output.prefix){

#unlist objects
list.a<-m
for (h in 1:length(list.a)){
  n<-names(list.a)[h]
  g<-unlist(list.a[h], use.names=F, recursive = T)
  if(is.list(g))g<-g[[1]]
  if(!is.null(n))assign(x=n, value= g)
} 
  

if(dir.exists(output.folder)==F) dir.create(output.folder)

#export files
shapefile(
  p.sp,
  filename=paste0(output.folder,'\\',output.prefix, 'points.shp'),
  overwrite=T
)
shapefile(
  strat.sp,
  filename=paste0(output.folder,'\\',output.prefix, 'strata.shp'),
  overwrite=T
)
p.df<-data.frame(p.sp@coords, p.sp@data)
if(method=='dir'|method=='stratdir')p.df<-p.df[,c('x','y','no','mae')]
if(method=='grid'|method=='stratrand')p.df<-p.df[,c('x','y','no')]

write.table(
  x=p.df,
  file=paste0(output.folder,'\\',output.prefix, 'points.txt'),
  sep = "\t",
  row.names = F,
  col.names = T,
  quote=F
  )

#give feedback
print('-Please check your output folder. I have exported your results!')

#return objects
rm('list.a')
return(list(mget(c('p.df', 'p.sp', 'strat.sp')))[[1]])
}

