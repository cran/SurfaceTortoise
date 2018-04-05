#' @title SurfaceTortoise

#' @description Optimizing spatial sampling using the Surface Tortoise
#' algoritm. Grid sampling and random stratified sampling are also available.

#' @author Kristin Piikki & Mats Söderström ,  \email{kristin.piikki@@slu.se}

#' @param x1 Raster layer. Additionally four raster layers (x2, x3, x4 and x5)
#' can be used. These will be transferred to the same extent and resulution as x1.
#' All rasters will be centered and scaled. The first principal component of the
#' rasters will be used as the surface for wich to optimize the sampling.
#' @param x2 Raster layer, see x1.
#' @param x3 Raster layer, see x1.
#' @param x4 Raster layer, see x1.
#' @param x5 Raster layer, see x1.
#' @param y SpatialPolygonsDataframe delineating the area to be sampled.
#' The coordinate system shall be the same as for the rasters. If not completely
#'  overlapping, the itersection of x1 and y will be sampled.
#' @param epsg Epsg code (numeric) for the orthognonal spatial reference system onto
#' which x1 and y are projected, see e.g. http://spatialreference.org/ref/epsg/
#' @param out_folder Output folder (path) to which the output files shall be
#'  exported.
#' @param out_prefix Prefix (character) for the output filenames.
#' @param method Sampling alogorithm: 'dir' = directed (SurfaceTortoise
#' algorithm), 'stratdir' = stratified directed (SurfaceTortoise algorithm),
#' 'grid' = regular grid and 'stratrand' = random stratified.
#' @param ncell_meanfilter Optional. Side of the square window (number of
#' raster cells) used for mean filtering of the raster(s).
#' @param p_idw Power exponent used for idw-interpolation (SurfaceTortoise
#' algorithm)
#' @param nmax_idw Number of neighbouring samples used for idw-interpolation
#'  (SurfaceTortoise algorithm).
#' @param edge Buffer zone (metre) inside the sample area border, where sampling is
#' prohibited.
#' @param strat_size Cell side (metre) of a square staratification grid.
#' @param min_dist Minimum distance allowed between samples. Valid for the
#' 'dir' and the 'stratdir' methods.
#' @param stop_n A stopping criterium. No more samples will be added when this
#' number of samples has been reached. Valid for the 'dir',
#' 'stratdir' and 'stratrand' methods.
#' @param stop_dens1 A stopping criterium. No more samples will be added when
#' this number of samples per hectare has been reached. Valid
#' for 'dir', stratdir' and 'stratrand' methods.
#' @param stop_dens2 A stopping criterium. No more samples will be added when
#' this number of samples per stratum has been reached. Valid
#' for 'stratdir' and 'stratrand' methods.
#' @param stop_change A stopping criterium. No more samples will be added when
#' the improvement in mean absolute error obtained by the last added sample is
#' less than this percentage of the hitherto improvement since sample no 3.
#' Valid for the 'dir' and the 'stratdir' methods.
#' @param plot_results Logical. shall results be plotted.

#' @return  A list with 1) a dataframe of sample coordinates 2) a spatialPointsDataFrame
#' with sample locations and 3) a SpatialPolygonsDataFrame with the
#' stratification grid or sampling area.If an output folder is specified, the following files will be exported to the specified output folder:
#' A point shapefile with sample locations; a tab-separated text file with
#' sample coordinates; and a polygon shapefile with the stratification grid.

#' @details The Surface Tortoise algorithm (ST-algorithm) uses a
#' spatial covariate in raster format (or the first principal component of a maximum
#' of five covariates in raster format) to find optimal locations for sampling.
#' The sampling strategy is based on the principle that an interpolation of the
#' samples should results be a similar surface as the covariate. When sample
#' locations are assigned, first the raster cell with the maximum deviation from
#' the covariate raster mean is sampled. Then the raster cell with the maximum deviation
#' from the first sampled raster cell is sampled. From then on, the values of the
#' sampled raster cells are interpolated by inverse distance weighting (idw)
#' and the raster cell with the largest absolute difference to the covariate
#' (error) is sampled. A new idw interpolation is made and a new cell is sampled.
#' This is repeated until any of the specified stopping criteria is reached.
#' The sampling can be stratified, which means that the areas is split by a
#' square grid. When a sample has been located in a stratum, no more samples
#' can be placed in that stratum again until all other strata have been sampled.
#' The likelihood for a clipped stratum, e.g. at the edge of the area to be sampled,
#' is equal to the area of that stratum divided by the area of a full stratum.

#' @import raster
#' @import gstat
#' @import rgeos
#' @import sp
#' @import gtools
#' @import rgdal

#' @importFrom grDevices terrain.colors
#' @importFrom graphics par
#' @importFrom stats as.formula
#' @importFrom stats complete.cases
#' @importFrom stats prcomp
#' @importFrom stats runif
#' @importFrom utils write.table
#' @importFrom utils globalVariables

#' @references Olsson, D. 2002. A method to optimize soil sampling from
#' ancillary data. Poster presenterad at: NJF seminar no. 336,
#' Implementation of Precision Farming in Practical Agriculture, 10-12 June
#' 2002, Skara, Sweden.

#' @examples
#' data(boundary)
#' grid.sampling<-tortoise(y=boundary,method='grid',edge=30,strat_size=50,
#' min_dist=10,plot_results=FALSE)

#' @export
tortoise<-function(x1=NA,
  x2=NA,
  x3=NA,
  x4=NA,
  x5=NA,
  y=NA,
  epsg=3006,
  out_folder=NA,
  out_prefix='st_',
  method='stratdir',
  ncell_meanfilter=3,
  p_idw=2,
  nmax_idw=8,
  edge=0,
  strat_size=100,
  min_dist=0,
  stop_n=100,
  stop_dens1=100,
  stop_dens2=100,
  stop_change=0.0001,
  plot_results=T
  )
  {
  #import data
  list.a<-import(r_name1=x1, r_name2=x2, r_name3=x3, r_name4=x4, r_name5=x5,
    b_name=y,cc=epsg)

  #prepare data
  list.a<-prepare(m=list.a,method=method, mean.filter=ncell_meanfilter, buff= edge)

  #create stratification grid
  list.a<-stratify(m=list.a,  s=strat_size, method=method, stop.n=stop_n)

  #do stratified or unstratified directed sampling
  if(method=='dir'|method=='stratdir') {
    list.a<-dir_sample(m=list.a, p.idw=p_idw, nmax.idw=nmax_idw,
    mindist=min_dist, stop.n=stop_n,stop.dens1=stop_dens1, stop.dens2=stop_dens2,
    stop.change=stop_change, plot.results=plot_results)
  }

  #do grid sampling
  if(method=='grid') list.a<-grid_sample(m=list.a, plot.results=plot_results)

  #do random stratified sampling
  if(method=='stratrand') {
    list.a<-rand_sample(m=list.a, mindist=min_dist, plot.results=plot_results, stop.dens2=stop_dens2)
  }

  #export data
  if(!is.na(out_folder)) export(m=list.a, output.folder=out_folder, method=method, output.prefix=out_prefix)

  #return
  return(list.a)
  }
