library(EnvStats)
library(plyr)
library(ncdf4)
library(raster)
library(rasterVis)
library(maptools)
library(maps)
library(rgdal)
library(colorspace)
library(fields)
library(rgeos)

##########################
#this function extracts data from a raster for a grid
#layer is the raster that contains the data
#x is the number of polygon in the grid
#grid is the grid object
extract_data_grid<-function(layer,x,grid){
  cat('getting values for cell ',x,'\n')
  data<-as.data.frame(extract(layer,spTransform(grid[[x]],CRS=proj4string(layer))))
  data<-na.omit(data)
  data<-data[,1]
  return(data)
}

####present day elevation with ETOPO####
##downloaded from https://www.ngdc.noaa.gov/mgg/global/ (ETOPO1 Ice Surface)
##present, aggregating
get_present_elevation_grid_aggregate<-function(etopo1.path,grid){
  raster.ETOPO<-raster(etopo1.path)
  crs(raster.ETOPO)<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  #open land, merge and check geometry
  #land polygons including major islands downloaded from https://www.naturalearthdata.com/downloads/10m-physical-vectors/
  land<-readOGR('./raw_data/ne_10m_land/','ne_10m_land')
  #plot(land)
  land<-gUnaryUnion(land)
  land<-gBuffer(land, width=0)
  gIsValid(land)
  #plot(land)
  #mask with land raster to remove sea
  #this will get rid of negative elevations (sea values) which would affect the estimates of elevation in coastal cells
  raster.ETOPO.land<-mask(raster.ETOPO,land)
  #plot(raster.ETOPO.land)
  #aggregate at 1x1 (this speeds up calculations)
  raster.ETOPO.land.agg<-aggregate(raster.ETOPO.land,fact=60,fun=mean)
  system.time(grid.ETOPO.land.agg<-lapply(c(1:length(grid)),function(x)extract_data_grid(raster.ETOPO.land.agg,x,grid)))
  #get mean, range and tri within each cell
  grid.ETOPO.land.agg.elevation<-lapply(grid.ETOPO.land.agg,function(x)c(mean(x),range(x)[2]-range(x)[1]))
  #and output to table
  grid.ETOPO.land.agg.elevation.df<-as.data.frame(cbind(c(1:length(grid)),unlist(lapply(grid.ETOPO.land.agg.elevation,function(x)x[[1]])),unlist(lapply(grid.ETOPO.land.agg.elevation,function(x)x[[2]]))))
  colnames(grid.ETOPO.land.agg.elevation.df)<-c('cells','mean.elevation.ETOPO.land.agg','range.elevation.ETOPO.land.agg')
  return(grid.ETOPO.land.agg.elevation.df)
}

get_present_elevation_grid<-function(etopo1.path,grid){
  raster.ETOPO<-raster(etopo1.path)
  crs(raster.ETOPO)<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  #open land, merge and check geometry
  #land polygons including major islands downloaded from https://www.naturalearthdata.com/downloads/10m-physical-vectors/
  land<-readOGR('./raw_data/ne_10m_land/','ne_10m_land')
  #plot(land)
  land<-gUnaryUnion(land)
  land<-gBuffer(land, width=0)
  gIsValid(land)
  #plot(land)
  #mask with land raster to remove sea
  #this will get rid of negative elevations (sea values) which would affect the estimates of elevation in coastal cells
  raster.ETOPO.land<-mask(raster.ETOPO,land)
  #plot(raster.ETOPO.land)
  system.time(grid.ETOPO.land<-lapply(c(1:length(grid)),function(x)extract_data_grid(raster.ETOPO.land,x,grid)))
  #get mean, range and tri within each cell
  grid.ETOPO.land.elevation<-lapply(grid.ETOPO.land,function(x)c(mean(x),range(x)[2]-range(x)[1]))
  #and output to table
  grid.ETOPO.land.elevation.df<-as.data.frame(cbind(c(1:length(grid)),unlist(lapply(grid.ETOPO.land.elevation,function(x)x[[1]])),unlist(lapply(grid.ETOPO.land.elevation,function(x)x[[2]]))))
  colnames(grid.ETOPO.land.elevation.df)<-c('cells','mean.elevation.ETOPO.land','range.elevation.ETOPO.land')
  return(grid.ETOPO.land.elevation.df)
}

#PRISM4 elevation data downloaded from
#https://geology.er.usgs.gov/egpsc/prism/4_data.html
get_past_elevation_grid<-function(PRISM4.path,grid){
  raster.PRISM4<-raster(PRISM4.path)
  crs(raster.PRISM4)<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  #open land, merge and check geometry
  #land polygons including major islands downloaded from https://www.naturalearthdata.com/downloads/10m-physical-vectors/
  land.PRISM4<-raster('./raw_data/rasters/Plio_enh_LSM_v1.0.nc')
  plot(raster.PRISM4)
  plot(land.PRISM4)
  #mask the sea cells (value 0 in the topo raster)
  raster.PRISM4.land<-mask(raster.PRISM4,land.PRISM4,maskvalue=0)
  #extract data overlapping the grid and the raster
  #this takes ca 10 minutes
  system.time(grid.PRISM4<-lapply(c(1:length(grid)),function(x)extract_data_grid(raster.PRISM4.land,x,grid)))
  #get mean elevation in each cell
  grid.PRISM4.df<-as.data.frame(cbind(c(1:length(grid)),unlist(lapply(grid.PRISM4,function(x)mean(x))),unlist(lapply(grid.PRISM4,function(x)range(x)[2]-range(x)[1]))))
  colnames(grid.PRISM4.df)<-c('cells','mean.elevation.PRISM4.land','range.elevation.PRISM4.land')
  return(grid.PRISM4.df)
}
