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
##########################
substract_current_minus_past_bio1<-function(path.present,path.past){
  bio1_current<-raster(path.present)
  bio1_past<-raster(path.past)
  #compare values of extent of rasters
  xmin<-max(c(extent(bio1_current)@xmin,extent(bio1_past)@xmin))
  xmax<-min(c(extent(bio1_current)@xmax,extent(bio1_past)@xmax))
  ymin<-max(c(extent(bio1_current)@ymin,extent(bio1_past)@ymin))
  ymax<-min(c(extent(bio1_current)@ymax,extent(bio1_past)@ymax))
  extent.both<-extent(c(xmin,xmax,ymin,ymax))
  #crop to minimal common extent
  bio1_current_crop<-crop(bio1_current,extent.both,snap="near")
  bio1_past_crop<-crop(bio1_past,extent.both,snap="near")
  #resample because origin is different
  bio1_past_crop<-resample(bio1_past_crop,bio1_current_crop)
  #stack rasters
  bio1_both<-stack(bio1_current_crop,bio1_past_crop)
  #substract current - past to get change in bio1
  bio1_change<-overlay(bio1_both[[1]],bio1_both[[2]],fun=function(x,y){x-y})
  #divide by ten to get degrees centigrades
  bio1_change<-bio1_change/10
  return(bio1_change)
}


get_present_temperature_grid<-function(bio1raster.path,grid){
  bio1.raster<-raster(bio1raster.path)
  system.time(grid.currentT<-lapply(c(1:length(grid)),function(x)extract_data_grid(bio1_current_crop,x,grid)))
  grid.currentT.values<-lapply(grid.currentT,function(x)c(mean(x),range(x)[2]-range(x)[1]))
  grid.currentT.values.df<-as.data.frame(cbind(c(1:length(grid)),unlist(lapply(grid.currentT.values,function(x)x[[1]])),unlist(lapply(grid.currentT.values,function(x)x[[2]]))))
  colnames(grid.currentT.values.df)<-c('cells','mean.present.T','range.present.T')
  grid.currentT.values.df$range.present.T[!is.finite(grid.currentT.values.df$range.present.T)]<-0
  grid.currentT.values.df$mean.present.T[is.na(grid.currentT.values.df$mean.present.T)]<-NA
  return(grid.currentT.values.df)
}  

get_past_temperature_grid<-function(bio1raster.past.path,grid){
  bio1.past.raster<-raster(bio1raster.past.path)
  system.time(grid.pastT<-lapply(c(1:length(grid)),function(x)extract_data_grid(bio1.past.raster,x,grid)))
  grid.pastT.values<-lapply(grid.pastT,function(x)c(mean(x),range(x)[2]-range(x)[1]))
  grid.pastT.values.df<-as.data.frame(cbind(c(1:length(grid)),unlist(lapply(grid.pastT.values,function(x)x[[1]])),unlist(lapply(grid.pastT.values,function(x)x[[2]]))))
  colnames(grid.pastT.values.df)<-c('cells','mean.past.T','range.past.T')
  grid.pastT.values.df$range.past.T[!is.finite(grid.pastT.values.df$range.past.T)]<-0
  grid.pastT.values.df$mean.past.T[is.na(grid.pastT.values.df$mean.past.T)]<-NA
  return(grid.pastT.values.df)
}  
  
