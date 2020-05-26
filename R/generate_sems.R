
library(spdep)
library(ape)
library(piecewiseSEM)
library(spatialreg)
library(rgdal)

sem_spatial_pseudoposterior <- function(variables_table, replicate){
  #with wDR
  cells.table<-read.table(variables_table, sep='\t',header=T,stringsAsFactors = F)
  #replace DR metrics with pseudoposterior replicate
  mammals_replicateDR <- read.table(paste0('./output/mammals/tables/DR_posterior/mammals_DRpseudoposterior_', 
                                           replicate,
                                           '_cells_table.txt'), 
                                    sep='\t',
                                    header=T,
                                    stringsAsFactors = F)
  birds_replicateDR <- read.table(paste0('./output/birds/tables/DR_posterior/birds_DRpseudoposterior_',
                                         replicate,'_cells_table.txt'), 
                                  sep='\t',
                                  header=T,
                                  stringsAsFactors = F)
  
  cells.table$mammals.mean.DR <- mammals_replicateDR$mean.DR
  cells.table$mammals.mean.wDR <- mammals_replicateDR$mean.wDR
  cells.table$mammals.geomean.wDR <- mammals_replicateDR$geomean.wDR
  cells.table$birds.mean.DR <- birds_replicateDR$mean.DR
  cells.table$birds.mean.wDR <- birds_replicateDR$mean.wDR
  cells.table$birds.geomean.wDR <- birds_replicateDR$geomean.wDR
  
  #drop cells where speciation = 0 (there are no species)
  cells.table<-cells.table[!is.na(cells.table$mammals.mean.DR)&!is.na(cells.table$birds.mean.DR)&cells.table$mammals.mean.DR>0&cells.table$birds.mean.DR>0,]
  #drop cells where past or present elevation is NA or negative
  cells.table<-cells.table[!is.na(cells.table$mean.elevation.ETOPO.land)&!is.na(cells.table$mean.elevation.PRISM4.land)&cells.table$mean.elevation.ETOPO.land>0&cells.table$mean.elevation.PRISM4.land>0,]
  #drop unnecessary columns
  cells.table<-cells.table[,c('cells','mammals.mean.wDR','birds.mean.wDR','mean.elevation.ETOPO.land','mean.elevation.PRISM4.land','mean.present.T','mean.past.T','present.minus.past.elevation','present.minus.past.temperature')]
  cells.table<-cells.table[complete.cases(cells.table),]
  #log transform response variables
  cells.table$mammals.mean.wDR<-log(cells.table$mammals.mean.wDR)
  cells.table$birds.mean.wDR<-log(cells.table$birds.mean.wDR)
  #log transform other variables with positive values
  cells.table$mean.elevation.ETOPO.land<-log(cells.table$mean.elevation.ETOPO.land)
  cells.table$mean.elevation.PRISM4.land<-log(cells.table$mean.elevation.PRISM4.land)
  #scale predictors
  cells.table[,c(4:9)]<-apply(cells.table[,c(4:9)],2,function(x)scale(x))
  
  
  #load grid and create neighbourhood if neighbours.1000.w object does not exist
  if(exists("neighbours.1000.w") == F){
    grid.world<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
    #get coordinates
    grid.world.longlat<-lapply(grid.world,function(x) spTransform(x,CRS("+proj=longlat")))
    grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) sp::coordinates(x))
    grid.coordinates<-do.call("rbind", grid.coordinates)
    #get neighbours in 1000km distqnce
    neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
    neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)
    
  }
  #sem with sarlm for mammals wDR, elevation, temperature and change elevation and change temperature
  sarlm.sem.mammals.wDR.elevation.temp<-psem(errorsarlm(mammals.mean.wDR~mean.elevation.ETOPO.land+present.minus.past.elevation+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
  sarlm.sem.birds.wDR.elevation.temp<-psem(errorsarlm(birds.mean.wDR~mean.elevation.ETOPO.land+present.minus.past.elevation+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
  
  obj_sarlmsem <- list(sarlm.sem.mammals.wDR.elevation.temp, sarlm.sem.birds.wDR.elevation.temp)
  names(obj_sarlmsem) <- c('mammals', 'birds')
  saveRDS(obj_sarlmsem, 
          file = paste0('./output/world/sems/pseudoposterior/global/sarlm_sems_wDR_elevation_temp_global_replicate_', 
                        replicate,
                        '.rds'))
  
}


sem_gain_spatial_pseudoposterior <- function(variables_table, replicate){
  #with wDR
  cells.table<-read.table(variables_table, sep='\t',header=T,stringsAsFactors = F)
  #replace DR metrics with pseudoposterior replicate
  mammals_replicateDR <- read.table(paste0('./output/mammals/tables/DR_posterior/mammals_DRpseudoposterior_', 
                                           replicate,
                                           '_cells_table.txt'), 
                                    sep='\t',
                                    header=T,
                                    stringsAsFactors = F)
  birds_replicateDR <- read.table(paste0('./output/birds/tables/DR_posterior/birds_DRpseudoposterior_',
                                         replicate,'_cells_table.txt'), 
                                  sep='\t',
                                  header=T,
                                  stringsAsFactors = F)
  
  cells.table$mammals.mean.DR <- mammals_replicateDR$mean.DR
  cells.table$mammals.mean.wDR <- mammals_replicateDR$mean.wDR
  cells.table$mammals.geomean.wDR <- mammals_replicateDR$geomean.wDR
  cells.table$birds.mean.DR <- birds_replicateDR$mean.DR
  cells.table$birds.mean.wDR <- birds_replicateDR$mean.wDR
  cells.table$birds.geomean.wDR <- birds_replicateDR$geomean.wDR
  
  #drop cells where speciation = 0 (there are no species)
  cells.table<-cells.table[!is.na(cells.table$mammals.mean.DR)&!is.na(cells.table$birds.mean.DR)&cells.table$mammals.mean.DR>0&cells.table$birds.mean.DR>0,]
  #drop cells where past or present elevation is NA or negative
  cells.table<-cells.table[!is.na(cells.table$mean.elevation.ETOPO.land)&!is.na(cells.table$mean.elevation.PRISM4.land)&cells.table$mean.elevation.ETOPO.land>0&cells.table$mean.elevation.PRISM4.land>0,]
  #drop unnecessary columns
  cells.table<-cells.table[,c('cells','mammals.mean.wDR','birds.mean.wDR','mean.elevation.ETOPO.land','mean.present.T','elevation.gain','present.minus.past.temperature')]
  #get columns where elevation.gain >0
  cells.table<-cells.table[cells.table$elevation.gain>0,]
  cells.table<-cells.table[complete.cases(cells.table),]
  #log transform response variables
  cells.table$mammals.mean.wDR<-log(cells.table$mammals.mean.wDR)
  cells.table$birds.mean.wDR<-log(cells.table$birds.mean.wDR)
  #log transform other variables with positive values
  cells.table$mean.elevation.ETOPO.land<-log(cells.table$mean.elevation.ETOPO.land)
  cells.table$elevation.gain<-log(cells.table$elevation.gain)
  #scale predictors
  cells.table[,c(4:7)]<-apply(cells.table[,c(4:7)],2,function(x)scale(x))
  
  
  #load grid and create neighbourhood if neighbours.1000.w object does not exist
  if(exists("neighbours.1000.w") == F){
    grid.world<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
    #get coordinates
    grid.world.longlat<-lapply(grid.world,function(x) spTransform(x,CRS("+proj=longlat")))
    grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) sp::coordinates(x))
    grid.coordinates<-do.call("rbind", grid.coordinates)
    #get neighbours in 1000km distqnce
    neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
    neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)
    
  }
  #sem with sarlm for mammals wDR, elevation, temperature and change elevation and change temperature
  sarlm.sem.mammals.wDR.elevation.gain.temp<-psem(errorsarlm(mammals.mean.wDR~mean.elevation.ETOPO.land+elevation.gain+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
  sarlm.sem.birds.wDR.elevation.gain.temp<-psem(errorsarlm(birds.mean.wDR~mean.elevation.ETOPO.land+elevation.gain+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
  
  obj_sarlmsem <- list(sarlm.sem.mammals.wDR.elevation.gain.temp, sarlm.sem.birds.wDR.elevation.gain.temp)
  names(obj_sarlmsem) <- c('mammals', 'birds')
  saveRDS(obj_sarlmsem, 
          file = paste0('./output/world/sems/pseudoposterior/gain_elevation/sarlm_sems_wDR_gainelevation_temp_replicate_', 
                        replicate,
                        '.rds'))
  
}



sem_loss_spatial_pseudoposterior <- function(variables_table, replicate){
  #with wDR
  cells.table<-read.table(variables_table, sep='\t',header=T,stringsAsFactors = F)
  #replace DR metrics with pseudoposterior replicate
  mammals_replicateDR <- read.table(paste0('./output/mammals/tables/DR_posterior/mammals_DRpseudoposterior_', 
                                           replicate,
                                           '_cells_table.txt'), 
                                    sep='\t',
                                    header=T,
                                    stringsAsFactors = F)
  birds_replicateDR <- read.table(paste0('./output/birds/tables/DR_posterior/birds_DRpseudoposterior_',
                                         replicate,'_cells_table.txt'), 
                                  sep='\t',
                                  header=T,
                                  stringsAsFactors = F)
  
  cells.table$mammals.mean.DR <- mammals_replicateDR$mean.DR
  cells.table$mammals.mean.wDR <- mammals_replicateDR$mean.wDR
  cells.table$mammals.geomean.wDR <- mammals_replicateDR$geomean.wDR
  cells.table$birds.mean.DR <- birds_replicateDR$mean.DR
  cells.table$birds.mean.wDR <- birds_replicateDR$mean.wDR
  cells.table$birds.geomean.wDR <- birds_replicateDR$geomean.wDR
  
  #drop cells where speciation = 0 (there are no species)
  cells.table<-cells.table[!is.na(cells.table$mammals.mean.DR)&!is.na(cells.table$birds.mean.DR)&cells.table$mammals.mean.DR>0&cells.table$birds.mean.DR>0,]
  #drop cells where past or present elevation is NA or negative
  cells.table<-cells.table[!is.na(cells.table$mean.elevation.ETOPO.land)&!is.na(cells.table$mean.elevation.PRISM4.land)&cells.table$mean.elevation.ETOPO.land>0&cells.table$mean.elevation.PRISM4.land>0,]
  #drop unnecessary columns
  cells.table<-cells.table[,c('cells','mammals.mean.wDR','birds.mean.wDR','mean.elevation.ETOPO.land','mean.present.T','elevation.loss','present.minus.past.temperature')]
  #get columns where elevation.loss >0
  cells.table<-cells.table[cells.table$elevation.loss>0,]
  cells.table<-cells.table[complete.cases(cells.table),]
  #log transform response variables
  cells.table$mammals.mean.wDR<-log(cells.table$mammals.mean.wDR)
  cells.table$birds.mean.wDR<-log(cells.table$birds.mean.wDR)
  #log transform other variables with positive values
  cells.table$mean.elevation.ETOPO.land<-log(cells.table$mean.elevation.ETOPO.land)
  cells.table$elevation.loss<-log(cells.table$elevation.loss)
  #scale predictors
  cells.table[,c(4:7)]<-apply(cells.table[,c(4:7)],2,function(x)scale(x))
  
  
  #load grid and create neighbourhood if neighbours.1000.w object does not exist
  if(exists("neighbours.1000.w") == F){
    grid.world<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
    #get coordinates
    grid.world.longlat<-lapply(grid.world,function(x) spTransform(x,CRS("+proj=longlat")))
    grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) sp::coordinates(x))
    grid.coordinates<-do.call("rbind", grid.coordinates)
    #get neighbours in 1000km distqnce
    neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
    neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)
    
  }
  #sem with sarlm for mammals wDR, elevation, temperature and change elevation and change temperature
  sarlm.sem.mammals.wDR.elevation.loss.temp<-psem(errorsarlm(mammals.mean.wDR~mean.elevation.ETOPO.land+elevation.loss+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
  sarlm.sem.birds.wDR.elevation.loss.temp<-psem(errorsarlm(birds.mean.wDR~mean.elevation.ETOPO.land+elevation.loss+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
  
  obj_sarlmsem <- list(sarlm.sem.mammals.wDR.elevation.loss.temp, sarlm.sem.birds.wDR.elevation.loss.temp)
  names(obj_sarlmsem) <- c('mammals', 'birds')
  saveRDS(obj_sarlmsem, 
          file = paste0('./output/world/sems/pseudoposterior/loss_elevation/sarlm_sems_wDR_losselevation_temp_replicate_', 
                        replicate,
                        '.rds'))
  
}