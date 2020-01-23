
#########function to add offset to cells in a table
add_offset_to_cells<-function(table,offset){
  table<-read.table(table,header=T,sep='\t',stringsAsFactors = F)
  for (i in 1:nrow(table)){
    if(table[i,'range.cells']==0){
      next
    }else{
      cells<-table[i,'cells']
      cells<-as.numeric(unlist(strsplit(cells,' ')))
      cells<-cells+offset
      cells<-paste(cells,collapse=' ')
      table[i,'cells']<-cells
    }
    
  }
  return(table)
}
####

###########function to add offset to hotspot cells in a table
add_offset_to_hotspot<-function(table,offset){
  table<-read.table(table,header=T,sep='\t',stringsAsFactors = F)
  for (i in 1:nrow(table)){
    #cat(i,'\n')
    if(table[i,'hotspots']==0){
      next
    }else{
      hotspots<-table[i,'hotspots']
      if(class(hotspots)=='integer'){
        hotspots<-hotspots+offset
        table[i,'hotspots']<-as.character(hotspots)
        next;
      }
      hotspots<-as.numeric(unlist(strsplit(hotspots,',')))
      if (length(hotspots)==1){
        hotspots<-hotspots+offset
        table[i,'hotspots']<-as.character(hotspots)
      }else{
        hotspots<-hotspots+offset
        hotspots<-paste(hotspots,collapse=',')
        table[i,'hotspots']<-hotspots
      }
      
    }
    
  }
  return(table)
}

#####wrapper function 1: to merge species grid occurrence and richness table for world
#path is the path where *_species_gridoccurrence and *_richness_grid_table' are found
merge_realms_speciesdata_into_world<-function(path){
  species.grids<-list.files(path=path,pattern='*_species_gridoccurrence')
  ####delete all_realms_species just in case
  if(length(grep('all_realms',species.grids))>0){
    species.grids<-species.grids[-grep('all_realms',species.grids)]  
  }
  ####
  cell.list<-list.files(path=path,pattern='*_richness_grid_table')
  if(length(grep('all_realms',cell.list))>0){
    cell.list<-cell.list[-grep('all_realms',cell.list)]
  }
  
  
  ####merging species occurrences
  cell.list.counts<-sapply(cell.list,function(x)system(paste('wc -l ',path,'/',x,sep=''),intern=TRUE))
  cell.list.counts<-sub(cell.list.counts,pattern='./output/.+',replacement = '')
  cell.list.counts<-sub(cell.list.counts,pattern=' +',replacement=' ')
  cell.list.counts<-sub(cell.list.counts,pattern=' ',replacement='')
  cell.list.counts<-sub(cell.list.counts,pattern=' ',replacement='')
  cell.list.counts<-as.numeric(cell.list.counts)
  #calculate offset
  offset.list<-cell.list.counts-1
  offset.list[2]<-offset.list[1]+offset.list[2]
  offset.list[3]<-offset.list[2]+offset.list[3]
  offset.list[4]<-offset.list[3]+offset.list[4]
  offset.list[5]<-offset.list[4]+offset.list[5]
  offset.list[6]<-offset.list[5]+offset.list[6]
  new.tables<-list()
  new.tables[[1]]<-read.table(paste(path,'/',species.grids[1],sep=''),header=T,sep='\t',stringsAsFactors = F)
  colnames(new.tables[[1]])[c(2,3)]<-sapply(colnames(new.tables[[1]])[c(2,3)],function(x)paste(x,'.1',sep=''))
  #add offset to realms
  for (i in 2:length(species.grids)){
    new.tables[[i]]<-add_offset_to_cells(table=paste(path,'/',species.grids[i],sep=''),offset=offset.list[i-1])
    colnames(new.tables[[i]])[c(2,3)]<-sapply(colnames(new.tables[[i]])[c(2,3)],function(x)paste(x,'.',i,sep=''))
  }
  merged.table<-Reduce(function(...) merge(..., by='spp',all=T), new.tables)
  colnames(merged.table)[2:ncol(merged.table)]<-c('cells.realm1','range.cells1','cells.realm2','range.cells2','cells.realm3','range.cells3','cells.realm4','range.cells4','cells.realm5','range.cells5','cells.realm6','range.cells6','cells.realm7','range.cells7')
  merged.table$sum.range.cells<-merged.table$range.cells1+merged.table$range.cells2+merged.table$range.cells3+merged.table$range.cells4+merged.table$range.cells5+merged.table$range.cells6+merged.table$range.cells7
  merged.table$all.realm.cells<-paste(merged.table$cells.realm1,merged.table$cells.realm2,merged.table$cells.realm3,merged.table$cells.realm4,merged.table$cells.realm5,merged.table$cells.realm6,merged.table$cells.realm7,sep=' ')
  world.table<-merged.table[,c('spp','all.realm.cells','sum.range.cells')]
  colnames(world.table)<-c('spp','cells','range.cells')
  write.table(world.table,file=paste(path,'/100_all_realms_species_gridoccurrence_table.txt',sep=''),sep='\t',quote=F,row.names=F)
  ###merging the grid species richness
  new.tables.richness<-list()
  new.tables.richness[[1]]<-read.table(paste(path,cell.list[1],sep=''),header=T,sep='\t',stringsAsFactors = F)
  
  for (i in 2:length(cell.list)){
    new.tables.richness[[i]]<-read.table(paste(path,cell.list[i],sep=''),header=T,sep='\t',stringsAsFactors = F)
    new.tables.richness[[i]]$cell<-new.tables.richness[[i]]$cell+offset.list[i-1]
  }
  world.richness.table <- do.call("rbind", new.tables.richness)
  write.table(world.richness.table,file=paste(path,'/100_all_realms_richness_grid_table.txt',sep=''),sep='\t',quote=F,row.names=F)
}

#####wrapper function 1: to merge species grid occurrence and richness table for world
#path is the path where *_hotspot_grid_table are found
merge_hotspotrealms_into_world<-function(path){
  ###merging hotspot grids
  hotspot.list<-list.files(path=path,pattern='*_hotspot_grid_table')
  hotspot.number<-numeric()
  for (i in 1:length(hotspot.list)){
    table<-read.table(paste(path,'/',hotspot.list[i],sep=''),header=T,sep='\t',stringsAsFactors = F)
    hotspot.codes.sorted<-sort(as.numeric(unique(unlist(strsplit(paste(names(table(table$hotspots)),collapse=','),',')))))
    hotspot.number[i]<-hotspot.codes.sorted[length(hotspot.codes.sorted)]
  }
  
  offset.list.hotspot<-hotspot.number
  offset.list.hotspot[2]<-offset.list.hotspot[1]+offset.list.hotspot[2]
  offset.list.hotspot[3]<-offset.list.hotspot[2]+offset.list.hotspot[3]
  offset.list.hotspot[4]<-offset.list.hotspot[3]+offset.list.hotspot[4]
  offset.list.hotspot[5]<-offset.list.hotspot[4]+offset.list.hotspot[5]
  offset.list.hotspot[6]<-offset.list.hotspot[5]+offset.list.hotspot[6]
  
  new.tables.hotspot<-list()
  new.tables.hotspot[[1]]<-read.table(paste(path,'/',hotspot.list[1],sep=''),header=T,sep='\t',stringsAsFactors = F)
  
  
  for (i in 2:length(hotspot.list)){
    new.tables.hotspot[[i]]<-add_offset_to_hotspot(table=paste(path,'/',hotspot.list[i],sep=''),offset=offset.list.hotspot[i-1])
  }
  
  world.hotspots.table <- do.call("rbind", new.tables.hotspot)
  world.hotspots.table$cell<-seq(from=1,to=nrow(world.hotspots.table),by=1)
  write.table(world.hotspots.table,file=paste(path,'/100_all_realms_hotspot_grid_table.txt',sep=''),sep='\t',quote=F,row.names=F)
}

#####wrapper function 3: to merge species grids into world grid
#path is the path where grid_[A-Z].+_100.rds are found
merge_realms_grid_into_world<-function(path){
  ###this merges the realm-grids into a world grid
  species.grids<-list.files(path=path,pattern='grid_[A-Z].+_100.rds')
  #grid.list<-lapply(species.grids[1:3],function(x)load(paste('./output/mammals/',x,sep='')))
  grid.list<-list()
  for (i in 1:length(species.grids)){
    grid<-readRDS(paste(path,'/',species.grids[i],sep=''))
    grid.list[[i]]<-grid
  }
  grid.world<-unlist(grid.list)
  saveRDS(grid.world,file=paste(path,'/grid_World_RealmsMerged_100.rds',sep=''))
}









