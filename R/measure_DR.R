library(ape)
library(picante)

#this generates a table with DR metric for all mammal species in the tree
measure_DR_tree_table<-function(treefile){
  if(class(treefile)=='phylo'){
    tree<-treefile
  }else{
    tree<-read.tree(treefile)  
  }
  DR.results<-evol.distinct(tree,type='equal.splits')
  colnames(DR.results)[2]<-c('equal.splits')
  DR.results$DR<-1/DR.results$equal.splits
  DR.results<-DR.results[c(1,3)]
  return(DR.results)
}

#this generates a table with DR metric for a list of dataframes with DR metric 
get_pseudoposterior_median_DRtable<-function(DR.pseudoposterior.file,path,name){
  DR.pseudoposterior.file<-readRDS(DR.pseudoposterior.file)
  DR.pseudoposterior.species<-as.character(DR.pseudoposterior.file[[1]]$Species)
  DR.pseudoposterior<-lapply(DR.pseudoposterior.file,function(x)x[,2])
  DR.pseudoposterior<-as.data.frame(cbind(DR.pseudoposterior.species,do.call('cbind',DR.pseudoposterior)),stringsAsFactors = F)
  DR.pseudoposterior.median<-as.numeric(apply(DR.pseudoposterior,1,function(x) median(as.numeric(x[c(2:length(x))]))))
  DR.pseudoposterior.median<-as.data.frame(cbind(DR.pseudoposterior.species,DR.pseudoposterior.median),stringsAsFactors = F)
  colnames(DR.pseudoposterior.median)<-c('Species','median.pseudoposterior.DR')
  write.table(DR.pseudoposterior.median,file=paste(path,name,'_pseudoposteriorDRmedian.txt',sep=''),quote=F,sep='\t',row.names=F)
}

#this plots the DR with the median DR across the pseudoposterior
compare_DR_MCC_vs_medianpseudoposterior<-function(DR.tablefile,DR.pseudoposterior.tablefile){
  DR<-read.table(DR.tablefile,sep='\t',header=T,stringsAsFactors = F)
  DR.pseudoposterior<-read.table(DR.pseudoposterior.tablefile,sep='\t',header=T,stringsAsFactors = F)
  DR.all<-merge(DR,DR.pseudoposterior,by='Species')
  corr<-cor.test(DR.all$DR,DR.all$median.pseudoposterior.DR)
  plot(DR.all$DR,DR.all$median.pseudoposterior.DR,xlim=c(0,3),ylim=c(0,3),pch=16,cex=.5,main=paste('DR_MCC-DR pseudoposterior correlation (r=',round(corr$estimate,3),')',sep=''),cex.main=.8,col="#00000033",xlab='DR.MCC',ylab='median.DR.pseudoposterior')
  cat(paste('Pearsons r = ',round(corr$estimate,3),'\n',sep=''))
}

#this plots the histogram of correlations of DR measures in the MCC tree vs each tree in the 100 randoms
compare_DR_MCC_vs_allpseudoposterior<-function(DR.tablefile,DR.pseudoposterior.file){
  DR<-read.table(DR.tablefile,sep='\t',header=T,stringsAsFactors = F)
  DR.pseudoposterior<-readRDS(DR.pseudoposterior.file)
  DR.correlations.pseudoposterior<-lapply(DR.pseudoposterior,function(x){DR.all<-merge(DR,x,by='Species');return(cor.test(DR.all$DR.x,DR.all$DR.y)$estimate)})
  hist(unlist(DR.correlations.pseudoposterior),xlim=c(0,1),breaks=50,yaxs='i',xaxs='i',main='',xlab='Pearsons r correlation',ylab='')
}

DR_stats_grid<-function(list.species.ranges.DR,speciesDR){
  #for species
  #calculate meanDR per cell
  list.species.ranges.DR<-lapply(names(list.species.ranges),function(x){species.DR<-rep(speciesDR[speciesDR$Species==x,'DR'],times=length(list.species.ranges[[x]]));if(length(species.DR)==0){return(NA)};names(species.DR)<-list.species.ranges[[x]];species.DR})
  names(list.species.ranges.DR)<-names(list.species.ranges)
  #delete species with NA
  list.species.ranges.DR<-list.species.ranges.DR[!is.na(list.species.ranges.DR)]
  list.species.ranges.DR.flat<-list.species.ranges.DR
  names(list.species.ranges.DR.flat)<-NULL
  list.species.ranges.DR.flat<-unlist(list.species.ranges.DR.flat)
  names(list.species.ranges.DR.flat)<-as.numeric(names(list.species.ranges.DR.flat))
  #get mean of wDR per cell
  cell.mean.species.DR<-lapply(c(1:max(as.numeric(names(list.species.ranges.DR.flat)))),function(x){cat(x,'\n'); list.species.ranges.DR.flat[grep(paste('^',x,'$',sep=''),names(list.species.ranges.DR.flat))]})
  #calculate weightedDR per cell
  #for each species, get the DR estimate + divide it by each range + create a named vector with each cell number
  list.species.ranges.wDR<-lapply(names(list.species.ranges),function(x){species.wDR<-rep(speciesDR[speciesDR$Species==x,'DR'],times=length(list.species.ranges[[x]]));if(length(species.wDR)==0){return(NA)};names(species.wDR)<-list.species.ranges[[x]];species.wDR/length(species.wDR)})
  names(list.species.ranges.wDR)<-names(list.species.ranges)
  #delete species with NA
  list.species.ranges.wDR<-list.species.ranges.wDR[!is.na(list.species.ranges.wDR)]
  list.species.ranges.wDR.flat<-list.species.ranges.wDR
  names(list.species.ranges.wDR.flat)<-NULL
  list.species.ranges.wDR.flat<-unlist(list.species.ranges.wDR.flat)
  names(list.species.ranges.wDR.flat)<-as.numeric(names(list.species.ranges.wDR.flat))
  #get mean of wDR per cell
  cell.mean.species.wDR<-lapply(c(1:max(as.numeric(names(list.species.ranges.wDR.flat)))),function(x){cat(x,'\n'); list.species.ranges.wDR.flat[grep(paste('^',x,'$',sep=''),names(list.species.ranges.wDR.flat))]})
  species.cells.DR<-as.data.frame(cbind(c(1:16517),unlist(lapply(cell.mean.species.DR,function(x)mean(x,na.rm=T))),unlist(lapply(cell.mean.species.wDR,function(x)mean(x,na.rm=T))),unlist(lapply(cell.mean.species.wDR,function(x)geoMean(x,na.rm=T)))),stringsAsFactors=F)
  colnames(species.cells.DR)<-c('cells','mean.DR','mean.wDR','geomean.wDR')
  return(species.cells.DR)
}


