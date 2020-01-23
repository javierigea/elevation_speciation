library(ape)
library(BAMMtools)
library(coda)
library(mvtnorm)
library(lattice)
library(plyr)
####this function calculates the sampling fraction and runs setBAMMpriors for a treefile
prepare_BAMM_input<-function(treefile,dictionaryfile,path,name){
  tree<-read.tree(treefile)
  table.synonyms<-read.csv(dictionaryfile,header=T)
  table.synonyms<-data.frame(table.synonyms$Genus,table.synonyms$Species,paste(table.synonyms$Genus,table.synonyms$Species,sep=' '),table.synonyms$Synonyms,stringsAsFactors = FALSE)
  colnames(table.synonyms)<-c('Genus','Species','Binomial','Synonyms')
  #clean extra characters
  table.synonyms$Synonyms<-sub(table.synonyms$Synonyms,pattern='  ',replacement='')
  table.synonyms<-table.synonyms[,c('Genus','Binomial')]
  table.synonyms<-unique(table.synonyms)
  #sampling fraction as number of tip labels / total entries in taxonomy
  sampling.fraction<-length(tree$tip.label)/nrow(table.synonyms)
  write.table(sampling.fraction,file=paste(path,name,'_BAMMsamplingfraction.txt',sep=''),quote=F,sep='\t',row.names=F)
  #checks
  cat('ultrametric check: ',is.ultrametric(tree),'\n')
  cat('binary check: ',is.binary.tree(tree),'\n')
  # Now to check min branch length:
  cat('min branch length: ',min(tree$edge.length),'\n')
  setBAMMpriors(phy=tree,total.taxa=nrow(table.synonyms),outfile = paste(path,name,'_BAMMpriors.txt',sep=''))
}

#mcmcmout:mcmcoutfile
analyse_BAMM_convergence<-function(mcmcout,burnin){
  mcmcout.div <- read.csv(mcmcout, header=T)
  #see if convergence is ok
  plot(mcmcout.div$logLik ~ mcmcout.div$generation)
  #remove initial 10% (burnin)
  burnstart.div <- floor(burnin * nrow(mcmcout.div))
  postburn.div <- mcmcout.div[burnstart.div:nrow(mcmcout.div), ]
  plot(postburn.div$logLik ~ postburn.div$generation)
  plot(postburn.div$logLik ~ postburn.div$N_shifts)
  #check ESS (should be >200)
  print(paste('ESS Lik',effectiveSize(postburn.div$logLik),sep=': '))
  print(paste('ESS Nshifts',effectiveSize(postburn.div$N_shifts),sep=': '))
  #compute posterior probabilities of models
  post_probs.div <- table(postburn.div$N_shifts) / nrow(postburn.div)
  names(post_probs.div)
  #plot the posterior distribution of Nrateshifts
  plot(names(post_probs.div),post_probs.div)
}

get_tipRates_BAMM<-function(treefile,eventfile,burnin,path,name){
  tree<-read.tree(treefile)
  cat('getting BAMMeventdata','\n')
  eventdata<-getEventData(phy=tree,eventdata = eventfile,burnin=burnin,nsamples = 1000,verbose = F,type='diversification')
  tipRates<-getTipRates(eventdata,returnNetDiv = F)
  tipRates.div<-getTipRates(eventdata,returnNetDiv = T)
  tipRates.df<-as.data.frame(cbind(names(tipRates$lambda.avg),tipRates$lambda.avg,tipRates$mu.avg,tipRates.div$netdiv.avg),stringsAsFactors = F)
  colnames(tipRates.df)<-c('spp','lambda.avg','mu.avg','netdiv.avg')
  row.names(tipRates.df)<-NULL
  write.table(tipRates.df,file=paste(path,name,'_BAMMTipRates.txt',sep=''),quote=F,row.names=F,sep='\t')
}

correlation_speciesDR_tipBAMM<-function(DRtable.file,BAMMtable.file){
  DRtable<-read.table(DRtable.file,header=T,sep='\t',stringsAsFactors = F)
  colnames(DRtable)[1]<-c('spp')
  BAMMtable<-read.table(BAMMtable.file,header=T,sep='\t',stringsAsFactors = F)
  colnames(BAMMtable)[1]<-c('spp')
  DRBAMMrates<-merge(DRtable,BAMMtable,by='spp')
  corr<-cor.test(log10(DRBAMMrates$DR),log10(DRBAMMrates$netdiv.avg))
  plot(log10(DRBAMMrates$DR),log10(DRBAMMrates$netdiv.avg),xlab='log10.DR',ylab='log10.BAMM.netdiv',pch=16,cex=.5,xaxt='n',yaxt='n',main=paste('BAMM-DR correlation (r=',round(corr$estimate,3),')',sep=''),cex.main=.8,col="#00000033")
  #smoothScatter(log10(DRBAMMrates$DR),log10(DRBAMMrates$netdiv.avg),xlab='log10.DR',ylab='log10.BAMM.netdiv',pch=16,cex=.5,xaxt='n',yaxt='n',main=paste('BAMM-DR correlation (r=',round(corr$estimate,3),')',sep=''),cex.main=.8)
  #plot(log10(DRBAMMrates$DR),log10(DRBAMMrates$netdiv.avg),xlab='log10.DR',ylab='log10.BAMM.netdiv',pch=16,cex=.5,main=paste('BAMM-DR correlation (r=',round(corr$estimate,3),')',sep=''),cex.main=.8)
  axis(1,at=log10(c(0.01,0.05,0.2,0.5,2)),labels=c(0.01,0.05,0.2,0.5,2))
  axis(2,at=log10(c(0.01,0.05,0.2)),labels=c(0.01,0.05,0.2))
}

BAMM_stats_grid<-function(list.species.ranges.BAMM,speciesBAMM){
  #for species
  #calculate meanBAMM per cell
  list.species.ranges.BAMM<-lapply(names(list.species.ranges),function(x){species.BAMM<-rep(speciesBAMM[speciesBAMM$Species==x,'lambda.avg'],times=length(list.species.ranges[[x]]));if(length(species.BAMM)==0){return(NA)};names(species.BAMM)<-list.species.ranges[[x]];species.BAMM})
  names(list.species.ranges.BAMM)<-names(list.species.ranges)
  #delete species with NA
  list.species.ranges.BAMM<-list.species.ranges.BAMM[!is.na(list.species.ranges.BAMM)]
  list.species.ranges.BAMM.flat<-list.species.ranges.BAMM
  names(list.species.ranges.BAMM.flat)<-NULL
  list.species.ranges.BAMM.flat<-unlist(list.species.ranges.BAMM.flat)
  names(list.species.ranges.BAMM.flat)<-as.numeric(names(list.species.ranges.BAMM.flat))
  #get mean of wBAMM per cell
  cell.mean.species.BAMM<-lapply(c(1:max(as.numeric(names(list.species.ranges.BAMM.flat)))),function(x){cat(x,'\n'); list.species.ranges.BAMM.flat[grep(paste('^',x,'$',sep=''),names(list.species.ranges.BAMM.flat))]})
  #calculate weightedBAMM per cell
  #for each species, get the BAMM estimate + divide it by each range + create a named vector with each cell number
  list.species.ranges.wBAMM<-lapply(names(list.species.ranges),function(x){species.wBAMM<-rep(speciesBAMM[speciesBAMM$Species==x,'lambda.avg'],times=length(list.species.ranges[[x]]));if(length(species.wBAMM)==0){return(NA)};names(species.wBAMM)<-list.species.ranges[[x]];species.wBAMM/length(species.wBAMM)})
  names(list.species.ranges.wBAMM)<-names(list.species.ranges)
  #delete species with NA
  list.species.ranges.wBAMM<-list.species.ranges.wBAMM[!is.na(list.species.ranges.wBAMM)]
  list.species.ranges.wBAMM.flat<-list.species.ranges.wBAMM
  names(list.species.ranges.wBAMM.flat)<-NULL
  list.species.ranges.wBAMM.flat<-unlist(list.species.ranges.wBAMM.flat)
  names(list.species.ranges.wBAMM.flat)<-as.numeric(names(list.species.ranges.wBAMM.flat))
  #get mean of wBAMM per cell
  cell.mean.species.wBAMM<-lapply(c(1:max(as.numeric(names(list.species.ranges.wBAMM.flat)))),function(x){cat(x,'\n'); list.species.ranges.wBAMM.flat[grep(paste('^',x,'$',sep=''),names(list.species.ranges.wBAMM.flat))]})
  species.cells.BAMM<-as.data.frame(cbind(c(1:16517),unlist(lapply(cell.mean.species.BAMM,function(x)mean(x,na.rm=T))),unlist(lapply(cell.mean.species.wBAMM,function(x)mean(x,na.rm=T))),unlist(lapply(cell.mean.species.wBAMM,function(x)geoMean(x,na.rm=T)))),stringsAsFactors=F)
  colnames(species.cells.BAMM)<-c('cells','mean.BAMM','mean.wBAMM','geomean.wBAMM')
  return(species.cells.BAMM)
}

