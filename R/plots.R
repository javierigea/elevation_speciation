library(colorspace)
library(fields)

####################
#this function creates a vector of diverge_hsvcolours from an initial vector of values 
create_colours_vector<-function(ncategories,vector){
  colours.vector<-diverge_hsv(n=ncategories)
  quantiles.colour<-quantile(vector, seq(from=0,to=1,by=(1/ncategories)),na.rm=TRUE)
  grid.colours<-vector()
  #assign colours to each element of vector 
  for (i in 1:length(vector)){
    if(is.na(vector[i])){
      grid.colours[i]<-'grey'
      next;
    }
    new.vector<-sort(c(quantiles.colour,vector[i]))
    interval<-which(new.vector==vector[i])-1
    if(length(interval)>1){
      interval<-interval[1]
    }
    if(interval==0){
      interval<-1
    }
    grid.colours[i]<-colours.vector[interval]
  }
  return(grid.colours)
}
####################





####################
#this function plots the global grid with cells coloured according to variable
#ncategories is the number of colours in the gradient
#positive.values is a logical vector; use only rows where variable >0

plot_grid_worldmap_variable<-function(table.env.file,variable,ncategories,positive.values){
  #open table
  table.env<-read.table(table.env.file,header=T,sep='\t',stringsAsFactors = F)
  #check presence of variable in colnames; exit if failed
  if(!variable%in%colnames(table.env)){
    stop('variable not present in table, check name','\n')
  }
  #subset table to variable + cells
  table.run<-table.env[,c('cells',variable)]
  table.run<-table.run[complete.cases(table.run),]
  #check if negative values have to be removed
  if(positive.values==T){
    table.run<-table.run[which(table.run[,2]>0),]  
  }
  colours.variable<-create_colours_vector(ncategories=ncategories,vector=table.run[,2])
  #load grid
  grid.world<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
  #this assumes that latlon is in kms
  plot(c(1,1),xlim=c(-20050,20050),ylim=c(-6342,6325),type='n',xaxt='n',yaxt='n',xlab='',ylab='',main=variable)
  lapply(grid.world,function(x) plot(x,add=T,lwd=0.05))
  lapply(c(1:length(table.run[,2])),function(x) plot(grid.world[[table.run[x,'cells']]],add=T,lwd=0.05,col=colours.variable[x]))
  quantiles.variable<-quantile(table.run[,2],seq(from=0,to=1,by=(1/10)),na.rm=TRUE)
  legend.text<-sapply(c(1:10),function(x)paste(round(quantiles.variable,5)[x],round(quantiles.variable,5)[x+1],sep='-'))
  legend(-15000,0,legend=rev(legend.text),fill=rev(diverge_hsv(n=10)),cex=.4,bty='n')
}
plot_grid_worldmap_variable_scalebar<-function(table.env.file,variable,ncategories,positive.values){
  #open table
  table.env<-read.table(table.env.file,header=T,sep='\t',stringsAsFactors = F)
  #check presence of variable in colnames; exit if failed
  if(!variable%in%colnames(table.env)){
    stop('variable not present in table, check name','\n')
  }
  #subset table to variable + cells
  table.run<-table.env[,c('cells',variable)]
  table.run<-table.run[complete.cases(table.run),]
  #check if negative values have to be removed
  if(positive.values==T){
    table.run<-table.run[which(table.run[,2]>0),]  
  }
  colours.variable<-create_colours_vector(ncategories=ncategories,vector=table.run[,2])
  #load grid
  grid.world<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
  #this assumes that latlon is in kms
  plot(c(1,1),xlim=c(-20050,20050),ylim=c(-6342,6325),type='n',xaxt='n',yaxt='n',xlab='',ylab='',main=variable)
  lapply(grid.world,function(x) plot(x,add=T,lwd=0.05))
  lapply(c(1:length(table.run[,2])),function(x) plot(grid.world[[table.run[x,'cells']]],add=T,lwd=0.05,col=colours.variable[x]))
  quantiles.variable<-quantile(table.run[,2],seq(from=0,to=1,by=(1/10)),na.rm=TRUE)
  #legend.text<-sapply(c(1:10),function(x)paste(round(quantiles.variable,5)[x],round(quantiles.variable,5)[x+1],sep='-'))
  #legend(-15000,0,legend=rev(legend.text),fill=rev(diverge_hsv(n=10)),cex=.4,bty='n')
  colorbar.plot(-15000,-5000,as.matrix(c(1:ncategories)),col=diverge_hsv(n=ncategories),horizontal=FALSE,strip.width=1/50)
  legend(-19000,-4000,max(table.run[,2],na.rm=T),bty='n',cex=.8)
  legend(-19000,-5000,min(table.run[,2],na.rm=T),bty='n',cex=.8)
  legend(-19000,-4500,median(table.run[,2],na.rm=T),bty='n',cex=.8)

}


