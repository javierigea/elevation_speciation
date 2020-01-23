library(rgdal)
library(gdalUtils)
library(parallel)
all.features<-ogrInfo(dsn = './raw_data/birds/BOTW/BOTW.gdb/',layer = 'All_Species','so')
nfeatures<-all.features$nrows

for (i in seq(from=1,to=nfeatures,by=500)){
  start<-i
  end<-i+499
  if(end>nfeatures){
    end<-nfeatures
  }
  cat(i,'\n')
  mclapply(mc.cores=2,c(start:end),function(x) ogr2ogr('./raw_data/birds/BOTW/BOTW.gdb/','All_Species',fid=x,append=TRUE,dst_datasource_name = paste('./raw_data/birds/birds_ranges_2016/',start,'/',sep=''),nln=paste('All_Species_',x,sep='')))
}




