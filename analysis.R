#this is the analysis script for the elevation_speciation repo

#load libraries
library(hexbin)
library(spdep)
library(ape)
library(piecewiseSEM)
library(corrplot)
library(spatialreg)
library(rgdal)
####TO DO: delete local path here####
setwd('~/Dropbox/Work_in_progress/elevation/')

####TO DO: write raw_data folder contents####
#raw_data folder contents
#

#folder structure etc
dir.create('./plots/')
dir.create('./output/')
dir.create('./output/mammals/')
dir.create('./output/birds/')
dir.create('./output/world/')
dir.create('./output/mammals/trees/')
dir.create('./output/birds/trees/')
dir.create('./output/mammals/tables/')
dir.create('./output/mammals/tables/DR_posterior/')
dir.create('./output/mammals/tables/DR_posterior/input_tables/')
dir.create('./output/birds/tables/')
dir.create('./output/birds/tables/DR_posterior/')
dir.create('./output/birds/tables/DR_posterior/input_tables/')
dir.create('./output/world/tables/')
dir.create('./output/world/sems/')
dir.create('./output/world/sems/pseudoposterior/')
dir.create('./output/world/sems/pseudoposterior/global/')
dir.create('./output/world/sems/pseudoposterior/gain_elevation/')
dir.create('./output/world/sems/pseudoposterior/loss_elevation/')

####---A) DATA PREP---####
####1) PHYLOGENETIC DATA PREP####
####this is taken from github.com/javierigea/hotspots_mambird_paper/
####prepare MCCtree for mammals####
####Rolland et al Plos Biol: they build a MCC tree with the 100 trees in Fritz file, but they use -keep with the nodes. MCCtree is almost equal to tree #23 in the Fritz 100 file
####system('/Applications2/BEAST/1.8.2/bin/treeannotator ./raw_data/mammals/trees/FritzTree.rs200k.100trees.tre ./output/mammals/trees/mammals_MCC_keep.nexus')
#for mammals
dir.create('./output/mammals/trees/')
system('/Applications2/BEAST/1.8.2/bin/treeannotator -heights ca ./raw_data/mammals/trees/FritzTree.rs200k.100trees.tre ./output/mammals/trees/mammals_MCC.nexus')
#converting MCC from Treeannotator (nexus) to Newick
nexus.tree<-read.nexus('./output/mammals/trees/mammals_MCC.nexus')
write.tree(nexus.tree,'./output/mammals/trees/mammals_MCC.tree')
#calibrate the MCC tree with Meredith dates (as in Rolland Plos Biol)
source('./R/calibrate_mammals_trees.R')
calibrated.tree.MCC<-calibrate_tree_Meredith(treefile='./output/mammals/trees/mammals_MCC.tree')
write.tree(calibrated.tree.MCC,file='./output/mammals/trees/mammals_MCC_calibrated.tree')

####use taxonomy to deal with synonyms for mammals####
#read in trees and use IUCN taxonomy to remove synonyms
mammal.tree<-read.tree('./output/mammals/trees/mammals_MCC_calibrated.tree')
source('./R/taxonomy_IUCN_birdlife.R')
mammal.tree.IUCN<-mammal.tree
mammal.tree.IUCN$tip.label<-iucn.taxonomy.synonyms(species=mammal.tree$tip.label,dictionaryfile = './raw_data/IUCNTaxonomy_Mammalia_30082017_version2017-1.csv')
#remove duplicated tips after correcting synonyms
tip.duplicated.mammal.tree.IUCN<-return.duplicate.tips(mammal.tree.IUCN)
mammal.tree.IUCN<-drop.tip(mammal.tree.IUCN,tip.duplicated.mammal.tree.IUCN)
write.tree(mammal.tree.IUCN,'./output/mammals/trees/mammals_tree_IUCN.tree')

####prepare MCCtree for birds####
dir.create('./output/birds/trees/')
#get 100 random trees from pseudoposterior
trees<-scan('./raw_data/birds/trees/BirdzillaEricsonAllTrees.tre',sep='\n',what='char')
#gets 100 trees evenly distributed from the pseudoposterior
trees.sub<-trees[round(seq(1,length(trees),length.out=100))]
writeLines(trees.sub,'./output/birds/trees/BirdzillaEricsonAllTrees_100.trees')
trees100<-read.tree('./output/birds/trees/BirdzillaEricsonAllTrees_100.trees')
writeNexus(trees100,'./output/birds/trees//BirdzillaEricsonAllTrees_100.nex')
system('/Applications2/BEAST/1.8.2/bin/treeannotator -heights ca ./output/birds/trees/BirdzillaEricsonAllTrees_100.nex ./output/birds/trees/birds_MCC.nexus')
#converting MCC from Treeannotator (nexus) to Newick
nexus.tree<-read.nexus('./output/birds/trees/birds_MCC.nexus')
write.tree(nexus.tree,'./output/birds/trees/birds_MCC.tree')

####use taxonomy to deal with synonyms for birds####
bird.tree<-read.tree('./output/birds/trees/birds_MCC.tree')
source('./R/taxonomy_IUCN_birdlife.R')
bird.tree.IUCN<-bird.tree
bird.tree.IUCN$tip.label<-iucn.taxonomy.synonyms(species=bird.tree$tip.label,dictionaryfile = './raw_data/IUCNTaxonomy_Aves_30082017_version2017-1.csv')
#remove duplicated tips after correcting synonyms
tip.duplicated.bird.tree.IUCN<-return.duplicate.tips(bird.tree.IUCN)
bird.tree.IUCN<-drop.tip(bird.tree.IUCN,tip.duplicated.bird.tree.IUCN)
write.tree(bird.tree.IUCN,'./output/birds/trees/birds_tree_IUCN.tree')

####mammal pseudoposterior time calibration####
#100 trees from the mammal pseudoposterior####
source('./R/calibrate_mammals_trees.R')
calibrate_pseudoposterior_Meredith(treesfile='./raw_data/mammals/trees/FritzTree.rs200k.100trees.tre')
#use IUCN taxonomy on posterior calibrated trees
trees100<-read.nexus('./output/mammals/trees/posterior_calibrated/FritzTree.rs200k.100trees_Meredithdates.trees')
source('./R/taxonomy_IUCN_birdlife.R')
trees100.IUCN<-lapply(trees100,function(x){x$tip.label<-iucn.taxonomy.synonyms(species=x$tip.label,dictionaryfile = './raw_data/IUCNTaxonomy_Mammalia_30082017_version2017-1.csv');tip.duplicated.mammal.tree.IUCN<-return.duplicate.tips(x);x<-drop.tip(x,tip.duplicated.mammal.tree.IUCN);return(x)})
write.nexus(trees100.IUCN,file='./output/mammals/trees/posterior_calibrated/FritzTree.rs200k.100trees_Meredithdates_IUCN.trees')

####bird pseudoposterior####
####use IUCN taxonomy on Jetz pseudoposterior 100 trees
bird100trees<-read.tree('./output/birds/trees/BirdzillaEricsonAllTrees_100.trees')
bird100trees.IUCN<-lapply(bird100trees,function(x){x$tip.label<-iucn.taxonomy.synonyms(species=x$tip.label,dictionaryfile = './raw_data/IUCNTaxonomy_Aves_30082017_version2017-1.csv');tip.duplicated.bird.tree.IUCN<-return.duplicate.tips(x);x<-drop.tip(x,tip.duplicated.bird.tree.IUCN);return(x)})
write.nexus(bird100trees.IUCN,file='./output/birds/trees/BirdzillaEricsonAllTrees_100_IUCN.trees')

####2)SPECIES SPATIAL DATA PREP####
####this is taken from github.com/javierigea/hotspots_mambird_paper/
#####overlap mammal IUCN ranges layers with grid####
#run './R/overlap_realms_mammals_grid_cluster.R' on hydrogen
#/scripts/conscriptoR /home/ji247/hotspots_vertebrates/overlap_realms_mammals_grid_cluster.R -p32 (#for 32 cores,it takes ~20 hours)
#move the *_realms_species_gridoccurrence_table.txt and *_richness_grid_table.txt to ./output/mammals/tables/
#merge realm based data to world level
source('./R/merging_all_realm_occurrences.R')
merge_realms_speciesdata_into_world(path='./output/mammals/tables/')
#copy grid_*100.rds to /output/grids/
merge_realms_grid_into_world(path='./output/grids/')
#then use IUCN taxonomy on them (this is redundant)
source('./R/taxonomy_IUCN_birdlife.R')
mammals_grid.table<-read.table('./output/mammals/tables/100_all_realms_species_gridoccurrence_table.txt',header=T,sep='\t',stringsAsFactors = F)
mammals_grid.table.IUCN<-iucn.taxonomy.synonyms(species=mammals_grid.table$spp,dictionaryfile = './raw_data/IUCNTaxonomy_Mammalia_30082017_version2017-1.csv')
mammals_grid.table$spp<-mammals_grid.table.IUCN
#deal with duplicates: merge the ranges
duplicated.mammals_grid.table.spp<-mammals_grid.table$spp[duplicated(mammals_grid.table$spp)]
if(length(duplicated.birds_grid.table.spp)>0){
  mammals_grid.table<-merge_duplicates_in_speciesgridoccurrence_table(duplicated.species = duplicated.mammals_grid.table.spp,species.gridoccurrence.table = mammals_grid.table)  
}
write.table(mammals_grid.table,'./output/mammals/tables/mammals_100_all_realms_species_gridoccurrence_table.txt',sep='\t',quote=F,row.names=F)

#####overlap bird birdlife range layers with grid####
#BOTW.gbd has to be on './raw_data/birds/BOTW/BOTW.gdb/'; run './R/run_split_birds_gb_into_shp.R
#this will create folders 1,501,1001,1501...18001, each with 500 shapefiles
source('./R/run_split_birds_gdb_into_shp.R')
#zip + copy those folders to cluster
#run './R/overlap_realms_birds_grid_cluster.R' on hydrogen
#/scripts/conscriptoR /home/ji247/hotspots_vertebrates/overlap_realms_birds_grid_cluster.R -p48 (#for 48 cores,it takes 4-5 days)
#and then run commented final parts of overlap_realms_birds_grid_cluster.R
#move the *_realms_species_gridoccurrence_table.txt and *_richness_grid_table.txt to ./output/birds/tables/
#merge realm based data to world level
source('./R/merging_all_realm_occurrences.R')
merge_realms_speciesdata_into_world(path='./output/birds/tables/')
#copy grid_*100.rds to /output/grids/
merge_realms_grid_into_world(path='./output/grids/')
#then use IUCN taxonomy on them (this is redundant)
source('./R/taxonomy_IUCN_birdlife.R')
birds_grid.table<-read.table('./output/birds/tables/100_all_realms_species_gridoccurrence_table.txt',header=T,sep='\t',stringsAsFactors = F)
birds_grid.table.IUCN<-iucn.taxonomy.synonyms(species=birds_grid.table$spp,dictionaryfile = './raw_data/IUCNTaxonomy_Aves_30082017_version2017-1.csv')
birds_grid.table$spp<-birds_grid.table.IUCN
#deal with duplicates: merge the ranges
duplicated.birds_grid.table.spp<-birds_grid.table$spp[duplicated(birds_grid.table$spp)]
if(length(duplicated.birds_grid.table.spp)>0){
  birds_grid.table<-merge_duplicates_in_speciesgridoccurrence_table(duplicated.species = duplicated.birds_grid.table.spp,species.gridoccurrence.table = birds_grid.table)  
}
write.table(birds_grid.table,'./output/birds/tables/birds_100_all_realms_species_gridoccurrence_table.txt',sep='\t',quote=F,row.names=F)

#correct a mistake with Strix_butleri 131632469 should be 13163 2469
birds.ranges<-read.table('./output/birds/tables/birds_100_all_realms_species_gridoccurrence_table.txt',header=T,sep='\t',stringsAsFactors = F)
birds.ranges[birds.ranges$spp=='Strix_butleri','cells']<-gsub(birds.ranges[birds.ranges$spp=='Strix_butleri','cells'],pattern='131632469',replacement='13163 2469')
write.table(birds.ranges, './output/birds/tables/birds_100_all_realms_species_gridoccurrence_table.txt',sep='\t',quote=F,row.names=F)

####3) DIVERSIFICATION RATE METRICS (DR,BAMM)#####
####this is taken from github.com/javierigea/hotspots_mambird_paper/
####measure DR in mammal trees####
source('./R/measure_DR.R')
#for mammals
DR.mammals<-measure_DR_tree_table(treefile = "./output/mammals/mammals_tree_IUCN.tree")
#select terrestrial
DR.mammals.marine<-get.marine.species(species = DR.mammals$Species,dictionaryfile = './raw_data/IUCNTaxonomy_Mammalia_30082017_version2017-1.csv')
DR.mammals.marine<-c(DR.mammals.marine,'Platanista_minor')
DR.mammals.terrestrial<-DR.mammals[!(DR.mammals$Species%in%DR.mammals.marine),]
write.table(DR.mammals,'./output/mammals/DR_mammals_tree_IUCN.txt',sep='\t',quote=F,row.names=F)
write.table(DR.mammals.terrestrial,'./output/mammals/DR_mammals_terrestrial_tree_IUCN.txt',sep='\t',quote=F,row.names=F)
#for mammals pseudoposterior
mammals100trees<-read.nexus('./output/mammals/trees/posterior_calibrated/FritzTree.rs200k.100trees_Meredithdates_IUCN.trees')
DR.mammals.100<-lapply(mammals100trees,function(x)measure_DR_tree_table(x))
#select terrestrial
DR.mammals.marine<-get.marine.species(species = DR.mammals.100[[1]]$Species,dictionaryfile = './raw_data/IUCNTaxonomy_Mammalia_30082017_version2017-1.csv')
DR.mammals.marine<-c(DR.mammals.marine,'Platanista_minor')
DR.mammals.terrestrial.100<-lapply(DR.mammals.100,function(x)x[!(x$Species%in%DR.mammals.marine),])
DR.mammals.terrestrial<-DR.mammals[!(DR.mammals$Species%in%DR.mammals.marine),]
#save object
saveRDS(DR.mammals.terrestrial.100,file='./output/mammals/trees/DR.mammals.terrestrial.100.rds')
#write table with 100 tables
lapply(c(1:100), function(x) write.table(DR.mammals.terrestrial.100[[x]], file = paste0('./output/mammals/tables/DR_posterior/input_tables/DR_mammals_terrestrial_tree_IUCN_', x, '.txt'), sep = '\t', quote = F, row.names = F))

#get median DR across the pseudoposterior for each species
get_pseudoposterior_median_DRtable(DR.pseudoposterior.file='./output/mammals/trees/DR.mammals.terrestrial.100.rds',path='./output/mammals/tables/',name='mammals')
#compare DR from the MCC tree with median of the pseudoposterior
#this is part of Fig. S2
source('./R/measure_DR.R')
pdf('./output/mammals/plots/mammals_DRMCC_vs_pseudoposterior.pdf')
compare_DR_MCC_vs_medianpseudoposterior(DR.tablefile='./output/mammals/DR_mammals_terrestrial_tree_IUCN.txt',DR.pseudoposterior.tablefile='./all_realms_new_pseudoposteriorDRmedian.txt')
compare_DR_MCC_vs_allpseudoposterior(DR.tablefile = './output/mammals/DR_mammals_terrestrial_tree_IUCN.txt',DR.pseudoposterior.file ='./output/mammals/trees/DR.mammals.terrestrial.100.rds' )
dev.off()

####measure DR in bird trees####
source('./R/measure_DR.R')
#for birds
DR.birds<-measure_DR_tree_table(treefile = "./output/birds/birds_tree_IUCN.tree")
write.table(DR.birds,'./output/birds/DR_birds_tree_IUCN.txt',sep='\t',quote=F,row.names=F)
#for birds pseudoposterior
birds100trees<-read.nexus('./output/birds/trees/BirdzillaEricsonAllTrees_100_IUCN.trees')
counter<-0
DR.birds.100<-lapply(birds100trees,function(x){counter<<-counter+1;cat(counter,'\n');measure_DR_tree_table(x)})
#save object
saveRDS(DR.birds.100,file='./output/birds/trees/DR.birds.100.rds')
#write table with 100 tables
lapply(c(1:100), function(x) write.table(DR.birds.100[[x]], file = paste0('./output/birds/tables/DR_posterior/input_tables/DR_birds_tree_IUCN_', x, '.txt'), sep = '\t', quote = F, row.names = F))

#get median DR across the pseudoposterior for each species
get_pseudoposterior_median_DRtable(DR.pseudoposterior.file='./output/birds/trees/DR.birds.100.rds',path='./output/birds/tables/',name='all_realms_birds')
#compare DR from the MCC tree with median of the pseudoposterior
#this is part of Fig. S2
source('./R/measure_DR.R')
pdf('./output/birds/plots/birds_DRMCC_vs_pseudoposterior.pdf')
compare_DR_MCC_vs_pseudoposterior(DR.tablefile='./output/birds/DR_birds_tree_IUCN.txt',DR.pseudoposterior.tablefile='./output/birds/tables/all_realms_birds_pseudoposteriorDRmedian.txt')
compare_DR_MCC_vs_allpseudoposterior(DR.tablefile = './output/birds/DR_birds_tree_IUCN.txt',DR.pseudoposterior.file ='./output/birds/trees/DR.birds.100.rds')
dev.off()

####BAMM analyses####
#prepare BAMM input
source('./R/BAMM_functions.R')
prepare_BAMM_input(treefile='./output/mammals/trees/mammals_tree_IUCN.tree',dictionaryfile = './raw_data/IUCNTaxonomy_Mammalia_30082017_version2017-1.csv',path='./output/mammals/trees/',name='mammals_IUCN_BAMM')
prepare_BAMM_input(treefile='./output/birds/birds_tree_IUCN.tree',dictionaryfile = './raw_data/IUCNTaxonomy_Aves_30082017_version2017-1.csv',path='./output/birds/trees/',name='birds_IUCN_BAMM')  
#for mammals it takes ~40 hours (50 million generations), more than enough for convergence (ESS:1200, could work with much less generations) (it took ~80 hours on 4 cores on node 8 - hydrogen)
#for mammals it takes ~21 hours for 30 million generations; for birds it takes 64 hours for 30 million generations
#check convergence
analyse_BAMM_convergence(mcmcout = './output/mammals/trees/mcmc_out_mammals_30m.txt',burnin=0.25)
analyse_BAMM_convergence(mcmcout = './output/birds/trees/mcmc_out_birds_30m.txt',burnin=0.25)
#get TipRates
get_tipRates_BAMM(treefile='./output/mammals/trees/mammals_tree_IUCN.tree',eventfile='./output/mammals/trees/event_data_mammals_30m.txt',burnin=0.25,path='./output/mammals/tables/',name='mammals_all_realms')
get_tipRates_BAMM(treefile='./output/birds/trees/birds_tree_IUCN.tree',eventfile='./output/birds/trees/event_data_birds_30m.txt',burnin=0.25,path='./output/birds/tables/',name='birds_all_realms')

####4) DIVERSIFICATION RATES IN SPACE####
####mammal DR metrics in space####
source('./R/measure_DR.R')
#read mammals ranges
mammals.ranges<-read.table('./output/mammals/tables/mammals_100_all_realms_species_gridoccurrence_table.txt',header=T,sep='\t',stringsAsFactors = F)
list.mammals.ranges<-lapply(mammals.ranges$cells,function(x){char<-unlist(strsplit(as.character(x),' '));char<-char[char!=''];as.numeric(char)})
names(list.mammals.ranges)<-as.character(mammals.ranges$spp)
#read mammals DR table
mammalsDR<-read.table('./output/mammals/tables/DR_mammals_terrestrial_tree_IUCN.txt',header=T,sep='\t',stringsAsFactors = F)
#this takes ca 2 hours
mammals.DR.grid.table<-DR_stats_grid(list.species.ranges=list.mammals.ranges,speciesDR=mammalsDR)
write.table(mammals.DR.grid.table,file='./output/mammals/tables/mammals_DR_cells_table.txt',sep='\t',quote=F,row.names=F)

####mammal pseudoposterior DR metrics in space####
#this is to be run on hydrogen
#source('./R/measure_DR.R')
#run './R/measure_DRmammals_pseudo.R 1' on hydrogen
#copy tables to './output/mammals/tables/DR_posterior/'

#checking the correlation of wDR with MCC tree and wDR across the pseudoposterior
DRmammals <- read.table('./output/mammals/tables/mammals_DR_cells_table.txt',
                 header = T, 
                 sep = '\t', 
                 stringsAsFactors = F)

DRpseudomammals <- lapply (list.files('./output/mammals/tables/DR_posterior/', '_cells_table.txt'), function(x) read.table(paste0('../../elevation/output/mammals/tables/DR_posterior/',x), sep = '\t', header = T, stringsAsFactors = F))
hist(unlist(lapply(DRpseudomammals, function(x) cor.test(x$mean.wDR, DRmammals$mean.wDR, method = 's')$estimate)))
#very strong positive correlation


####bird DR metrics in space####
source('./R/measure_DR.R')
#read birds ranges
birds.ranges<-read.table('./output/birds/tables/birds_100_all_realms_species_gridoccurrence_table.txt',header=T,sep='\t',stringsAsFactors = F)
list.birds.ranges<-lapply(birds.ranges$cells,function(x){char<-unlist(strsplit(as.character(x),' '));char<-char[char!=''];as.numeric(char)})
names(list.birds.ranges)<-as.character(birds.ranges$spp)
#read birds DR table
birdsDR<-read.table('./output/birds/tables/DR_birds_terrestrial_tree_IUCN.txt',header=T,sep='\t',stringsAsFactors = F)
#this takes ca 2 hours
birds.DR.grid.table<-DR_stats_grid(list.species.ranges=list.birds.ranges,speciesDR=birdsDR)
write.table(birds.DR.grid.table,file='./output/birds/tables/birds_DR_cells_table.txt',sep='\t',quote=F,row.names=F)

####bird pseudoposterior DR metrics in space####
#this is to be run on hydrogen
#source('./R/measure_DR.R')
#run './R/measure_DRbirds_pseudo.R 1' on hydrogen
#copy tables to './output/birds/tables/DR_posterior/'

#checking the correlation of wDR with MCC tree and wDR across the pseudoposterior
DRbirds <- read.table('./output/birds/tables/birds_DR_cells_table.txt',
                      header = T, 
                      sep = '\t', 
                      stringsAsFactors = F)

DRpseudobirds <- lapply (list.files('./output/birds/tables/DR_posterior/', '_cells_table.txt'), function(x) read.table(paste0('../../elevation/output/birds/tables/DR_posterior/',x), sep = '\t', header = T, stringsAsFactors = F))
hist(unlist(lapply(DRpseudobirds, function(x) cor.test(x$mean.wDR, DRbirds$mean.wDR, method = 's')$estimate)))
#very strong positive correlation

####mammal BAMM metrics in space####
source('./R/BAMM_functions.R')
#read mammals ranges
mammals.ranges<-read.table('./output/mammals/tables/mammals_100_all_realms_species_gridoccurrence_table.txt',header=T,sep='\t',stringsAsFactors = F)
list.mammals.ranges<-lapply(mammals.ranges$cells,function(x){char<-unlist(strsplit(as.character(x),' '));char<-char[char!=''];as.numeric(char)})
names(list.mammals.ranges)<-as.character(mammals.ranges$spp)
#read mammals BAMM table
mammalsBAMM<-read.table('./output/mammals/tables/BAMM_mammals_terrestrial_tree_IUCN.txt',header=T,sep='\t',stringsAsFactors = F)
#this takes ca 2 hours
mammals.BAMM.grid.table<-BAMM_stats_grid(list.species.ranges=list.mammals.ranges,speciesBAMM=mammalsBAMM)
write.table(mammals.BAMM.grid.table,file='./output/mammals/tables/mammals_BAMM_cells_table.txt',sep='\t',quote=F,row.names=F)

####bird BAMM metrics in space####
source('./R/BAMM_functions.R')
#read birds ranges
birds.ranges<-read.table('./output/birds/tables/birds_100_all_realms_species_gridoccurrence_table.txt',header=T,sep='\t',stringsAsFactors = F)
#correct a mistake with Strix_butleri 131632469 should be 13163 2469
birds.ranges[birds.ranges$spp=='Strix_butleri','cells']<-gsub(birds.ranges[birds.ranges$spp=='Strix_butleri','cells'],pattern='131632469',replacement='13163 2469')

list.birds.ranges<-lapply(birds.ranges$cells,function(x){char<-unlist(strsplit(as.character(x),' '));char<-char[char!=''];as.numeric(char)})
names(list.birds.ranges)<-as.character(birds.ranges$spp)
#read birds BAMM table
birdsBAMM<-read.table('./output/birds/tables/BAMM_birds_terrestrial_tree_IUCN.txt',header=T,sep='\t',stringsAsFactors = F)
#this takes ca 2 hours
birds.BAMM.grid.table<-BAMM_stats_grid(list.species.ranges=list.birds.ranges,speciesBAMM=birdsBAMM)
write.table(birds.BAMM.grid.table,file='./output/birds/tables/birds_BAMM_cells_table.txt',sep='\t',quote=F,row.names=F)

####5) ELEVATION DATA####
source('./R/elevation_get_data.R')
####present day elevation with ETOPO####
##downloaded from https://www.ngdc.noaa.gov/mgg/global/ (ETOPO1 Ice Surface)
gridWorld<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
present.elevation.table<-get_present_elevation_grid(etopo1.path='./raw_data/rasters/ETOPO1_Ice_g_geotiff.tif',grid=gridWorld)
write.table(present.elevation.table,file='./output/world/tables/grid_present_elevation_table.txt',sep='\t',quote=F,row.names=F)

####past elevation data with PRISM4####
gridWorld<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
PRISM4.elevation.table<-get_past_elevation_grid(PRISM4.path='./raw_data/rasters/Plio_enh_topo_v1.0_PRISM4.nc',grid=gridWorld)
write.table(PRISM4.elevation.table,file='./output/world/tables/grid_PRISM4_agg_elevation_table.txt',sep='\t',quote=F,row.names=F)

####6) TEMPERATURE DATA####
source('./R/temperature_get_data.R')

####present temperature data####
#get present temperature values for grid
gridWorld<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
present.temperature.table<-get_present_temperature_grid(bio1raster.path='./raw_data/rasters/current_2_5m/bio1.bil',grid=gridWorld)
write.table(present.temperature.table,file='./output/world/tables/grid_present_temperature_table.txt',sep='\t',quote=F,row.names=F)

####past temperature data####
#get past temperature values for grid
gridWorld<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
past.temperature.table<-get_past_temperature_grid(bio1raster.past.path='./raw_data/rasters/paleoclim_2_5min/bio_1.tif',grid=gridWorld)
write.table(past.temperature.table,file='./output/world/tables/grid_past_temperature_table.txt',sep='\t',quote=F,row.names=F)

#####TO DO: delete this below####
#substract bio1 present minus bio1 past
#bio1_change<-substract_current_minus_past_bio1(path.present='./raw_data/rasters/current_2_5m/bio1.bil',path.past='./raw_data/rasters/paleoclim_2_5min/bio_1.tif')
#writeRaster(bio1_change,'./output/world/bio1_present_minus_past_degreescent')
#get change in temperature values for grid
#gridWorld<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
#present.minus.past.temperature.table<-get_change_temperature_grid(raster.path='./output/world/bio1_present_minus_past_degreescent',grid=gridWorld)
#write.table(past.temperature.table,file='./output/world/tables/grid_past_temperature_table.txt',sep='\t',quote=F,row.names=F)



####7)AGGREGATE ALL DIVERSIFICATION AND ELEVATION DATA####
##open DR & BAMM tables for mammals and birds
DRmammals<-read.table('./output/mammals/tables/mammals_DR_cells_table.txt',header=T,sep='\t',stringsAsFactors = F)
colnames(DRmammals)[c(2,3,4)]<-paste0('mammals.',colnames(DRmammals)[c(2,3,4)])
DRbirds<-read.table('./output/birds/tables/birds_DR_cells_table.txt',header=T,sep='\t',stringsAsFactors = F)
colnames(DRbirds)[c(2,3,4)]<-paste0('birds.',colnames(DRbirds)[c(2,3,4)])
BAMMmammals<-read.table('./output/mammals/tables/mammals_BAMM_cells_table.txt',header=T,sep='\t',stringsAsFactors = F)
colnames(BAMMmammals)[c(2,3,4)]<-paste0('mammals.',colnames(BAMMmammals)[c(2,3,4)])
BAMMbirds<-read.table('./output/birds/tables/birds_BAMM_cells_table.txt',header=T,sep='\t',stringsAsFactors = F)
colnames(BAMMbirds)[c(2,3,4)]<-paste0('birds.',colnames(BAMMbirds)[c(2,3,4)])
#open elevation and temperature tables
presentelevation<-read.table('./output/world/tables/grid_present_elevation_table.txt',header=T,sep='\t',stringsAsFactors = F)
pastelevation<-read.table('./output/world/tables/grid_PRISM4_elevation_table.txt',header=T,sep='\t',stringsAsFactors = F)
presenttemp<-read.table('./output/world/tables/grid_present_temperature_table.txt',header=T,sep='\t',stringsAsFactors = F)
pasttemp<-read.table('./output/world/tables/grid_past_temperature_table.txt',header=T,sep='\t',stringsAsFactors = F)

#combine all tables
Div.all<-Reduce(function(x,y) merge(x=x,y=y,by = "cells"),list(DRmammals, DRbirds, BAMMmammals,BAMMbirds,presentelevation,pastelevation,presenttemp,pasttemp))
#divide temperature by 10 (200 = 20degrees)
Div.all$mean.present.T<-Div.all$mean.present.T/10
Div.all$range.present.T<-Div.all$range.present.T/10
Div.all$mean.past.T<-Div.all$mean.past.T/10
Div.all$range.past.T<-Div.all$range.past.T/10

#calculate changes in elevation and in temperature and gains and losses
Div.all$present.minus.past.elevation<-Div.all$mean.elevation.ETOPO.land-Div.all$mean.elevation.PRISM4.land
Div.all$present.minus.past.temperature<-Div.all$mean.present.T-Div.all$mean.past.T

#calculate gains and losses in elevation
Div.all$elevation.gain<-Div.all$present.minus.past.elevation
Div.all[!is.na(Div.all$elevation.gain)&Div.all$elevation.gain<0,'elevation.gain']<-0
Div.all$elevation.loss<-Div.all$present.minus.past.elevation
Div.all[!is.na(Div.all$elevation.loss)&Div.all$elevation.loss>0,'elevation.loss']<-0
Div.all$elevation.loss<-abs(Div.all$elevation.loss)
write.table(Div.all,file='./output/all_variables_grid_table.txt',sep='\t',quote=F,row.names = F)
####---B) ANALYSES ----####

####B1) linear/spatial autocorrelation models####
####TO DO: put this in a function####
library(corrplot)
#prepare data
#linear models of speciation with present/past/change elevation/temperature
#with wDR
cells.table<-read.table('./output/all_variables_grid_table.txt',sep='\t',header=T,stringsAsFactors = F)
#drop cells where speciation = 0 (there are no species)
cells.table<-cells.table[!is.na(cells.table$mammals.mean.DR)&!is.na(cells.table$birds.mean.DR)&cells.table$mammals.mean.DR>0&cells.table$birds.mean.DR>0,]
#drop cells where past or present elevation is NA or negative
cells.table<-cells.table[!is.na(cells.table$mean.elevation.ETOPO.land)&!is.na(cells.table$mean.elevation.PRISM4.land)&cells.table$mean.elevation.ETOPO.land>0&cells.table$mean.elevation.PRISM4.land>0,]
#drop unnecessary columns
cells.table<-cells.table[,c('cells','mammals.mean.wDR','birds.mean.wDR','mean.elevation.ETOPO.land','mean.elevation.PRISM4.land','mean.present.T','mean.past.T','present.minus.past.elevation','present.minus.past.temperature')]
cells.table<-cells.table[complete.cases(cells.table),]
#check correlations among predictors
corrplot(cor(cells.table[,c(2:9)]),method = 'number')
#log transform response variables
cells.table$mammals.mean.wDR<-log(cells.table$mammals.mean.wDR)
cells.table$birds.mean.wDR<-log(cells.table$birds.mean.wDR)
#log transform other variables with positive values
cells.table$mean.elevation.ETOPO.land<-log(cells.table$mean.elevation.ETOPO.land)
cells.table$mean.elevation.PRISM4.land<-log(cells.table$mean.elevation.PRISM4.land)
#scale predictors
cells.table[,c(4:9)]<-apply(cells.table[,c(4:9)],2,function(x)scale(x))


####multiple linear regressions for wDR####
#predicting wDR in mammals and birds using present elevation, change in elevation, present temperature and change in temperature
lm.mammals.wDR.elevation.temp<-lm(mammals.mean.wDR~mean.elevation.ETOPO.land+present.minus.past.elevation+mean.present.T+present.minus.past.temperature,data=cells.table)
lm.birds.wDR.elevation.temp<-lm(birds.mean.wDR~mean.elevation.ETOPO.land+present.minus.past.elevation+mean.present.T+present.minus.past.temperature,data=cells.table)

#linear models of speciation with present/past/change elevation/temperature
#with wlambda.avg
cells.table<-read.table('./output/all_variables_grid_table.txt',sep='\t',header=T,stringsAsFactors = F)
#drop cells where speciation = 0 (there are no species)
cells.table<-cells.table[!is.na(cells.table$mammals.mean.DR)&!is.na(cells.table$birds.mean.DR)&cells.table$mammals.mean.DR>0&cells.table$birds.mean.DR>0,]
#drop cells where past or present elevation is NA or negative
cells.table<-cells.table[!is.na(cells.table$mean.elevation.ETOPO.land)&!is.na(cells.table$mean.elevation.PRISM4.land)&cells.table$mean.elevation.ETOPO.land>0&cells.table$mean.elevation.PRISM4.land>0,]
#drop unnecessary columns
cells.table<-cells.table[,c('cells','mammals.mean.wlambda.avg','birds.mean.wlambda.avg','mean.elevation.ETOPO.land','mean.elevation.PRISM4.land','mean.present.T','mean.past.T','present.minus.past.elevation','present.minus.past.temperature')]
cells.table<-cells.table[complete.cases(cells.table),]
#check correlations among predictors
corrplot(cor(cells.table[,c(2:9)]),method = 'number')
#log transform response variables
cells.table$mammals.mean.wlambda.avg<-log(cells.table$mammals.mean.wlambda.avg)
cells.table$birds.mean.wlambda.avg<-log(cells.table$birds.mean.wlambda.avg)
#log transform other variables with positive values
cells.table$mean.elevation.ETOPO.land<-log(cells.table$mean.elevation.ETOPO.land)
cells.table$mean.elevation.PRISM4.land<-log(cells.table$mean.elevation.PRISM4.land)
#scale predictors
cells.table[,c(4:9)]<-apply(cells.table[,c(4:9)],2,function(x)scale(x))


####multiple linear regressions for wlambda.avg####
#predicting wlambda.avg in mammals and birds using present elevation, change in elevation, present temperature and change in temperature
lm.mammals.wlambda.avg.elevation.temp<-lm(mammals.mean.wlambda.avg~mean.elevation.ETOPO.land+present.minus.past.elevation+mean.present.T+present.minus.past.temperature,data=cells.table)
lm.birds.wlambda.avg.elevation.temp<-lm(birds.mean.wlambda.avg~mean.elevation.ETOPO.land+present.minus.past.elevation+mean.present.T+present.minus.past.temperature,data=cells.table)

####assessing spatial autocorrelation in residuals####
#for wDR
#load grid
gridWorld<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
#get coordinates
grid.world.longlat<-lapply(gridWorld,function(x) spTransform(x,CRS("+proj=longlat")))
grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) coordinates(x))
grid.coordinates<-do.call("rbind", grid.coordinates)
#get neighbours in 1000km distqnce
neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)
#significant spatial autocorrelation in the residuals (Moran's I)
moran.lm.mammals.wDR.elevation.temp<-moran.test(lm.mammals.wDR.elevation.temp$residuals,listw = neighbours.1000.w,zero.policy = TRUE)
moran.lm.birds.wDR.elevation.temp<-moran.test(lm.birds.wDR.elevation.temp$residuals,listw = neighbours.1000.w,zero.policy = TRUE)
moran.mc.lm.mammals.wDR.elevation.temp<-moran.mc(lm.mammals.wDR.elevation.temp$residuals,listw = neighbours.1000.w,nsim = 1000,zero.policy = TRUE)
moran.mc.lm.birds.wDR.elevation.temp<-moran.mc(lm.birds.wDR.elevation.temp$residuals,listw = neighbours.1000.w,nsim = 1000,zero.policy = TRUE)

#for wlambda.avg
#load grid
gridWorld<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
#get coordinates
grid.world.longlat<-lapply(gridWorld,function(x) spTransform(x,CRS("+proj=longlat")))
grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) coordinates(x))
grid.coordinates<-do.call("rbind", grid.coordinates)
#get neighbours in 1000km distqnce
neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)
#significant spatial autocorrelation in the residuals (Moran's I)
moran.lm.mammals.wlambda.avg.elevation.temp<-moran.test(lm.mammals.wlambda.avg.elevation.temp$residuals,listw = neighbours.1000.w,zero.policy = TRUE)
moran.lm.birds.wlambda.avg.elevation.temp<-moran.test(lm.birds.wlambda.avg.elevation.temp$residuals,listw = neighbours.1000.w,zero.policy = TRUE)
moran.mc.lm.mammals.wlambda.avg.elevation.temp<-moran.mc(lm.mammals.wlambda.avg.elevation.temp$residuals,listw = neighbours.1000.w,nsim = 1000,zero.policy = TRUE)
moran.mc.lm.birds.wlambda.avg.elevation.temp<-moran.mc(lm.birds.wlambda.avg.elevation.temp$residuals,listw = neighbours.1000.w,nsim = 1000,zero.policy = TRUE)


####spatial autocorrelation models####
#load grid
grid.world<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
#get coordinates
grid.world.longlat<-lapply(grid.world,function(x) spTransform(x,CRS("+proj=longlat")))
grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) sp::coordinates(x))
grid.coordinates<-do.call("rbind", grid.coordinates)
#get neighbours in 1000km distqnce
neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)

#for wDR
sarlm.mammals.wDR.elevation.temp<-errorsarlm(mammals.mean.wDR~mean.elevation.ETOPO.land+present.minus.past.elevation+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam')
sarlm.birds.wDR.elevation.temp<-errorsarlm(birds.mean.wDR~mean.elevation.ETOPO.land+present.minus.past.elevation+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam')

#for wlambda.avg
sarlm.mammals.wlambda.avg.elevation.temp<-errorsarlm(mammals.mean.wlambda.avg~mean.elevation.ETOPO.land+present.minus.past.elevation+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam')
sarlm.birds.wlambda.avg.elevation.temp<-errorsarlm(birds.mean.wlambda.avg~mean.elevation.ETOPO.land+present.minus.past.elevation+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam')

####B2) sems####
####assessing spatial autocorrelation in residuals####
#for wDR
#load grid
gridWorld<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
#get coordinates
grid.world.longlat<-lapply(gridWorld,function(x) spTransform(x,CRS("+proj=longlat")))
grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) coordinates(x))
grid.coordinates<-do.call("rbind", grid.coordinates)
#get neighbours in 1000km distqnce
neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)
#built sem models without sarlm
lm.sem.mammals.wDR.elevation.temp<-psem(lm(mammals.mean.wDR~mean.elevation.ETOPO.land+present.minus.past.elevation+mean.present.T+present.minus.past.temperature,data=cells.table),lm(mean.elevation.ETOPO.land~present.minus.past.elevation,data=cells.table),lm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table),lm(present.minus.past.temperature~present.minus.past.elevation,data=cells.table))
lm.sem.birds.wDR.elevation.temp<-psem(lm(birds.mean.wDR~mean.elevation.ETOPO.land+present.minus.past.elevation+mean.present.T+present.minus.past.temperature,data=cells.table),lm(mean.elevation.ETOPO.land~present.minus.past.elevation,data=cells.table),lm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table),lm(present.minus.past.temperature~present.minus.past.elevation,data=cells.table))
#moran's I (for residuals on wDR)
moran.lm.sem.mammals.wDR.elevation.temp<-moran.test(residuals(lm.sem.mammals.wDR.elevation.temp)[,1],listw = neighbours.1000.w,zero.policy = TRUE)
moran.lm.birds.wDR.elevation.temp<-moran.test(residuals(lm.sem.birds.wDR.elevation.temp)[,1],listw = neighbours.1000.w,zero.policy = TRUE)

####B2_1) sems with all cells####

####with wDR
#prepare data
cells.table<-read.table('./output/all_variables_grid_table.txt',sep='\t',header=T,stringsAsFactors = F)
#drop cells where speciation = 0 (there are no species)
cells.table<-cells.table[!is.na(cells.table$mammals.mean.DR)&!is.na(cells.table$birds.mean.DR)&cells.table$mammals.mean.DR>0&cells.table$birds.mean.DR>0,]
#drop cells where past or present elevation is NA or negative
cells.table<-cells.table[!is.na(cells.table$mean.elevation.ETOPO.land)&!is.na(cells.table$mean.elevation.PRISM4.land)&cells.table$mean.elevation.ETOPO.land>0&cells.table$mean.elevation.PRISM4.land>0,]
#drop unnecessary columns
cells.table<-cells.table[,c('cells','mammals.mean.wDR','birds.mean.wDR','mean.elevation.ETOPO.land','mean.elevation.PRISM4.land','mean.present.T','mean.past.T','present.minus.past.elevation','present.minus.past.temperature')]
cells.table<-cells.table[complete.cases(cells.table),]
#check correlations among predictors
corrplot(cor(cells.table[,c(2:9)]),method = 'number')
#log transform response variables
cells.table$mammals.mean.wDR<-log(cells.table$mammals.mean.wDR)
cells.table$birds.mean.wDR<-log(cells.table$birds.mean.wDR)
#log transform other variables with positive values
cells.table$mean.elevation.ETOPO.land<-log(cells.table$mean.elevation.ETOPO.land)
cells.table$mean.elevation.PRISM4.land<-log(cells.table$mean.elevation.PRISM4.land)
#scale predictors
cells.table[,c(4:9)]<-apply(cells.table[,c(4:9)],2,function(x)scale(x))

#load grid
grid.world<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
#get coordinates
grid.world.longlat<-lapply(grid.world,function(x) spTransform(x,CRS("+proj=longlat")))
grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) sp::coordinates(x))
grid.coordinates<-do.call("rbind", grid.coordinates)
#get neighbours in 1000km distqnce
neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)

#sem with sarlm for mammals wDR, elevation, temperature and change elevation and change temperature
sarlm.sem.mammals.wDR.elevation.temp<-psem(errorsarlm(mammals.mean.wDR~mean.elevation.ETOPO.land+present.minus.past.elevation+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.mammals.wDR.elevation.temp,standardize = 'none')
sarlm.sem.birds.wDR.elevation.temp<-psem(errorsarlm(birds.mean.wDR~mean.elevation.ETOPO.land+present.minus.past.elevation+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.birds.wDR.elevation.temp,standardize = 'none')

#for pseudoposteriors: #sem with sarlm for mammals wDR, elevation, temperature and change elevation and change temperature
source('./R/generate_sems.R')
lapply(c(1:100), function(x) sem_spatial_pseudoposterior(variables_table = './output/all_variables_grid_table.txt', 
                            replicate = x))
####with wlambda.avg
#prepare data
cells.table<-read.table('./output/all_variables_grid_table.txt',sep='\t',header=T,stringsAsFactors = F)
#drop cells where speciation = 0 (there are no species)
cells.table<-cells.table[!is.na(cells.table$mammals.mean.DR)&!is.na(cells.table$birds.mean.DR)&cells.table$mammals.mean.DR>0&cells.table$birds.mean.DR>0,]
#drop cells where past or present elevation is NA or negative
cells.table<-cells.table[!is.na(cells.table$mean.elevation.ETOPO.land)&!is.na(cells.table$mean.elevation.PRISM4.land)&cells.table$mean.elevation.ETOPO.land>0&cells.table$mean.elevation.PRISM4.land>0,]
#drop unnecessary columns
cells.table<-cells.table[,c('cells','mammals.mean.wlambda.avg','birds.mean.wlambda.avg','mean.elevation.ETOPO.land','mean.elevation.PRISM4.land','mean.present.T','mean.past.T','present.minus.past.elevation','present.minus.past.temperature')]
cells.table<-cells.table[complete.cases(cells.table),]
#check correlations among predictors
corrplot(cor(cells.table[,c(2:9)]),method = 'number')
#log transform response variables
cells.table$mammals.mean.wlambda.avg<-log(cells.table$mammals.mean.wlambda.avg)
cells.table$birds.mean.wlambda.avg<-log(cells.table$birds.mean.wlambda.avg)
#log transform other variables with positive values
cells.table$mean.elevation.ETOPO.land<-log(cells.table$mean.elevation.ETOPO.land)
cells.table$mean.elevation.PRISM4.land<-log(cells.table$mean.elevation.PRISM4.land)
#scale predictors
cells.table[,c(4:9)]<-apply(cells.table[,c(4:9)],2,function(x)scale(x))

#load grid
grid.world<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
#get coordinates
grid.world.longlat<-lapply(grid.world,function(x) spTransform(x,CRS("+proj=longlat")))
grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) sp::coordinates(x))
grid.coordinates<-do.call("rbind", grid.coordinates)
#get neighbours in 1000km distqnce
neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)

#sem with sarlm for mammals wlambda.avg, elevation, temperature and change elevation and change temperature
sarlm.sem.mammals.wlambda.avg.elevation.temp<-psem(errorsarlm(mammals.mean.wlambda.avg~mean.elevation.ETOPO.land+present.minus.past.elevation+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.mammals.wlambda.avg.elevation.temp,standardize = 'none')
sarlm.sem.birds.wlambda.avg.elevation.temp<-psem(errorsarlm(birds.mean.wlambda.avg~mean.elevation.ETOPO.land+present.minus.past.elevation+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.birds.wlambda.avg.elevation.temp,standardize = 'none')




####B2_2) sems looking at elevation gain only####

####with wDR
#prepare data
cells.table<-read.table('./output/all_variables_grid_table.txt',sep='\t',header=T,stringsAsFactors = F)
#drop cells where speciation = 0 (there are no species)
cells.table<-cells.table[!is.na(cells.table$mammals.mean.DR)&!is.na(cells.table$birds.mean.DR)&cells.table$mammals.mean.DR>0&cells.table$birds.mean.DR>0,]
#drop cells where past or present elevation is NA or negative
cells.table<-cells.table[!is.na(cells.table$mean.elevation.ETOPO.land)&!is.na(cells.table$mean.elevation.PRISM4.land)&cells.table$mean.elevation.ETOPO.land>0&cells.table$mean.elevation.PRISM4.land>0,]
#drop unnecessary columns
cells.table<-cells.table[,c('cells','mammals.mean.wDR','birds.mean.wDR','mean.elevation.ETOPO.land','mean.present.T','elevation.gain','present.minus.past.temperature')]
#get columns where elevation.gain >0
cells.table<-cells.table[cells.table$elevation.gain>0,]
cells.table<-cells.table[complete.cases(cells.table),]
#check correlations among predictors
corrplot(cor(cells.table[,c(2:7)]),method = 'number')
#log transform response variables
cells.table$mammals.mean.wDR<-log(cells.table$mammals.mean.wDR)
cells.table$birds.mean.wDR<-log(cells.table$birds.mean.wDR)
#log transform other variables with positive values
cells.table$mean.elevation.ETOPO.land<-log(cells.table$mean.elevation.ETOPO.land)
cells.table$elevation.gain<-log(cells.table$elevation.gain)
#scale predictors
cells.table[,c(4:7)]<-apply(cells.table[,c(4:7)],2,function(x)scale(x))

#load grid
grid.world<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
#get coordinates
grid.world.longlat<-lapply(grid.world,function(x) spTransform(x,CRS("+proj=longlat")))
grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) sp::coordinates(x))
grid.coordinates<-do.call("rbind", grid.coordinates)
#get neighbours in 1000km distqnce
neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)

#sem with sarlm for mammals wDR, elevation, temperature and change elevation and change temperature
sarlm.sem.mammals.wDR.elevation.gain.temp<-psem(errorsarlm(mammals.mean.wDR~mean.elevation.ETOPO.land+elevation.gain+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.mammals.wDR.elevation.gain.temp,standardize = 'none')
sarlm.sem.birds.wDR.elevation.gain.temp<-psem(errorsarlm(birds.mean.wDR~mean.elevation.ETOPO.land+elevation.gain+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.birds.wDR.elevation.gain.temp,standardize = 'none')

#for pseudoposteriors: #sem with sarlm for mammals wDR, elevation, temperature and change elevation and change temperature
source('./R/generate_sems.R')
lapply(c(1:100), function(x) sem_gain_spatial_pseudoposterior(variables_table = './output/all_variables_grid_table.txt', 
                                                         replicate = x))

####with wlambda.avg
cells.table<-read.table('./output/all_variables_grid_table.txt',sep='\t',header=T,stringsAsFactors = F)
#drop cells where speciation = 0 (there are no species)
cells.table<-cells.table[!is.na(cells.table$mammals.mean.DR)&!is.na(cells.table$birds.mean.DR)&cells.table$mammals.mean.DR>0&cells.table$birds.mean.DR>0,]
#drop cells where past or present elevation is NA or negative
cells.table<-cells.table[!is.na(cells.table$mean.elevation.ETOPO.land)&!is.na(cells.table$mean.elevation.PRISM4.land)&cells.table$mean.elevation.ETOPO.land>0&cells.table$mean.elevation.PRISM4.land>0,]
#drop unnecessary columns
cells.table<-cells.table[,c('cells','mammals.mean.wlambda.avg','birds.mean.wlambda.avg','mean.elevation.ETOPO.land','mean.present.T','elevation.gain','present.minus.past.temperature')]
#get columns where elevation.gain >0
cells.table<-cells.table[cells.table$elevation.gain>0,]
cells.table<-cells.table[complete.cases(cells.table),]
#check correlations among predictors
corrplot(cor(cells.table[,c(2:7)]),method = 'number')
#log transform response variables
cells.table$mammals.mean.wlambda.avg<-log(cells.table$mammals.mean.wlambda.avg)
cells.table$birds.mean.wlambda.avg<-log(cells.table$birds.mean.wlambda.avg)
#log transform other variables with positive values
cells.table$mean.elevation.ETOPO.land<-log(cells.table$mean.elevation.ETOPO.land)
cells.table$elevation.gain<-log(cells.table$elevation.gain)
#scale predictors
cells.table[,c(4:7)]<-apply(cells.table[,c(4:7)],2,function(x)scale(x))

#load grid
grid.world<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
#get coordinates
grid.world.longlat<-lapply(grid.world,function(x) spTransform(x,CRS("+proj=longlat")))
grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) sp::coordinates(x))
grid.coordinates<-do.call("rbind", grid.coordinates)
#get neighbours in 1000km distqnce
neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)

#sem with sarlm for mammals wlambda.avg, elevation, temperature and change elevation and change temperature
sarlm.sem.mammals.wlambda.avg.elevation.gain.temp<-psem(errorsarlm(mammals.mean.wlambda.avg~mean.elevation.ETOPO.land+elevation.gain+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.mammals.wlambda.avg.elevation.gain.temp,standardize = 'none')
sarlm.sem.birds.wlambda.avg.elevation.gain.temp<-psem(errorsarlm(birds.mean.wlambda.avg~mean.elevation.ETOPO.land+elevation.gain+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.birds.wlambda.avg.elevation.gain.temp,standardize = 'none')


####B2_3) sems looking at elevation loss only####

####with wDR
#prepare data
cells.table<-read.table('./output/all_variables_grid_table.txt',sep='\t',header=T,stringsAsFactors = F)
#drop cells where speciation = 0 (there are no species)
cells.table<-cells.table[!is.na(cells.table$mammals.mean.DR)&!is.na(cells.table$birds.mean.DR)&cells.table$mammals.mean.DR>0&cells.table$birds.mean.DR>0,]
#drop cells where past or present elevation is NA or negative
cells.table<-cells.table[!is.na(cells.table$mean.elevation.ETOPO.land)&!is.na(cells.table$mean.elevation.PRISM4.land)&cells.table$mean.elevation.ETOPO.land>0&cells.table$mean.elevation.PRISM4.land>0,]
#drop unnecessary columns
cells.table<-cells.table[,c('cells','mammals.mean.wDR','birds.mean.wDR','mean.elevation.ETOPO.land','mean.present.T','elevation.loss','present.minus.past.temperature')]
#get columns where elevation.loss >0
cells.table<-cells.table[cells.table$elevation.loss>0,]
cells.table<-cells.table[complete.cases(cells.table),]
#check correlations among predictors
corrplot(cor(cells.table[,c(2:7)]),method = 'number')
#log transform response variables
cells.table$mammals.mean.wDR<-log(cells.table$mammals.mean.wDR)
cells.table$birds.mean.wDR<-log(cells.table$birds.mean.wDR)
#log transform other variables with positive values
cells.table$mean.elevation.ETOPO.land<-log(cells.table$mean.elevation.ETOPO.land)
cells.table$elevation.loss<-log(cells.table$elevation.loss)
#scale predictors
cells.table[,c(4:7)]<-apply(cells.table[,c(4:7)],2,function(x)scale(x))

#load grid
grid.world<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
#get coordinates
grid.world.longlat<-lapply(grid.world,function(x) spTransform(x,CRS("+proj=longlat")))
grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) sp::coordinates(x))
grid.coordinates<-do.call("rbind", grid.coordinates)
#get neighbours in 1000km distqnce
neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)

#sem with sarlm for mammals wDR, elevation, temperature and change elevation and change temperature
sarlm.sem.mammals.wDR.elevation.loss.temp<-psem(errorsarlm(mammals.mean.wDR~mean.elevation.ETOPO.land+elevation.loss+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.mammals.wDR.elevation.loss.temp,standardize = 'none')
sarlm.sem.birds.wDR.elevation.loss.temp<-psem(errorsarlm(birds.mean.wDR~mean.elevation.ETOPO.land+elevation.loss+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.birds.wDR.elevation.loss.temp,standardize = 'none')

#for pseudoposteriors: #sem with sarlm for mammals wDR, elevation, temperature and change elevation and change temperature
source('./R/generate_sems.R')
lapply(c(1:100), function(x) sem_loss_spatial_pseudoposterior(variables_table = './output/all_variables_grid_table.txt', 
                                                              replicate = x))

####with wlambda.avg
#prepare data
cells.table<-read.table('./output/all_variables_grid_table.txt',sep='\t',header=T,stringsAsFactors = F)
#drop cells where speciation = 0 (there are no species)
cells.table<-cells.table[!is.na(cells.table$mammals.mean.DR)&!is.na(cells.table$birds.mean.DR)&cells.table$mammals.mean.DR>0&cells.table$birds.mean.DR>0,]
#drop cells where past or present elevation is NA or negative
cells.table<-cells.table[!is.na(cells.table$mean.elevation.ETOPO.land)&!is.na(cells.table$mean.elevation.PRISM4.land)&cells.table$mean.elevation.ETOPO.land>0&cells.table$mean.elevation.PRISM4.land>0,]
#drop unnecessary columns
cells.table<-cells.table[,c('cells','mammals.mean.wlambda.avg','birds.mean.wlambda.avg','mean.elevation.ETOPO.land','mean.present.T','elevation.loss','present.minus.past.temperature')]
#get columns where elevation.loss >0
cells.table<-cells.table[cells.table$elevation.loss>0,]
cells.table<-cells.table[complete.cases(cells.table),]
#check correlations among predictors
corrplot(cor(cells.table[,c(2:7)]),method = 'number')
#log transform response variables
cells.table$mammals.mean.wlambda.avg<-log(cells.table$mammals.mean.wlambda.avg)
cells.table$birds.mean.wlambda.avg<-log(cells.table$birds.mean.wlambda.avg)
#log transform other variables with positive values
cells.table$mean.elevation.ETOPO.land<-log(cells.table$mean.elevation.ETOPO.land)
cells.table$elevation.loss<-log(cells.table$elevation.loss)
#scale predictors
cells.table[,c(4:7)]<-apply(cells.table[,c(4:7)],2,function(x)scale(x))

#load grid
grid.world<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
#get coordinates
grid.world.longlat<-lapply(grid.world,function(x) spTransform(x,CRS("+proj=longlat")))
grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) sp::coordinates(x))
grid.coordinates<-do.call("rbind", grid.coordinates)
#get neighbours in 1000km distqnce
neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)

#sem with sarlm for mammals wlambda.avg, elevation, temperature and change elevation and change temperature
sarlm.sem.mammals.wlambda.avg.elevation.loss.temp<-psem(errorsarlm(mammals.mean.wlambda.avg~mean.elevation.ETOPO.land+elevation.loss+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.mammals.wlambda.avg.elevation.loss.temp,standardize = 'none')
sarlm.sem.birds.wlambda.avg.elevation.loss.temp<-psem(errorsarlm(birds.mean.wlambda.avg~mean.elevation.ETOPO.land+elevation.loss+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.birds.wlambda.avg.elevation.loss.temp,standardize = 'none')

####---C) PLOTS----####

####C2) sem plots#####
####C2_1) plot sems with all cells####
source ('./R/plot_sems.R')

#point estimates (main figures) for mammals
#This is Fig. 2a in the Main Text
sarlm.sem.mammals.wDR.elevation.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.mammals.wDR.elevation.temp,standardize = 'none'),
                                                               mode = 'single')
grViz(sarlm.sem.mammals.wDR.elevation.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_mammals_all_main.pdf")

#This is Fig. S1a in the Main Text
#estimates with CIs (supp figures) for mammals
sarlm.sem.mammals.wDR.elevation.temp.CI.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.mammals.wDR.elevation.temp,standardize = 'none'), 
                                                                mode = 'with_CI')
grViz(sarlm.sem.mammals.wDR.elevation.temp.CI.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_mammals_all_CI.pdf")

#point estimates (main figures) for birds
sarlm.sem.birds.wDR.elevation.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.birds.wDR.elevation.temp,standardize = 'none'),
                                                             mode = 'single')
grViz(sarlm.sem.birds.wDR.elevation.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_birds_all_main.pdf")

#estimates with CIs (supp figures) for birds
sarlm.sem.birds.wDR.elevation.temp.CI.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.birds.wDR.elevation.temp,standardize = 'none'), 
                                                                mode = 'with_CI')
grViz(sarlm.sem.birds.wDR.elevation.temp.CI.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_birds_all_CI.pdf")

####C2_1A) plot sems of pseudoposteriors with all cells####
source('./R/generate_sems.R')
#read in the sems generated in B2_1
sems_all <- lapply(list.files('./output/world/sems/pseudoposterior/global/',
                              '.rds'), function(x) readRDS(paste0('./output/world/sems/pseudoposterior/global/',
                                                                  x)))
#summarise the sems for mammals
posterior_all_sems_mammals_df <- get_CI_sems_pseudoposteriors(sems_object = sems_all,
                                                              taxa = 'mammals')

#plot median estimates + CI across pseudoposterior for mammals
grViz_all_pseudomammals <- coefsdf_to_grViz(coefs.df = posterior_all_sems_mammals_df,mode = 'pseudoposterior')
grViz(grViz_all_pseudomammals)%>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_mammals_all_pseudopos.pdf")

#summarise the sems for birds
posterior_all_sems_birds_df <- get_CI_sems_pseudoposteriors(sems_object = sems_all,
                                                            taxa = 'birds')

#plot median estimates + CI across pseudoposterior for birds
grViz_all_pseudobirds <- coefsdf_to_grViz(coefs.df = posterior_all_sems_birds_df,mode = 'pseudoposterior')
grViz(grViz_all_pseudobirds)%>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_birds_all_pseudopos.pdf")


####C2_2) plot sems with cells with elevation gain only####
source ('./R/plot_sems.R')

#point estimates (main figures) for mammals
sarlm.sem.mammals.wDR.elevation.gain.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.mammals.wDR.elevation.gain.temp,standardize = 'none'),
                                                                    mode = 'single')
grViz(sarlm.sem.mammals.wDR.elevation.gain.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_mammals_gain_main.pdf")

#estimates with CIs (supp figures) for mammals
sarlm.sem.mammals.wDR.elevation.gain.temp.CI.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.mammals.wDR.elevation.gain.temp,standardize = 'none'), 
                                                                       mode = 'with_CI')
grViz(sarlm.sem.mammals.wDR.elevation.gain.temp.CI.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_mammals_gain_CI.pdf")

#point estimates (main figures) for birds
sarlm.sem.birds.wDR.elevation.gain.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.birds.wDR.elevation.gain.temp,standardize = 'none'),
                                                                  mode = 'single')
grViz(sarlm.sem.birds.wDR.elevation.gain.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_birds_gain_main.pdf")

#estimates with CIs (supp figures) for birds
sarlm.sem.birds.wDR.elevation.gain.temp.CI.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.birds.wDR.elevation.gain.temp,standardize = 'none'), 
                                                                     mode = 'with_CI')
grViz(sarlm.sem.birds.wDR.elevation.gain.temp.CI.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_birds_gain_CI.pdf")

####C2_2A) plot seems of pseudoposterior with cells with elevation gain only####
source('./R/generate_sems.R')
#read in the sems generated in B2_2
sems_gain <- lapply(list.files('./output/world/sems/pseudoposterior/gain_elevation/',
                               '.rds'), function(x) readRDS(paste0('./output/world/sems/pseudoposterior/gain_elevation/',
                                                                   x)))
#summarise the sems for mammals
posterior_gain_sems_mammals_df <- get_CI_sems_pseudoposteriors(sems_object = sems_gain,
                                                               taxa = 'mammals')

#plot median estimates + CI across pseudoposterior for mammals
grViz_gain_pseudomammals <- coefsdf_to_grViz(coefs.df = posterior_gain_sems_mammals_df,mode = 'pseudoposterior')
grViz(grViz_gain_pseudomammals)%>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_mammals_gain_pseudopos.pdf")

#summarise the sems for birds
posterior_gain_sems_birds_df <- get_CI_sems_pseudoposteriors(sems_object = sems_gain,
                                                             taxa = 'birds')

#plot median estimates + CI across pseudoposterior for birds
grViz_gain_pseudobirds <- coefsdf_to_grViz(coefs.df = posterior_gain_sems_birds_df,mode = 'pseudoposterior')
grViz(grViz_gain_pseudobirds)%>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_birds_gain_pseudopos.pdf")



####C2_3) plot sems with cells with elevation loss only####
source ('./R/plot_sems.R')

#point estimates (main figures) for mammals
sarlm.sem.mammals.wDR.elevation.loss.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.mammals.wDR.elevation.loss.temp,standardize = 'none'),
                                                                    mode = 'single')
grViz(sarlm.sem.mammals.wDR.elevation.loss.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_mammals_loss_main.pdf")

#estimates with CIs (supp figures) for mammals
sarlm.sem.mammals.wDR.elevation.loss.temp.CI.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.mammals.wDR.elevation.loss.temp,standardize = 'none'), 
                                                                       mode = 'with_CI')
grViz(sarlm.sem.mammals.wDR.elevation.loss.temp.CI.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_mammals_loss_CI.pdf")

#point estimates (main figures) for birds
sarlm.sem.birds.wDR.elevation.loss.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.birds.wDR.elevation.loss.temp,standardize = 'none'),
                                                                  mode = 'single')
grViz(sarlm.sem.birds.wDR.elevation.loss.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_birds_loss_main.pdf")

#estimates with CIs (supp figures) for birds
sarlm.sem.birds.wDR.elevation.loss.temp.CI.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.birds.wDR.elevation.loss.temp,standardize = 'none'), 
                                                                     mode = 'with_CI')
grViz(sarlm.sem.birds.wDR.elevation.loss.temp.CI.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_birds_loss_CI.pdf")
####C2_3A) plot seems of pseudoposterior with cells with elevation loss only####
source('./R/generate_sems.R')
#read in the sems generated in B2_3
sems_loss <- lapply(list.files('./output/world/sems/pseudoposterior/loss_elevation/',
                               '.rds'), function(x) readRDS(paste0('./output/world/sems/pseudoposterior/loss_elevation/',
                                                                   x)))
#summarise the sems for mammals
posterior_loss_sems_mammals_df <- get_CI_sems_pseudoposteriors(sems_object = sems_loss,
                                                               taxa = 'mammals')

#plot median estimates + CI across pseudoposterior for mammals
grViz_loss_pseudomammals <- coefsdf_to_grViz(coefs.df = posterior_loss_sems_mammals_df,mode = 'pseudoposterior')
grViz(grViz_loss_pseudomammals)%>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_mammals_loss_pseudopos.pdf")

#summarise the sems for birds
posterior_loss_sems_birds_df <- get_CI_sems_pseudoposteriors(sems_object = sems_loss,
                                                             taxa = 'birds')

#plot median estimates + CI across pseudoposterior for birds
grViz_loss_pseudobirds <- coefsdf_to_grViz(coefs.df = posterior_loss_sems_birds_df,mode = 'pseudoposterior')
grViz(grViz_loss_pseudobirds)%>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_birds_loss_pseudopos.pdf")



#####SUPPLEMENTARY ANALYSES#####
####TO DO: delete spatial variation of wDR plots - not used at the moment####
####D**) plots of spatial variation of wDR####
source('./R/plots.R')
#maps of wDR for mammals and birds
pdf('./plots/mammalsmeanwDR_gridmap_scale_NEW.pdf',width=11.69,height=8.27)
plot_grid_worldmap_variable_scalebar(table.env.file='./output/all_variables_grid_table.txt',variable='mammals.mean.wDR',ncategories=10,positive.values=TRUE)
dev.off()

pdf('./plots/birdsmeanwDR_gridmap_scale_NEW.pdf',width=11.69,height=8.27)
plot_grid_worldmap_variable_scalebar(table.env.file='./output/all_variables_grid_table.txt',variable='birds.mean.wDR',ncategories=10,positive.values=TRUE)
dev.off()


####D1 analyses with meanDR####
####D1A_1) sems with all cells####

####with DR
#prepare data
cells.table<-read.table('./output/all_variables_grid_table.txt',sep='\t',header=T,stringsAsFactors = F)
#drop cells where speciation = 0 (there are no species)
cells.table<-cells.table[!is.na(cells.table$mammals.mean.DR)&!is.na(cells.table$birds.mean.DR)&cells.table$mammals.mean.DR>0&cells.table$birds.mean.DR>0,]
#drop cells where past or present elevation is NA or negative
cells.table<-cells.table[!is.na(cells.table$mean.elevation.ETOPO.land)&!is.na(cells.table$mean.elevation.PRISM4.land)&cells.table$mean.elevation.ETOPO.land>0&cells.table$mean.elevation.PRISM4.land>0,]
#drop unnecessary columns
cells.table<-cells.table[,c('cells','mammals.mean.DR','birds.mean.DR','mean.elevation.ETOPO.land','mean.elevation.PRISM4.land','mean.present.T','mean.past.T','present.minus.past.elevation','present.minus.past.temperature')]
cells.table<-cells.table[complete.cases(cells.table),]
#check correlations among predictors
corrplot(cor(cells.table[,c(2:9)]),method = 'number')
#log transform response variables
cells.table$mammals.mean.DR<-log(cells.table$mammals.mean.DR)
cells.table$birds.mean.DR<-log(cells.table$birds.mean.DR)
#log transform other variables with positive values
cells.table$mean.elevation.ETOPO.land<-log(cells.table$mean.elevation.ETOPO.land)
cells.table$mean.elevation.PRISM4.land<-log(cells.table$mean.elevation.PRISM4.land)
#scale predictors
cells.table[,c(4:9)]<-apply(cells.table[,c(4:9)],2,function(x)scale(x))

#load grid
grid.world<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
#get coordinates
grid.world.longlat<-lapply(grid.world,function(x) spTransform(x,CRS("+proj=longlat")))
grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) sp::coordinates(x))
grid.coordinates<-do.call("rbind", grid.coordinates)
#get neighbours in 1000km distqnce
neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)

#sem with sarlm for mammals DR, elevation, temperature and change elevation and change temperature
sarlm.sem.mammals.DR.elevation.temp<-psem(errorsarlm(mammals.mean.DR~mean.elevation.ETOPO.land+present.minus.past.elevation+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.mammals.DR.elevation.temp,standardize = 'none')
sarlm.sem.birds.DR.elevation.temp<-psem(errorsarlm(birds.mean.DR~mean.elevation.ETOPO.land+present.minus.past.elevation+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.birds.DR.elevation.temp,standardize = 'none')


####D1A_2) sems looking at elevation gain only####

####with DR
#prepare data
cells.table<-read.table('./output/all_variables_grid_table.txt',sep='\t',header=T,stringsAsFactors = F)
#drop cells where speciation = 0 (there are no species)
cells.table<-cells.table[!is.na(cells.table$mammals.mean.DR)&!is.na(cells.table$birds.mean.DR)&cells.table$mammals.mean.DR>0&cells.table$birds.mean.DR>0,]
#drop cells where past or present elevation is NA or negative
cells.table<-cells.table[!is.na(cells.table$mean.elevation.ETOPO.land)&!is.na(cells.table$mean.elevation.PRISM4.land)&cells.table$mean.elevation.ETOPO.land>0&cells.table$mean.elevation.PRISM4.land>0,]
#drop unnecessary columns
cells.table<-cells.table[,c('cells','mammals.mean.DR','birds.mean.DR','mean.elevation.ETOPO.land','mean.present.T','elevation.gain','present.minus.past.temperature')]
#get columns where elevation.gain >0
cells.table<-cells.table[cells.table$elevation.gain>0,]
cells.table<-cells.table[complete.cases(cells.table),]
#check correlations among predictors
corrplot(cor(cells.table[,c(2:7)]),method = 'number')
#log transform response variables
cells.table$mammals.mean.DR<-log(cells.table$mammals.mean.DR)
cells.table$birds.mean.DR<-log(cells.table$birds.mean.DR)
#log transform other variables with positive values
cells.table$mean.elevation.ETOPO.land<-log(cells.table$mean.elevation.ETOPO.land)
cells.table$elevation.gain<-log(cells.table$elevation.gain)
#scale predictors
cells.table[,c(4:7)]<-apply(cells.table[,c(4:7)],2,function(x)scale(x))

#load grid
grid.world<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
#get coordinates
grid.world.longlat<-lapply(grid.world,function(x) spTransform(x,CRS("+proj=longlat")))
grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) sp::coordinates(x))
grid.coordinates<-do.call("rbind", grid.coordinates)
#get neighbours in 1000km distqnce
neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)

#sem with sarlm for mammals DR, elevation, temperature and change elevation and change temperature
sarlm.sem.mammals.DR.elevation.gain.temp<-psem(errorsarlm(mammals.mean.DR~mean.elevation.ETOPO.land+elevation.gain+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.mammals.DR.elevation.gain.temp,standardize = 'none')
sarlm.sem.birds.DR.elevation.gain.temp<-psem(errorsarlm(birds.mean.DR~mean.elevation.ETOPO.land+elevation.gain+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.birds.DR.elevation.gain.temp,standardize = 'none')

####D1A_3) sems looking at elevation loss only####

####with DR
#prepare data
cells.table<-read.table('./output/all_variables_grid_table.txt',sep='\t',header=T,stringsAsFactors = F)
#drop cells where speciation = 0 (there are no species)
cells.table<-cells.table[!is.na(cells.table$mammals.mean.DR)&!is.na(cells.table$birds.mean.DR)&cells.table$mammals.mean.DR>0&cells.table$birds.mean.DR>0,]
#drop cells where past or present elevation is NA or negative
cells.table<-cells.table[!is.na(cells.table$mean.elevation.ETOPO.land)&!is.na(cells.table$mean.elevation.PRISM4.land)&cells.table$mean.elevation.ETOPO.land>0&cells.table$mean.elevation.PRISM4.land>0,]
#drop unnecessary columns
cells.table<-cells.table[,c('cells','mammals.mean.DR','birds.mean.DR','mean.elevation.ETOPO.land','mean.present.T','elevation.loss','present.minus.past.temperature')]
#get columns where elevation.loss >0
cells.table<-cells.table[cells.table$elevation.loss>0,]
cells.table<-cells.table[complete.cases(cells.table),]
#check correlations among predictors
corrplot(cor(cells.table[,c(2:7)]),method = 'number')
#log transform response variables
cells.table$mammals.mean.DR<-log(cells.table$mammals.mean.DR)
cells.table$birds.mean.DR<-log(cells.table$birds.mean.DR)
#log transform other variables with positive values
cells.table$mean.elevation.ETOPO.land<-log(cells.table$mean.elevation.ETOPO.land)
cells.table$elevation.loss<-log(cells.table$elevation.loss)
#scale predictors
cells.table[,c(4:7)]<-apply(cells.table[,c(4:7)],2,function(x)scale(x))

#load grid
grid.world<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
#get coordinates
grid.world.longlat<-lapply(grid.world,function(x) spTransform(x,CRS("+proj=longlat")))
grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) sp::coordinates(x))
grid.coordinates<-do.call("rbind", grid.coordinates)
#get neighbours in 1000km distqnce
neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)

#sem with sarlm for mammals DR, elevation, temperature and change elevation and change temperature
sarlm.sem.mammals.DR.elevation.loss.temp<-psem(errorsarlm(mammals.mean.DR~mean.elevation.ETOPO.land+elevation.loss+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.mammals.DR.elevation.loss.temp,standardize = 'none')
sarlm.sem.birds.DR.elevation.loss.temp<-psem(errorsarlm(birds.mean.DR~mean.elevation.ETOPO.land+elevation.loss+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.birds.DR.elevation.loss.temp,standardize = 'none')

####D1B_sem plots#####
####D1B_1) plot sems with all cells####
source ('./R/plot_sems.R')

#point estimates (main figures) for mammals
sarlm.sem.mammals.DR.elevation.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.mammals.DR.elevation.temp,standardize = 'none'),
                                                              mode = 'single')
grViz(sarlm.sem.mammals.DR.elevation.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_mammals_all_main_meanDR.pdf")

#point estimates (main figures) for birds
sarlm.sem.birds.DR.elevation.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.birds.DR.elevation.temp,standardize = 'none'),
                                                            mode = 'single')
grViz(sarlm.sem.birds.DR.elevation.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_birds_all_main_meanDR.pdf")

####D1B_2) plot sems with cells with elevation gain only####
source ('./R/plot_sems.R')

#point estimates (main figures) for mammals
sarlm.sem.mammals.DR.elevation.gain.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.mammals.DR.elevation.gain.temp,standardize = 'none'),
                                                                   mode = 'single')
grViz(sarlm.sem.mammals.DR.elevation.gain.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_mammals_gain_main_meanDR.pdf")

#point estimates (main figures) for birds
sarlm.sem.birds.DR.elevation.gain.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.birds.DR.elevation.gain.temp,standardize = 'none'),
                                                                 mode = 'single')
grViz(sarlm.sem.birds.DR.elevation.gain.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_birds_gain_main_meanDR.pdf")

####D1B_3) plot sems with cells with elevation loss only####
source ('./R/plot_sems.R')

#point estimates (main figures) for mammals
sarlm.sem.mammals.DR.elevation.loss.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.mammals.DR.elevation.loss.temp,standardize = 'none'),
                                                                   mode = 'single')
grViz(sarlm.sem.mammals.DR.elevation.loss.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_mammals_loss_main_meanDR.pdf")

#point estimates (main figures) for birds
sarlm.sem.birds.DR.elevation.loss.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.birds.DR.elevation.loss.temp,standardize = 'none'),
                                                                 mode = 'single')
grViz(sarlm.sem.birds.DR.elevation.loss.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_birds_loss_main_meanDR.pdf")


#####D2 analyses with weighted BAMM lambda avg####
####D2A_1) sems with all cells####
#prepare data
cells.table<-read.table('./output/all_variables_grid_table.txt',sep='\t',header=T,stringsAsFactors = F)
#drop cells where speciation = 0 (there are no species)
cells.table<-cells.table[!is.na(cells.table$mammals.mean.DR)&!is.na(cells.table$birds.mean.DR)&cells.table$mammals.mean.DR>0&cells.table$birds.mean.DR>0,]
#drop cells where past or present elevation is NA or negative
cells.table<-cells.table[!is.na(cells.table$mean.elevation.ETOPO.land)&!is.na(cells.table$mean.elevation.PRISM4.land)&cells.table$mean.elevation.ETOPO.land>0&cells.table$mean.elevation.PRISM4.land>0,]
#drop unnecessary columns
cells.table<-cells.table[,c('cells','mammals.mean.wlambda.avg','birds.mean.wlambda.avg','mean.elevation.ETOPO.land','mean.elevation.PRISM4.land','mean.present.T','mean.past.T','present.minus.past.elevation','present.minus.past.temperature')]
cells.table<-cells.table[complete.cases(cells.table),]
#check correlations among predictors
corrplot(cor(cells.table[,c(2:9)]),method = 'number')
#log transform response variables
cells.table$mammals.mean.wlambda.avg<-log(cells.table$mammals.mean.wlambda.avg)
cells.table$birds.mean.wlambda.avg<-log(cells.table$birds.mean.wlambda.avg)
#log transform other variables with positive values
cells.table$mean.elevation.ETOPO.land<-log(cells.table$mean.elevation.ETOPO.land)
cells.table$mean.elevation.PRISM4.land<-log(cells.table$mean.elevation.PRISM4.land)
#scale predictors
cells.table[,c(4:9)]<-apply(cells.table[,c(4:9)],2,function(x)scale(x))

#load grid
grid.world<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
#get coordinates
grid.world.longlat<-lapply(grid.world,function(x) spTransform(x,CRS("+proj=longlat")))
grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) sp::coordinates(x))
grid.coordinates<-do.call("rbind", grid.coordinates)
#get neighbours in 1000km distqnce
neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)

#sem with sarlm for mammals wlambda.avg, elevation, temperature and change elevation and change temperature
sarlm.sem.mammals.wlambda.avg.elevation.temp<-psem(errorsarlm(mammals.mean.wlambda.avg~mean.elevation.ETOPO.land+present.minus.past.elevation+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.mammals.wlambda.avg.elevation.temp,standardize = 'none')
sarlm.sem.birds.wlambda.avg.elevation.temp<-psem(errorsarlm(birds.mean.wlambda.avg~mean.elevation.ETOPO.land+present.minus.past.elevation+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.birds.wlambda.avg.elevation.temp,standardize = 'none')

####D2A_2) sems looking at elevation gain only####

#with BAMM lambda.avg
cells.table<-read.table('./output/all_variables_grid_table.txt',sep='\t',header=T,stringsAsFactors = F)
#drop cells where speciation = 0 (there are no species)
cells.table<-cells.table[!is.na(cells.table$mammals.mean.DR)&!is.na(cells.table$birds.mean.DR)&cells.table$mammals.mean.DR>0&cells.table$birds.mean.DR>0,]
#drop cells where past or present elevation is NA or negative
cells.table<-cells.table[!is.na(cells.table$mean.elevation.ETOPO.land)&!is.na(cells.table$mean.elevation.PRISM4.land)&cells.table$mean.elevation.ETOPO.land>0&cells.table$mean.elevation.PRISM4.land>0,]
#drop unnecessary columns
cells.table<-cells.table[,c('cells','mammals.mean.wlambda.avg','birds.mean.wlambda.avg','mean.elevation.ETOPO.land','mean.present.T','elevation.gain','present.minus.past.temperature')]
#get columns where elevation.gain >0
cells.table<-cells.table[cells.table$elevation.gain>0,]
cells.table<-cells.table[complete.cases(cells.table),]
#check correlations among predictors
corrplot(cor(cells.table[,c(2:7)]),method = 'number')
#log transform response variables
cells.table$mammals.mean.wlambda.avg<-log(cells.table$mammals.mean.wlambda.avg)
cells.table$birds.mean.wlambda.avg<-log(cells.table$birds.mean.wlambda.avg)
#log transform other variables with positive values
cells.table$mean.elevation.ETOPO.land<-log(cells.table$mean.elevation.ETOPO.land)
cells.table$elevation.gain<-log(cells.table$elevation.gain)
#scale predictors
cells.table[,c(4:7)]<-apply(cells.table[,c(4:7)],2,function(x)scale(x))

#load grid
grid.world<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
#get coordinates
grid.world.longlat<-lapply(grid.world,function(x) spTransform(x,CRS("+proj=longlat")))
grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) sp::coordinates(x))
grid.coordinates<-do.call("rbind", grid.coordinates)
#get neighbours in 1000km distqnce
neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)

#sem with sarlm for mammals wlambda.avg, elevation, temperature and change elevation and change temperature
sarlm.sem.mammals.wlambda.avg.elevation.gain.temp<-psem(errorsarlm(mammals.mean.wlambda.avg~mean.elevation.ETOPO.land+elevation.gain+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.mammals.wlambda.avg.elevation.gain.temp,standardize = 'none')
sarlm.sem.birds.wlambda.avg.elevation.gain.temp<-psem(errorsarlm(birds.mean.wlambda.avg~mean.elevation.ETOPO.land+elevation.gain+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.birds.wlambda.avg.elevation.gain.temp,standardize = 'none')


####D2A_3) sems looking at elevation loss only####

####with wlambda.avg
#prepare data
cells.table<-read.table('./output/all_variables_grid_table.txt',sep='\t',header=T,stringsAsFactors = F)
#drop cells where speciation = 0 (there are no species)
cells.table<-cells.table[!is.na(cells.table$mammals.mean.DR)&!is.na(cells.table$birds.mean.DR)&cells.table$mammals.mean.DR>0&cells.table$birds.mean.DR>0,]
#drop cells where past or present elevation is NA or negative
cells.table<-cells.table[!is.na(cells.table$mean.elevation.ETOPO.land)&!is.na(cells.table$mean.elevation.PRISM4.land)&cells.table$mean.elevation.ETOPO.land>0&cells.table$mean.elevation.PRISM4.land>0,]
#drop unnecessary columns
cells.table<-cells.table[,c('cells','mammals.mean.wlambda.avg','birds.mean.wlambda.avg','mean.elevation.ETOPO.land','mean.present.T','elevation.loss','present.minus.past.temperature')]
#get columns where elevation.loss >0
cells.table<-cells.table[cells.table$elevation.loss>0,]
cells.table<-cells.table[complete.cases(cells.table),]
#check correlations among predictors
corrplot(cor(cells.table[,c(2:7)]),method = 'number')
#log transform response variables
cells.table$mammals.mean.wlambda.avg<-log(cells.table$mammals.mean.wlambda.avg)
cells.table$birds.mean.wlambda.avg<-log(cells.table$birds.mean.wlambda.avg)
#log transform other variables with positive values
cells.table$mean.elevation.ETOPO.land<-log(cells.table$mean.elevation.ETOPO.land)
cells.table$elevation.loss<-log(cells.table$elevation.loss)
#scale predictors
cells.table[,c(4:7)]<-apply(cells.table[,c(4:7)],2,function(x)scale(x))

#load grid
grid.world<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
#get coordinates
grid.world.longlat<-lapply(grid.world,function(x) spTransform(x,CRS("+proj=longlat")))
grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) sp::coordinates(x))
grid.coordinates<-do.call("rbind", grid.coordinates)
#get neighbours in 1000km distqnce
neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)

#sem with sarlm for mammals wlambda.avg, elevation, temperature and change elevation and change temperature
sarlm.sem.mammals.wlambda.avg.elevation.loss.temp<-psem(errorsarlm(mammals.mean.wlambda.avg~mean.elevation.ETOPO.land+elevation.loss+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.mammals.wlambda.avg.elevation.loss.temp,standardize = 'none')
sarlm.sem.birds.wlambda.avg.elevation.loss.temp<-psem(errorsarlm(birds.mean.wlambda.avg~mean.elevation.ETOPO.land+elevation.loss+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.birds.wlambda.avg.elevation.loss.temp,standardize = 'none')


####D2B_sem plots####
####D2B_1) plot sems with all cells####
source ('./R/plot_sems.R')

#point estimates (main figures) for mammals
sarlm.sem.mammals.wlambda.avg.elevation.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.mammals.wlambda.avg.elevation.temp,standardize = 'none'),
                                                                       mode = 'single')
grViz(sarlm.sem.mammals.wlambda.avg.elevation.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_mammals_lambdaavg_all_main.pdf")

#point estimates (main figures) for birds
sarlm.sem.birds.wlambda.avg.elevation.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.birds.wlambda.avg.elevation.temp,standardize = 'none'),
                                                                     mode = 'single')
grViz(sarlm.sem.birds.wlambda.avg.elevation.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_birds_lambdaavg_all_main.pdf")

####D2B_2) plot sems with cells with elevation gain only####
source ('./R/plot_sems.R')

#point estimates (main figures) for mammals
sarlm.sem.mammals.wlambda.avg.elevation.gain.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.mammals.wlambda.avg.elevation.gain.temp,standardize = 'none'),
                                                                            mode = 'single')
grViz(sarlm.sem.mammals.wlambda.avg.elevation.gain.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_mammals_gain_main_wlambda.avg.pdf")

#point estimates (main figures) for birds
sarlm.sem.birds.wlambda.avg.elevation.gain.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.birds.wlambda.avg.elevation.gain.temp,standardize = 'none'),
                                                                          mode = 'single')
grViz(sarlm.sem.birds.wlambda.avg.elevation.gain.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_birds_gain_main_wlambda.avg.pdf")

####D2B_3) plot sems with cells with elevation loss only####
source ('./R/plot_sems.R')

#point estimates (main figures) for mammals
sarlm.sem.mammals.wlambda.avg.elevation.loss.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.mammals.wlambda.avg.elevation.loss.temp,standardize = 'none'),
                                                                            mode = 'single')
grViz(sarlm.sem.mammals.wlambda.avg.elevation.loss.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_mammals_loss_main_wlambda.avg.pdf")

#point estimates (main figures) for birds
sarlm.sem.birds.wlambda.avg.elevation.loss.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.birds.wlambda.avg.elevation.loss.temp,standardize = 'none'),
                                                                          mode = 'single')
grViz(sarlm.sem.birds.wlambda.avg.elevation.loss.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_birds_loss_main_wlambda.avg.pdf")

#####D3 analyses with weighted geometric mean BAMM lambda avg####
####D3A_1) sems with all cells####
#prepare data
cells.table<-read.table('./output/all_variables_grid_table.txt',sep='\t',header=T,stringsAsFactors = F)
#drop cells where speciation = 0 (there are no species)
cells.table<-cells.table[!is.na(cells.table$mammals.mean.DR)&!is.na(cells.table$birds.mean.DR)&cells.table$mammals.mean.DR>0&cells.table$birds.mean.DR>0,]
#drop cells where past or present elevation is NA or negative
cells.table<-cells.table[!is.na(cells.table$mean.elevation.ETOPO.land)&!is.na(cells.table$mean.elevation.PRISM4.land)&cells.table$mean.elevation.ETOPO.land>0&cells.table$mean.elevation.PRISM4.land>0,]
#drop unnecessary columns
cells.table<-cells.table[,c('cells','mammals.geomean.wlambda.avg','birds.geomean.wlambda.avg','mean.elevation.ETOPO.land','mean.elevation.PRISM4.land','mean.present.T','mean.past.T','present.minus.past.elevation','present.minus.past.temperature')]
cells.table<-cells.table[complete.cases(cells.table),]
#check correlations among predictors
corrplot(cor(cells.table[,c(2:9)]),method = 'number')
#log transform response variables
cells.table$mammals.geomean.wlambda.avg<-log(cells.table$mammals.geomean.wlambda.avg)
cells.table$birds.geomean.wlambda.avg<-log(cells.table$birds.geomean.wlambda.avg)
#log transform other variables with positive values
cells.table$mean.elevation.ETOPO.land<-log(cells.table$mean.elevation.ETOPO.land)
cells.table$mean.elevation.PRISM4.land<-log(cells.table$mean.elevation.PRISM4.land)
#scale predictors
cells.table[,c(4:9)]<-apply(cells.table[,c(4:9)],2,function(x)scale(x))

#load grid
grid.world<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
#get coordinates
grid.world.longlat<-lapply(grid.world,function(x) spTransform(x,CRS("+proj=longlat")))
grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) sp::coordinates(x))
grid.coordinates<-do.call("rbind", grid.coordinates)
#get neighbours in 1000km distqnce
neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)

#sem with sarlm for mammals wlambda.avg, elevation, temperature and change elevation and change temperature
sarlm.sem.mammals.geomean.wlambda.avg.elevation.temp<-psem(errorsarlm(mammals.geomean.wlambda.avg~mean.elevation.ETOPO.land+present.minus.past.elevation+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.mammals.geomean.wlambda.avg.elevation.temp,standardize = 'none')
sarlm.sem.birds.geomean.wlambda.avg.elevation.temp<-psem(errorsarlm(birds.geomean.wlambda.avg~mean.elevation.ETOPO.land+present.minus.past.elevation+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.birds.geomean.wlambda.avg.elevation.temp,standardize = 'none')

####D3A_2) sems looking at elevation gain only####

#with BAMM lambda.avg
cells.table<-read.table('./output/all_variables_grid_table.txt',sep='\t',header=T,stringsAsFactors = F)
#drop cells where speciation = 0 (there are no species)
cells.table<-cells.table[!is.na(cells.table$mammals.mean.DR)&!is.na(cells.table$birds.mean.DR)&cells.table$mammals.mean.DR>0&cells.table$birds.mean.DR>0,]
#drop cells where past or present elevation is NA or negative
cells.table<-cells.table[!is.na(cells.table$mean.elevation.ETOPO.land)&!is.na(cells.table$mean.elevation.PRISM4.land)&cells.table$mean.elevation.ETOPO.land>0&cells.table$mean.elevation.PRISM4.land>0,]
#drop unnecessary columns
cells.table<-cells.table[,c('cells','mammals.geomean.wlambda.avg','birds.geomean.wlambda.avg','mean.elevation.ETOPO.land','mean.present.T','elevation.gain','present.minus.past.temperature')]
#get columns where elevation.gain >0
cells.table<-cells.table[cells.table$elevation.gain>0,]
cells.table<-cells.table[complete.cases(cells.table),]
#check correlations among predictors
corrplot(cor(cells.table[,c(2:7)]),method = 'number')
#log transform response variables
cells.table$mammals.geomean.wlambda.avg<-log(cells.table$mammals.geomean.wlambda.avg)
cells.table$birds.geomean.wlambda.avg<-log(cells.table$birds.geomean.wlambda.avg)
#log transform other variables with positive values
cells.table$mean.elevation.ETOPO.land<-log(cells.table$mean.elevation.ETOPO.land)
cells.table$elevation.gain<-log(cells.table$elevation.gain)
#scale predictors
cells.table[,c(4:7)]<-apply(cells.table[,c(4:7)],2,function(x)scale(x))

#load grid
grid.world<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
#get coordinates
grid.world.longlat<-lapply(grid.world,function(x) spTransform(x,CRS("+proj=longlat")))
grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) sp::coordinates(x))
grid.coordinates<-do.call("rbind", grid.coordinates)
#get neighbours in 1000km distqnce
neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)

#sem with sarlm for mammals wlambda.avg, elevation, temperature and change elevation and change temperature
sarlm.sem.mammals.geomean.wlambda.avg.elevation.gain.temp<-psem(errorsarlm(mammals.geomean.wlambda.avg~mean.elevation.ETOPO.land+elevation.gain+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.mammals.geomean.wlambda.avg.elevation.gain.temp,standardize = 'none')
sarlm.sem.birds.geomean.wlambda.avg.elevation.gain.temp<-psem(errorsarlm(birds.geomean.wlambda.avg~mean.elevation.ETOPO.land+elevation.gain+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.birds.geomean.wlambda.avg.elevation.gain.temp,standardize = 'none')


####D3A_3) sems looking at elevation loss only####

####with geomean.wlambda.avg
#prepare data
cells.table<-read.table('./output/all_variables_grid_table.txt',sep='\t',header=T,stringsAsFactors = F)
#drop cells where speciation = 0 (there are no species)
cells.table<-cells.table[!is.na(cells.table$mammals.mean.DR)&!is.na(cells.table$birds.mean.DR)&cells.table$mammals.mean.DR>0&cells.table$birds.mean.DR>0,]
#drop cells where past or present elevation is NA or negative
cells.table<-cells.table[!is.na(cells.table$mean.elevation.ETOPO.land)&!is.na(cells.table$mean.elevation.PRISM4.land)&cells.table$mean.elevation.ETOPO.land>0&cells.table$mean.elevation.PRISM4.land>0,]
#drop unnecessary columns
cells.table<-cells.table[,c('cells','mammals.geomean.wlambda.avg','birds.geomean.wlambda.avg','mean.elevation.ETOPO.land','mean.present.T','elevation.loss','present.minus.past.temperature')]
#get columns where elevation.loss >0
cells.table<-cells.table[cells.table$elevation.loss>0,]
cells.table<-cells.table[complete.cases(cells.table),]
#check correlations among predictors
corrplot(cor(cells.table[,c(2:7)]),method = 'number')
#log transform response variables
cells.table$mammals.geomean.wlambda.avg<-log(cells.table$mammals.geomean.wlambda.avg)
cells.table$birds.geomean.wlambda.avg<-log(cells.table$birds.geomean.wlambda.avg)
#log transform other variables with positive values
cells.table$mean.elevation.ETOPO.land<-log(cells.table$mean.elevation.ETOPO.land)
cells.table$elevation.loss<-log(cells.table$elevation.loss)
#scale predictors
cells.table[,c(4:7)]<-apply(cells.table[,c(4:7)],2,function(x)scale(x))

#load grid
grid.world<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
#get coordinates
grid.world.longlat<-lapply(grid.world,function(x) spTransform(x,CRS("+proj=longlat")))
grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) sp::coordinates(x))
grid.coordinates<-do.call("rbind", grid.coordinates)
#get neighbours in 1000km distqnce
neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)

#sem with sarlm for mammals wlambda.avg, elevation, temperature and change elevation and change temperature
sarlm.sem.mammals.geomean.wlambda.avg.elevation.loss.temp<-psem(errorsarlm(mammals.geomean.wlambda.avg~mean.elevation.ETOPO.land+elevation.loss+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.mammals.geomean.wlambda.avg.elevation.loss.temp,standardize = 'none')
sarlm.sem.birds.geomean.wlambda.avg.elevation.loss.temp<-psem(errorsarlm(birds.geomean.wlambda.avg~mean.elevation.ETOPO.land+elevation.loss+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.birds.geomean.wlambda.avg.elevation.loss.temp,standardize = 'none')


####D3B_sem plots####
####D3B_1) plot sems with all cells####
source ('./R/plot_sems.R')

#point estimates (main figures) for mammals
sarlm.sem.mammals.geomean.wlambda.avg.elevation.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.mammals.geomean.wlambda.avg.elevation.temp,standardize = 'none'),
                                                                               mode = 'single')
grViz(sarlm.sem.mammals.geomean.wlambda.avg.elevation.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_mammals_lambdaavg_all_main.pdf")

#point estimates (main figures) for birds
sarlm.sem.birds.geomean.wlambda.avg.elevation.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.birds.geomean.wlambda.avg.elevation.temp,standardize = 'none'),
                                                                             mode = 'single')
grViz(sarlm.sem.birds.geomean.wlambda.avg.elevation.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_birds_lambdaavg_all_main.pdf")

####D3B_2) plot sems with cells with elevation gain only####
source ('./R/plot_sems.R')

#point estimates (main figures) for mammals
sarlm.sem.mammals.geomean.wlambda.avg.elevation.gain.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.mammals.geomean.wlambda.avg.elevation.gain.temp,standardize = 'none'),
                                                                                    mode = 'single')
grViz(sarlm.sem.mammals.geomean.wlambda.avg.elevation.gain.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_mammals_gain_main_geomean.wlambda.avg.pdf")

#point estimates (main figures) for birds
sarlm.sem.birds.geomean.wlambda.avg.elevation.gain.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.birds.geomean.wlambda.avg.elevation.gain.temp,standardize = 'none'),
                                                                                  mode = 'single')
grViz(sarlm.sem.birds.geomean.wlambda.avg.elevation.gain.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_birds_gain_main_geomean.wlambda.avg.pdf")

####DB_3) plot sems with cells with elevation loss only####
source ('./R/plot_sems.R')

#point estimates (main figures) for mammals
sarlm.sem.mammals.geomean.wlambda.avg.elevation.loss.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.mammals.geomean.wlambda.avg.elevation.loss.temp,standardize = 'none'),
                                                                                    mode = 'single')
grViz(sarlm.sem.mammals.geomean.wlambda.avg.elevation.loss.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_mammals_loss_main_geomean.wlambda.avg.pdf")

#point estimates (main figures) for birds
sarlm.sem.birds.geomean.wlambda.avg.elevation.loss.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.birds.geomean.wlambda.avg.elevation.loss.temp,standardize = 'none'),
                                                                                  mode = 'single')
grViz(sarlm.sem.birds.geomean.wlambda.avg.elevation.loss.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_birds_loss_main_geomean.wlambda.avg.pdf")

#####D4 analyses with weighted geometric mean wDR####
####D4A_1) sems with all cells####
#prepare data
cells.table<-read.table('./output/all_variables_grid_table.txt',sep='\t',header=T,stringsAsFactors = F)
#drop cells where speciation = 0 (there are no species)
cells.table<-cells.table[!is.na(cells.table$mammals.mean.DR)&!is.na(cells.table$birds.mean.DR)&cells.table$mammals.mean.DR>0&cells.table$birds.mean.DR>0,]
#drop cells where past or present elevation is NA or negative
cells.table<-cells.table[!is.na(cells.table$mean.elevation.ETOPO.land)&!is.na(cells.table$mean.elevation.PRISM4.land)&cells.table$mean.elevation.ETOPO.land>0&cells.table$mean.elevation.PRISM4.land>0,]
#drop unnecessary columns
cells.table<-cells.table[,c('cells','mammals.geomean.wDR','birds.geomean.wDR','mean.elevation.ETOPO.land','mean.elevation.PRISM4.land','mean.present.T','mean.past.T','present.minus.past.elevation','present.minus.past.temperature')]
cells.table<-cells.table[complete.cases(cells.table),]
#check correlations among predictors
corrplot(cor(cells.table[,c(2:9)]),method = 'number')
#log transform response variables
cells.table$mammals.geomean.wDR<-log(cells.table$mammals.geomean.wDR)
cells.table$birds.geomean.wDR<-log(cells.table$birds.geomean.wDR)
#log transform other variables with positive values
cells.table$mean.elevation.ETOPO.land<-log(cells.table$mean.elevation.ETOPO.land)
cells.table$mean.elevation.PRISM4.land<-log(cells.table$mean.elevation.PRISM4.land)
#scale predictors
cells.table[,c(4:9)]<-apply(cells.table[,c(4:9)],2,function(x)scale(x))

#load grid
grid.world<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
#get coordinates
grid.world.longlat<-lapply(grid.world,function(x) spTransform(x,CRS("+proj=longlat")))
grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) sp::coordinates(x))
grid.coordinates<-do.call("rbind", grid.coordinates)
#get neighbours in 1000km distqnce
neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)

#sem with sarlm for mammals geomeanwDR, elevation, temperature and change elevation and change temperature
sarlm.sem.mammals.geomean.wDR.elevation.temp<-psem(errorsarlm(mammals.geomean.wDR~mean.elevation.ETOPO.land+present.minus.past.elevation+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.mammals.geomean.wDR.elevation.temp,standardize = 'none')
sarlm.sem.birds.geomean.wDR.elevation.temp<-psem(errorsarlm(birds.geomean.wDR~mean.elevation.ETOPO.land+present.minus.past.elevation+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~present.minus.past.elevation,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.birds.geomean.wDR.elevation.temp,standardize = 'none')

####D4A_2) sems looking at elevation gain only####

#with geomean wDR
cells.table<-read.table('./output/all_variables_grid_table.txt',sep='\t',header=T,stringsAsFactors = F)
#drop cells where speciation = 0 (there are no species)
cells.table<-cells.table[!is.na(cells.table$mammals.mean.DR)&!is.na(cells.table$birds.mean.DR)&cells.table$mammals.mean.DR>0&cells.table$birds.mean.DR>0,]
#drop cells where past or present elevation is NA or negative
cells.table<-cells.table[!is.na(cells.table$mean.elevation.ETOPO.land)&!is.na(cells.table$mean.elevation.PRISM4.land)&cells.table$mean.elevation.ETOPO.land>0&cells.table$mean.elevation.PRISM4.land>0,]
#drop unnecessary columns
cells.table<-cells.table[,c('cells','mammals.geomean.wDR','birds.geomean.wDR','mean.elevation.ETOPO.land','mean.present.T','elevation.gain','present.minus.past.temperature')]
#get columns where elevation.gain >0
cells.table<-cells.table[cells.table$elevation.gain>0,]
cells.table<-cells.table[complete.cases(cells.table),]
#check correlations among predictors
corrplot(cor(cells.table[,c(2:7)]),method = 'number')
#log transform response variables
cells.table$mammals.geomean.wDR<-log(cells.table$mammals.geomean.wDR)
cells.table$birds.geomean.wDR<-log(cells.table$birds.geomean.wDR)
#log transform other variables with positive values
cells.table$mean.elevation.ETOPO.land<-log(cells.table$mean.elevation.ETOPO.land)
cells.table$elevation.gain<-log(cells.table$elevation.gain)
#scale predictors
cells.table[,c(4:7)]<-apply(cells.table[,c(4:7)],2,function(x)scale(x))

#load grid
grid.world<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
#get coordinates
grid.world.longlat<-lapply(grid.world,function(x) spTransform(x,CRS("+proj=longlat")))
grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) sp::coordinates(x))
grid.coordinates<-do.call("rbind", grid.coordinates)
#get neighbours in 1000km distqnce
neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)

#sem with sarlm for mammals geomeanwDR, elevation, temperature and change elevation and change temperature
sarlm.sem.mammals.geomean.wDR.elevation.gain.temp<-psem(errorsarlm(mammals.geomean.wDR~mean.elevation.ETOPO.land+elevation.gain+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.mammals.geomean.wDR.elevation.gain.temp,standardize = 'none')
sarlm.sem.birds.geomean.wDR.elevation.gain.temp<-psem(errorsarlm(birds.geomean.wDR~mean.elevation.ETOPO.land+elevation.gain+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.gain,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.birds.geomean.wDR.elevation.gain.temp,standardize = 'none')


####D4A_3) sems looking at elevation loss only####

####with geomean.wDR
#prepare data
cells.table<-read.table('./output/all_variables_grid_table.txt',sep='\t',header=T,stringsAsFactors = F)
#drop cells where speciation = 0 (there are no species)
cells.table<-cells.table[!is.na(cells.table$mammals.mean.DR)&!is.na(cells.table$birds.mean.DR)&cells.table$mammals.mean.DR>0&cells.table$birds.mean.DR>0,]
#drop cells where past or present elevation is NA or negative
cells.table<-cells.table[!is.na(cells.table$mean.elevation.ETOPO.land)&!is.na(cells.table$mean.elevation.PRISM4.land)&cells.table$mean.elevation.ETOPO.land>0&cells.table$mean.elevation.PRISM4.land>0,]
#drop unnecessary columns
cells.table<-cells.table[,c('cells','mammals.geomean.wDR','birds.geomean.wDR','mean.elevation.ETOPO.land','mean.present.T','elevation.loss','present.minus.past.temperature')]
#get columns where elevation.loss >0
cells.table<-cells.table[cells.table$elevation.loss>0,]
cells.table<-cells.table[complete.cases(cells.table),]
#check correlations among predictors
corrplot(cor(cells.table[,c(2:7)]),method = 'number')
#log transform response variables
cells.table$mammals.geomean.wDR<-log(cells.table$mammals.geomean.wDR)
cells.table$birds.geomean.wDR<-log(cells.table$birds.geomean.wDR)
#log transform other variables with positive values
cells.table$mean.elevation.ETOPO.land<-log(cells.table$mean.elevation.ETOPO.land)
cells.table$elevation.loss<-log(cells.table$elevation.loss)
#scale predictors
cells.table[,c(4:7)]<-apply(cells.table[,c(4:7)],2,function(x)scale(x))

#load grid
grid.world<-readRDS('./raw_data/grid_World_RealmsMerged_100.rds')
#get coordinates
grid.world.longlat<-lapply(grid.world,function(x) spTransform(x,CRS("+proj=longlat")))
grid.coordinates<-lapply(grid.world[cells.table$cells],function(x) sp::coordinates(x))
grid.coordinates<-do.call("rbind", grid.coordinates)
#get neighbours in 1000km distqnce
neighbours.1000<-dnearneigh(grid.coordinates,d1=0,d2=1000)
neighbours.1000.w<-nb2listw(neighbours.1000,style="W",zero.policy = TRUE)

#sem with sarlm for mammals geomeanwDR, elevation, temperature and change elevation and change temperature
sarlm.sem.mammals.geomean.wDR.elevation.loss.temp<-psem(errorsarlm(mammals.geomean.wDR~mean.elevation.ETOPO.land+elevation.loss+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.mammals.geomean.wDR.elevation.loss.temp,standardize = 'none')
sarlm.sem.birds.geomean.wDR.elevation.loss.temp<-psem(errorsarlm(birds.geomean.wDR~mean.elevation.ETOPO.land+elevation.loss+mean.present.T+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.elevation.ETOPO.land~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(mean.present.T~mean.elevation.ETOPO.land+present.minus.past.temperature,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'),errorsarlm(present.minus.past.temperature~elevation.loss,data=cells.table,listw = neighbours.1000.w,zero.policy = TRUE,quiet=FALSE,method='spam'))
coefs(sarlm.sem.birds.geomean.wDR.elevation.loss.temp,standardize = 'none')


####D4B_sem plots####
####D4B_1) plot sems with all cells####
source ('./R/plot_sems.R')

#point estimates (main figures) for mammals
sarlm.sem.mammals.geomean.wDR.elevation.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.mammals.geomean.wDR.elevation.temp,standardize = 'none'),
                                                                       mode = 'single')
grViz(sarlm.sem.mammals.geomean.wDR.elevation.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_mammals_geomeanwDR_all_main.pdf")

#point estimates (main figures) for birds
sarlm.sem.birds.geomean.wDR.elevation.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.birds.geomean.wDR.elevation.temp,standardize = 'none'),
                                                                     mode = 'single')
grViz(sarlm.sem.birds.geomean.wDR.elevation.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_birds_geomeanwDR_all_main.pdf")

####D4B_2) plot sems with cells with elevation gain only####
source ('./R/plot_sems.R')

#point estimates (main figures) for mammals
sarlm.sem.mammals.geomean.wDR.elevation.gain.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.mammals.geomean.wDR.elevation.gain.temp,standardize = 'none'),
                                                                            mode = 'single')
grViz(sarlm.sem.mammals.geomean.wDR.elevation.gain.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_mammals_gain_main_geomeanwDR.pdf")

#point estimates (main figures) for birds
sarlm.sem.birds.geomean.wDR.elevation.gain.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.birds.geomean.wDR.elevation.gain.temp,standardize = 'none'),
                                                                          mode = 'single')
grViz(sarlm.sem.birds.geomean.wDR.elevation.gain.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_birds_gain_main_geomeanwDR.pdf")

####D4B_3) plot sems with cells with elevation loss only####
source ('./R/plot_sems.R')

#point estimates (main figures) for mammals
sarlm.sem.mammals.geomean.wDR.elevation.loss.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.mammals.geomean.wDR.elevation.loss.temp,standardize = 'none'),
                                                                            mode = 'single')
grViz(sarlm.sem.mammals.geomean.wDR.elevation.loss.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_mammals_loss_main_geomeanwDR.pdf")

#point estimates (main figures) for birds
sarlm.sem.birds.geomean.wDR.elevation.loss.temp.grViz <- coefsdf_to_grViz(coefs.df = coefs(sarlm.sem.birds.geomean.wDR.elevation.loss.temp,standardize = 'none'),
                                                                          mode = 'single')
grViz(sarlm.sem.birds.geomean.wDR.elevation.loss.temp.grViz) %>%
  export_svg %>% charToRaw %>% rsvg_pdf("./plots/sem_birds_loss_main_geomeanwDR.pdf")




