#this is the analysis script for the elevation_speciation repo

#load libraries
library(hexbin)
library(spdep)
library(ape)
library(piecewiseSEM)

####TO DO: delete local path here####
setwd('~/Dropbox/Work_in_progress/elevation/')

#raw_data folder contents
#

#folder structure etc
dir.create('./plots/')
dir.create('./output/')
dir.create('./output/mammals/')
dir.create('./output/birds/')
dir.create('./output/mammals/trees/')
dir.create('./output/birds/trees/')

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

####2)diversification rates####
#DR

#prepare BAMM input

####3)species occurrence data, grid####

####4)diversification measures in grid####

####5)elevation and historic changes in elevation####

####6)temperature and past temperature####



