library(Hmisc)
library(stringr)
library(ape)
library(geiger)
################################################################################################################
################function to check species vector against IUCN mammals taxonomy
#dictionary file is the path to iucnmammaltaxonomy.csv (downloaded from IUCN)
#species is a vector with species to synonymise
iucn.taxonomy.synonyms<-function(species,dictionaryfile){
  #read file and delete extra columns
  table.synonyms<-read.csv(dictionaryfile,header=T)
  table.synonyms<-data.frame(table.synonyms$Genus,table.synonyms$Species,paste(table.synonyms$Genus,table.synonyms$Species,sep=' '),table.synonyms$Synonyms,stringsAsFactors = FALSE)
  colnames(table.synonyms)<-c('Genus','Species','Binomial','Synonyms')
  #clean extra characters
  table.synonyms$Synonyms<-sub(table.synonyms$Synonyms,pattern='  ',replacement='')
  #change _ in inputvector to ' '
  species<-sub(species,pattern='_',replacement=' ')
  species.result<-vector('character',length=length(species))
  #separate all synonyms using a pipe
  table.synonyms$Synonyms<-sub(table.synonyms$Synonyms,pattern=',',replacement='\\|')
  #split all synonyms into different columns
  all.synonyms<-str_split_fixed(table.synonyms$Synonyms,'\\|',n=100)
  #delete empty columns
  all.synonyms<-as.data.frame(all.synonyms[, colSums(all.synonyms != "") != 0])
  #calculate max number of columns (=max number of synonyms)
  col.synonyms<-ncol(all.synonyms)
  table.synonyms<-cbind(table.synonyms,all.synonyms)
  table.synonyms$Synonyms<-NULL
  #check match for each species in input vector and get corresponding Binomial from IUCN taxonomy
  for (i in 1:length(species)){
    index<-0
    #this is for species which are correct in inputvector
    if (species[i] %in% as.character(table.synonyms$Binomial)){
      index<-which(table.synonyms$Binomial==species[i])
      species.result[i]<-as.character(table.synonyms$Binomial[index])
      next
    }
    #check if there is a synonym present in table
    else for (a in 4:(col.synonyms+3)){
      
      if (species[i] %in% as.character(table.synonyms[,a])){
        index<-which(table.synonyms[,a]==species[i])
        #if there are two or more synonyms take the last one
        if (length(index)>1){
          #print (species[i])
          index<-length(index)
        }
        species.result[i]<-as.character(table.synonyms$Binomial[index])
      }
    }
    if (index==0){
      species.result[i]<-species[i]
    }
  }
  species.result<-sub(species.result,pattern=' ',replacement='_')
  return(species.result)
}

################function to check species vector against BIRDLIFE taxonomy
#dictionary file is the path to birdlife_taxonomy_v8_27062016 (downloaded from Birdlife)
#species is a vector with species to synonymise
birdlife.taxonomy.synonyms<-function(species,dictionaryfile){
  #read file and delete extra columns
  table.synonyms<-read.table(dictionaryfile,header=T,sep='\t')
  table.synonyms<-data.frame(table.synonyms$Family.name,table.synonyms$Scientific.name,table.synonyms$BirdLife.taxonomic.treatmen,table.synonyms$Synonyms,stringsAsFactors = FALSE)
  colnames(table.synonyms)<-c('Family','Binomial','Taxonomic.Status','Synonyms')
  #Taxonomic status "BirdLife taxonomic treatment:R = recognised as a species;NR = not recognised as a species;UR = under review"
  #select R only
  table.synonyms<-table.synonyms[table.synonyms$Taxonomic.Status=='R',]
  #change _ in inputvector to ' '
  species<-sub(species,pattern='_',replacement=' ')
  species.result<-vector('character',length=length(species))
  #separate all synonyms using a pipe
  table.synonyms$Synonyms<-gsub(table.synonyms$Synonyms,pattern='; ',replacement='\\|')
  #split all synonyms into different columns
  all.synonyms<-str_split_fixed(table.synonyms$Synonyms,'\\|',n=100)
  #delete empty columns
  all.synonyms<-as.data.frame(all.synonyms[, colSums(all.synonyms != "") != 0])
  #calculate max number of columns (=max number of synonyms)
  col.synonyms<-ncol(all.synonyms)
  table.synonyms<-cbind(table.synonyms,all.synonyms)
  table.synonyms$Synonyms<-NULL
  #check match for each species in input vector and get corresponding Binomial from Birdlife taxonomy
  for (i in 1:length(species)){
    index<-0
    #this is for species which are correct in inputvector
    if (species[i] %in% as.character(table.synonyms$Binomial)){
      index<-which(table.synonyms$Binomial==species[i])
      species.result[i]<-as.character(table.synonyms$Binomial[index])
      next
    }
    #check if there is a synonym present in table
    else for (a in 4:(col.synonyms+3)){
      
      if (species[i] %in% as.character(table.synonyms[,a])){
        index<-which(table.synonyms[,a]==species[i])
        #if there are two or more synonyms take the last one
        if (length(index)>1){
          print (species[i])
          index<-length(index)
        }
        species.result[i]<-as.character(table.synonyms$Binomial[index])
      }
    }
    if (index==0){
      species.result[i]<-species[i]
    }
  }
  species.result<-sub(species.result,pattern=' ',replacement='_')
  return(species.result)
}

################function to return duplicated tips from tree (will keep one instance of duplicated taxa)
#tree is a tree object
return.duplicate.tips<-function(tree){
  tree.duplicated<-tree$tip.label[duplicated(tree$tip.label)]
  tip.duplicated<-lapply(tree.duplicated, function(x) {nodes<-which(tree$tip.label==x); nodes<-nodes[-1];return(nodes)})
  tip.duplicated<-unlist(tip.duplicated)
  return(tip.duplicated)
}
################function to check species vector and pull out marine species (families listed inside)
#dictionary file is the path to iucnmammaltaxonomy.csv (downloaded from IUCN)
#species is a vector with species to check for marine species
get.marine.species<-function(species,dictionaryfile){
  table.synonyms<-read.csv(dictionaryfile,header=T)
  table.synonyms<-data.frame(table.synonyms$Genus,table.synonyms$Species,table.synonyms$Family,paste(table.synonyms$Genus,table.synonyms$Species,sep=' '),stringsAsFactors = FALSE)
  colnames(table.synonyms)<-c('Genus','Species','Family','Binomial')
  #this vector stores the families to be pulled out from species vector
  marine.families<-c('Balaenidae','Neobalaenidae','Balaenopteridae','Eschrichtiidae','Physeteridae','Kogiidae','Monodontidae','Ziphiidae','Delphinidae','Phocoenidae','Platanistidae','Iniidae','Pontoporiidae','Trichechidae','Dugongidae','Otariidae','Odobenidae','Phocidae')
  species<-sub(species,pattern='_',replacement=' ')
  species.marine<-vector('character',length=length(species))
  for (i in 1:length(species)){
    if (species[i] %in% as.character(table.synonyms$Binomial)){
      index<-which(table.synonyms$Binomial==species[i])
      if(capitalize(tolower(table.synonyms[index,'Family'])) %in% as.character(marine.families)){
        species.marine[i]<-species[i]
      }
    }
  }
  species.marine<-sub(species.marine,pattern=' ',replacement='_')
  species.marine<-species.marine[species.marine != ""]
  return(species.marine)
}
#this deals with duplicate species in a *_species_gridoccurrencetable by merging the ranges of the duplicates in a single entry
merge_duplicates_in_speciesgridoccurrence_table<-function(duplicated.species,species.gridoccurrence.table){
  for (i in 1:length(duplicated.species)){
    new.cells<-vector()
    new.range<-vector()
    new.cells<-unlist(strsplit(paste(species.gridoccurrence.table[grep(duplicated.species[i],species.gridoccurrence.table$spp),'cells'],collapse=''),split=' '))
    new.cells<-new.cells[-which(new.cells=='')]
    new.range<-length(new.cells)
    new.cells<-c(' ',new.cells)
    new.cells<-c(new.cells,' ')
    new.cells<-paste(new.cells,collapse=' ')
    species.gridoccurrence.table[grep(duplicated.species[i],species.gridoccurrence.table$spp),'cells']<-rep(new.cells,length(grep(duplicated.species[i],species.gridoccurrence.table$spp)))
    species.gridoccurrence.table[grep(duplicated.species[i],species.gridoccurrence.table$spp),'range.cells']<-rep(new.range,length(grep(duplicated.species[i],species.gridoccurrence.table$spp)))
    species.gridoccurrence.table[!duplicated(species.gridoccurrence.table),]
  }
  species.gridoccurrence.table<-species.gridoccurrence.table[!duplicated(species.gridoccurrence.table),]
  return(species.gridoccurrence.table)
}







###########wrapper for mammals to standardise tree + trait in iucn

taxonomy_tree_trait_mammals<-function(treefile,traitfile,dictionaryfile){
  tree<-read.tree(treefile)
  #run tree tips through IUCN
  tree$tip.label<-iucn.taxonomy.synonyms(tree$tip.label,dictionaryfile)
  #remove duplicate tips
  tip.duplicated<-remove.duplicate.tips(tree)
  tree<-drop.tip(tree,tip.duplicated)
  #run trait dataset through IUCN
  trait<-read.table(traitfile,header=T,sep='\t')
  colnames(trait)<-c('species','trait')
  trait$species<-iucn.taxonomy.synonyms(trait$species,dictionaryfile)
  #remove duplicate taxa
  trait<-unique(trait)
  #clean marine mammals from tree
  species.marine.tree<-get.marine.species(tree$tip.label,dictionaryfile)
  tree.terrestrial<-drop.tip(tree,species.marine.tree)
  #remove duplicate tips
  tip.terrestrial.duplicated<-remove.duplicate.tips(tree.terrestrial)
  tree.terrestrial<-drop.tip(tree.terrestrial,tip.terrestrial.duplicated)
  #clean marine mammals from trait dataset (not necessary, IUCN data was downloaded for terrestrial only)
  #species.marine.trait<-get.marine.species(trait$species,dictionaryfile)
  #check trait tree correspondence
  trait.vector<-trait$trait
  names(trait.vector)<-trait$species
  name.check<-name.check(tree.terrestrial,trait.vector)
  tree.not.data<-name.check$tree_not_data
  data.not.tree<-name.check$data_not_tree
  #drop tips that are in tree but not in data
  tree.terrestrial<-drop.tip(tree.terrestrial,tree.not.data)
  name.check<-name.check(tree.terrestrial,trait.vector)
  #write tree and table to output files
  write.tree(tree.terrestrial,'./output/mammals/mammals_Rolland_terrestrial_IUCN.tree')
  #replace 2 with 0, create a binary trait (0 = not endemic to hotspot, 1 = endemic)
  write.table(trait,'./output/mammals/all_spp_input_terrestrial_IUCN.txt',row.names = F,quote=F,sep='\t')
  #calculate sampling proportions (f)
  #f.hisse<-calculate.sampling.proportion(tree.terrestrial,trait)
  #f.hisse<-unname(f.hisse)
  
}

###########wrapper for birds to standardise tree + trait in iucn

taxonomy_tree_trait_birds<-function(treefile,traitfile,dictionaryfile){
  tree<-read.tree(treefile)
  #run tree tips through IUCN
  tree$tip.label<-birdlife.taxonomy.synonyms(tree$tip.label,dictionaryfile)
  #remove duplicate tips
  tip.duplicated<-remove.duplicate.tips(tree)
  tree<-drop.tip(tree,tip.duplicated)
  #run trait dataset through IUCN
  trait<-read.table(traitfile,header=T,sep='\t')
  colnames(trait)<-c('species','trait')
  trait$species<-birdlife.taxonomy.synonyms(trait$species,dictionaryfile)
  #remove duplicate taxa
  trait<-unique(trait)
  ##clean marine mammals from tree
  ##species.marine.tree<-get.marine.species(tree$tip.label,dictionaryfile)
  ##tree.terrestrial<-drop.tip(tree,species.marine.tree)
  ##remove duplicate tips
  ##tip.terrestrial.duplicated<-remove.duplicate.tips(tree.terrestrial)
  ##tree.terrestrial<-drop.tip(tree.terrestrial,tip.terrestrial.duplicated)
  #clean marine mammals from trait dataset (not necessary, IUCN data was downloaded for terrestrial only)
  #species.marine.trait<-get.marine.species(trait$species,dictionaryfile)
  #check trait tree correspondence
  trait.vector<-trait$trait
  names(trait.vector)<-trait$species
  name.check<-name.check(tree,trait.vector)
  tree.not.data<-name.check$tree_not_data
  data.not.tree<-name.check$data_not_tree
  #drop tips that are in tree but not in data
  tree<-drop.tip(tree,tree.not.data)
  name.check<-name.check(tree,trait.vector)
  #write tree and table to output files
  write.tree(tree,'./output/birds/BirdzillaEricson_1000Trees_MCC_birdlife.tree')
  #replace 2 with 0, create a binary trait (0 = not endemic to hotspot, 1 = endemic)
  write.table(trait,'./output/birds/all_spp_hisse_input_birdlife.txt',row.names = F,quote=F,sep='\t')
  #calculate sampling proportions (f)
  #f.hisse<-calculate.sampling.proportion(tree.terrestrial,trait)
  #f.hisse<-unname(f.hisse)
  
}


taxonomy_tree_mammals<-function(tree,dictionaryfile){
  #run tree tips through IUCN
  tree$tip.label<-iucn.taxonomy.synonyms(tree$tip.label,dictionaryfile)
  #remove duplicate tips
  tip.duplicated<-remove.duplicate.tips(tree)
  tree<-drop.tip(tree,tip.duplicated)
  #clean marine mammals from tree
  species.marine.tree<-get.marine.species(tree$tip.label,dictionaryfile)
  tree.terrestrial<-drop.tip(tree,species.marine.tree)
  #remove duplicate tips
  tip.terrestrial.duplicated<-remove.duplicate.tips(tree.terrestrial)
  tree.terrestrial<-drop.tip(tree.terrestrial,tip.terrestrial.duplicated)
  return(tree.terrestrial)
}



#######################function to calculate sampling proportions ('f') for HiSSE, BiSSE, etc
#tree is a tree object
#trait is a dataframe with species and trait column
calculate.sampling.proportion<-function(treefile,traitfile){
  tree<-read.tree(treefile)
  trait<-read.table(traitfile,header=T,sep='\t')
  colnames(trait)<-c('species','trait')
  # trait[trait$trait==2,]$trait<-0
  tree.trait<-merge(data.frame(tree$tip.label),trait,by.x='tree.tip.label',by.y='species')
  f1<-table(tree.trait$trait)[1]/table(trait$trait)[1]
  f2<-table(tree.trait$trait)[2]/table(trait$trait)[2]
  return(c(f1,f2))
}



