library(ape)
library(picante)

#this script calibrates a single tree with the Meredith dates
calibrate_tree_Meredith<-function(treefile){
  if(class(treefile)=='phylo'){
    tree.MCC<-treefile
  }else{
    tree.MCC<-read.tree(treefile)
  }
  tree<-read.tree('./raw_data/mammals_Rolland2014_tree.txt')
  #this gets a table with all the nodes in the Rolland tree (useful to check what nodes they calibrated)
  node.age(tree)->phy.age
  cbind(phy.age$edge,phy.age$age, tree$edge.length)->BL.position
  max(phy.age$age)-BL.position[,3]->dist.tip
  cbind(BL.position,dist.tip)->BL.positions
  BL.positions[,5]+BL.positions[,4]->ages
  cbind(BL.positions,ages)->BL.positions
  as.data.frame(BL.positions)->node.ages
  names(node.ages)<-c("parental.node","daughter.node","dist.root","BL","dist.tip","mrca.age")
  #I use the same nodes that they used in Rolland et al PlosBiol 2014 to calibrate the BinindaEmonds supertree with Meredith et al 2011 dates
  #define the nodes and create a dataframe (Table 1 in Meredith et al)
  nodes<-list()
  #Mammalia
  nodes[[1]]<-paste('mrca: ','Pan_paniscus,',' Ornithorhynchus_anatinus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Pan_paniscus','Ornithorhynchus_anatinus')),'mrca.age'][1],1),';',sep='')
  #Theria
  nodes[[2]]<-paste('mrca: ','Pan_paniscus,',' Phascolarctos_cinereus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Pan_paniscus','Phascolarctos_cinereus')),'mrca.age'][1],1),';',sep='')
  #Placentalia
  nodes[[3]]<-paste('mrca: ','Pan_paniscus,',' Elephas_maximus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Pan_paniscus','Elephas_maximus')),'mrca.age'][1],1),';',sep='')
  #Boreoeutheria
  nodes[[4]]<-paste('mrca: ','Pan_paniscus,',' Canis_lupus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Pan_paniscus','Canis_lupus')),'mrca.age'][1],1),';',sep='')
  #Laurasiatheria
  nodes[[5]]<-paste('mrca: ','Galemys_pyrenaicus,',' Canis_lupus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Galemys_pyrenaicus','Canis_lupus')),'mrca.age'][1],1),';',sep='')
  #Euarchontoglires
  nodes[[6]]<-paste('mrca: ','Pan_paniscus,',' Mus_musculus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Pan_paniscus','Mus_musculus')),'mrca.age'][1],1),';',sep='')
  #Marsupialia
  nodes[[7]]<-paste('mrca: ','Monodelphis_domestica,',' Phascolarctos_cinereus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Monodelphis_domestica','Phascolarctos_cinereus')),'mrca.age'][1],1),';',sep='')
  #Afrotheria
  nodes[[8]]<-paste('mrca: ','Orycteropus_afer,',' Elephas_maximus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Orycteropus_afer','Elephas_maximus')),'mrca.age'][1],1),';',sep='')
  #Glires
  nodes[[9]]<-paste('mrca: ','Oryctolagus_cuniculus,',' Mus_musculus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Oryctolagus_cuniculus','Mus_musculus')),'mrca.age'][1],1),';',sep='')
  #Primates (in Rolland et al they do not fix the Primate node, they fix Haplorrhini?)
  #nodes[[10]]<-paste('mrca: ','Pan_paniscus,',' Tarsius_tarsier,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Pan_paniscus','Lemur_catta')),'mrca.age'][1],1),';',sep='')
  nodes[[10]]<-paste('mrca: ','Pan_paniscus,',' Lemur_catta,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Pan_paniscus','Lemur_catta')),'mrca.age'][1],1),';',sep='')
  #Rodentia
  nodes[[11]]<-paste('mrca: ','Cavia_porcellus,',' Mus_musculus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Cavia_porcellus','Mus_musculus')),'mrca.age'][1],1),';',sep='')
  #Xenarthra
  nodes[[12]]<-paste('mrca: ','Dasypus_novemcinctus,',' Myrmecophaga_tridactyla,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Dasypus_novemcinctus','Myrmecophaga_tridactyla')),'mrca.age'][1],1),';',sep='')
  #Scandentia
  nodes[[13]]<-paste('mrca: ','Tupaia_belangeri,',' Ptilocercus_lowii,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Tupaia_belangeri','Ptilocercus_lowii')),'mrca.age'][1],1),';',sep='')
  #Lagomorpha
  nodes[[14]]<-paste('mrca: ','Ochotona_princeps,',' Oryctolagus_cuniculus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Ochotona_princeps','Oryctolagus_cuniculus')),'mrca.age'][1],1),';',sep='')
  #Monotremata
  nodes[[15]]<-paste('mrca: ','Ornithorhynchus_anatinus,',' Tachyglossus_aculeatus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Ornithorhynchus_anatinus','Tachyglossus_aculeatus')),'mrca.age'][1],1),';',sep='')
  #Sirenia
  nodes[[16]]<-paste('mrca: ','Dugong_dugon,',' Trichechus_manatus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Dugong_dugon','Trichechus_manatus')),'mrca.age'][1],1),';',sep='')
  PATHd8.path<-'/Applications2/PATHd8/PATHd8'
  infile<-'./output/mammals/trees/input_path8d.txt'
  outfile<-sub(infile,pattern='input_',replacement='output_')
  write.tree(tree.MCC,file=infile)
  write(unlist(nodes),file=infile,append=TRUE)
  system(paste(PATHd8.path,infile,outfile,sep=' '))
  calibrated.tree<-system(paste("grep ","'^d8 tree'",' ',outfile,sep=''),intern = TRUE)
  calibrated.tree<-sub(calibrated.tree,pattern='d8 tree   : ',replacement='')
  write(calibrated.tree,file=outfile)
  calibrated.tree<-read.tree(outfile)
  file.remove(infile)
  file.remove(outfile)
  class(calibrated.tree)<-'phylo'
  calibrated.tree<-ladderize(calibrated.tree,right=FALSE)
  return(calibrated.tree)
  
}


#this script calibrates the 100 trees from the Kuhn posterior
calibrate_pseudoposterior_Meredith<-function(treesfile){
  dir.create('./output/mammals/trees/posterior_calibrated/')
  if(class(treesfile)=='multiPhylo'){
    all100trees<-treesfile
  }else{
    all100trees<-read.nexus(treesfile)
  }
  tree<-read.tree('./raw_data/mammals_Rolland2014_tree.txt')
  #this gets a table with all the nodes in the Rolland tree (useful to check what nodes they calibrated)
  node.age(tree)->phy.age
  cbind(phy.age$edge,phy.age$age, tree$edge.length)->BL.position
  max(phy.age$age)-BL.position[,3]->dist.tip
  cbind(BL.position,dist.tip)->BL.positions
  BL.positions[,5]+BL.positions[,4]->ages
  cbind(BL.positions,ages)->BL.positions
  as.data.frame(BL.positions)->node.ages
  names(node.ages)<-c("parental.node","daughter.node","dist.root","BL","dist.tip","mrca.age")
  #I use the same nodes that they used in Rolland et al PlosBiol 2014 to calibrate the BinindaEmonds supertree with Meredith et al 2011 dates
  #define the nodes and create a dataframe (Table 1 in Meredith et al)
  nodes<-list()
  #Mammalia
  nodes[[1]]<-paste('mrca: ','Pan_paniscus,',' Ornithorhynchus_anatinus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Pan_paniscus','Ornithorhynchus_anatinus')),'mrca.age'][1],1),';',sep='')
  #Theria
  nodes[[2]]<-paste('mrca: ','Pan_paniscus,',' Phascolarctos_cinereus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Pan_paniscus','Phascolarctos_cinereus')),'mrca.age'][1],1),';',sep='')
  #Placentalia
  nodes[[3]]<-paste('mrca: ','Pan_paniscus,',' Elephas_maximus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Pan_paniscus','Elephas_maximus')),'mrca.age'][1],1),';',sep='')
  #Boreoeutheria
  nodes[[4]]<-paste('mrca: ','Pan_paniscus,',' Canis_lupus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Pan_paniscus','Canis_lupus')),'mrca.age'][1],1),';',sep='')
  #Laurasiatheria
  nodes[[5]]<-paste('mrca: ','Galemys_pyrenaicus,',' Canis_lupus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Galemys_pyrenaicus','Canis_lupus')),'mrca.age'][1],1),';',sep='')
  #Euarchontoglires
  nodes[[6]]<-paste('mrca: ','Pan_paniscus,',' Mus_musculus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Pan_paniscus','Mus_musculus')),'mrca.age'][1],1),';',sep='')
  #Marsupialia
  nodes[[7]]<-paste('mrca: ','Monodelphis_domestica,',' Phascolarctos_cinereus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Monodelphis_domestica','Phascolarctos_cinereus')),'mrca.age'][1],1),';',sep='')
  #Afrotheria
  nodes[[8]]<-paste('mrca: ','Orycteropus_afer,',' Elephas_maximus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Orycteropus_afer','Elephas_maximus')),'mrca.age'][1],1),';',sep='')
  #Glires
  nodes[[9]]<-paste('mrca: ','Oryctolagus_cuniculus,',' Mus_musculus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Oryctolagus_cuniculus','Mus_musculus')),'mrca.age'][1],1),';',sep='')
  #Primates (in Rolland et al they do not fix the Primate node, they fix Haplorrhini?)
  #nodes[[10]]<-paste('mrca: ','Pan_paniscus,',' Tarsius_tarsier,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Pan_paniscus','Lemur_catta')),'mrca.age'][1],1),';',sep='')
  nodes[[10]]<-paste('mrca: ','Pan_paniscus,',' Lemur_catta,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Pan_paniscus','Lemur_catta')),'mrca.age'][1],1),';',sep='')
  #Rodentia
  nodes[[11]]<-paste('mrca: ','Cavia_porcellus,',' Mus_musculus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Cavia_porcellus','Mus_musculus')),'mrca.age'][1],1),';',sep='')
  #Xenarthra
  nodes[[12]]<-paste('mrca: ','Dasypus_novemcinctus,',' Myrmecophaga_tridactyla,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Dasypus_novemcinctus','Myrmecophaga_tridactyla')),'mrca.age'][1],1),';',sep='')
  #Scandentia
  nodes[[13]]<-paste('mrca: ','Tupaia_belangeri,',' Ptilocercus_lowii,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Tupaia_belangeri','Ptilocercus_lowii')),'mrca.age'][1],1),';',sep='')
  #Lagomorpha
  nodes[[14]]<-paste('mrca: ','Ochotona_princeps,',' Oryctolagus_cuniculus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Ochotona_princeps','Oryctolagus_cuniculus')),'mrca.age'][1],1),';',sep='')
  #Monotremata
  nodes[[15]]<-paste('mrca: ','Ornithorhynchus_anatinus,',' Tachyglossus_aculeatus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Ornithorhynchus_anatinus','Tachyglossus_aculeatus')),'mrca.age'][1],1),';',sep='')
  #Sirenia
  nodes[[16]]<-paste('mrca: ','Dugong_dugon,',' Trichechus_manatus,',' fixage=',round(node.ages[node.ages$parental.node==getMRCA(tree,tip =c('Dugong_dugon','Trichechus_manatus')),'mrca.age'][1],1),';',sep='')
  #read in nexus of all trees
  #split multiphylo to one newick per tree (this will be the input for PATHd8)
  PATHd8.path<-'/Applications2/PATHd8/PATHd8'
  calibrated.trees<-list()
  for (i in 1:length(all100trees)){
    infile<-paste('./output/mammals/trees/posterior_calibrated/',i,'_input_path8d.txt',sep='')
    outfile<-sub(infile,pattern='_input_',replacement='_output_')
    write.tree(all100trees[[i]],file=infile)
    write(unlist(nodes),file=infile,append=TRUE)
    system(paste(PATHd8.path,infile,outfile,sep=' '))
    calibrated.treeline<-system(paste("grep ","'^d8 tree'",' ',outfile,sep=''),intern = TRUE)
    calibrated.treeline<-sub(calibrated.treeline,pattern='d8 tree   : ',replacement='')
    write(calibrated.treeline,file=outfile)
    calibrated.trees[[i]]<-read.tree(outfile)
    file.remove(infile)
    file.remove(outfile)
  }
  class(calibrated.trees)<-'multiPhylo'
  write.nexus(calibrated.trees,file='./output/mammals/trees/posterior_calibrated/FritzTree.rs200k.100trees_Meredithdates.trees')
}


