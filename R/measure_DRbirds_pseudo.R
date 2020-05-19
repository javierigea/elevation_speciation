#run model in hydrogen cluster
.libPaths('/home/ji247/R/x86_64-pc-linux-gnu-library/3.5')
library(ape)
library(picante)

source('./R/measure_DR.R')

args<-commandArgs(trailingOnly = TRUE)
print(args)
replicate<-as.numeric(args[1])

DRbirds_stats_grid_pseudoposterior_replicate(replicate = replicate)


