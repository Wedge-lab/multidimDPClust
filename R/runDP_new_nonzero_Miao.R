# setwd("/mnt/bmh01-rds/UoOxford_David_W/r44399mg/CRPM/Code/msDPClust_seperate_12/R")
# .libPaths("/gpfs0/users/wedge/naser/R/x86_64-pc-linux-gnu-library/3.3")
print(.libPaths())

source("../R/DirichletProcessClustering.R")
source("../R/subclone_Dirichlet_Gibbs_sampler_nD_binomial.R")
source("../R/AssignMutations.R")
source("../R/DensityEstimator.R")
source("../R/InformationCriterions.R")
source("../R/PlotDensities.R")
source("../R/SampleMutations.R")
source("../R/interconvertMutationBurdens.R")


#patch to cope with zero allele frequency in some samples
library(R.utils)

# source("subclone_Dirichlet_Gibbs_sampler_nD_binomial_DCW.R")
# source("plotnD.DCW.R")

#reassignInPackage("plotnD",pkgName="DPClust",plotnD.DCW)
#reassignInPackage("subclone.dirichlet.gibbs",pkgName="DPClust",subclone.dirichlet.gibbs.DCW)
#reassignInPackage("Gibbs.subclone.density.est",pkgName="DPClust",Gibbs.subclone.density.est.DCW)
#reassignInPackage("Gibbs.subclone.density.est.nd",pkgName="DPClust",Gibbs.subclone.density.est.nd.DCW)
####
####starting_point####

args <- commandArgs(trailingOnly = TRUE)
case = args[1]
outpath = args[2]

run_dir = paste0(outpath, "/", case, "/")

setwd(run_dir)

output_folder = "./04_GSoutput"
if(!file.exists(output_folder)){
  dir.create(output_folder)
}

# samplename="patientID"
# subsamplenames=read.table("samplenames.txt")


samplename=case 

sample_info = read.delim2("./DPinput/sample_info.txt")
samplenames = sample_info$samplenames
cellularity = as.numeric(sample_info$cellularity)
subsamplenames=  sapply(strsplit(samplenames, "_"), function(x) paste(x[-1], collapse = "_"))

#get cellularities
# purity_file = paste0("/mnt/bmh01-rds/UoOxford_David_W/r44399mg/CRPM/BB2210/",samplenames, "/" , samplenames, "_rho_and_psi.txt")
# no_samples = length(samplenames)

# for(s in 1:no_samples){
#   cellularity[s] = read.table(purity_file[s],header=T)$rho[2]
# }
print(samplenames)
print(cellularity)

#get DP info
DP.files = paste("./DPinput/",samplenames,"_allDirichletProcessInfo_fromMultipleSamplesWithoutChromosomeLosses.txt",sep="")

print(paste(DP.files))

#naser: exclude unwanted variants using the "non.zero" filter (as opposed to the original "-exclude.indices") as applied in the naser_multidimclustparallel_updated_final.R
data=list()
for (s in 1:length(DP.files)){
  print(s)
  data[[s]]=read.table(DP.files[s],sep="\t",header=T)
}
#check original dimension
print(paste(dim(data[[1]])))
#
# data[[1]] =read.table(DP.files[1],sep="\t",header=T)
# a= data[[1]]
WTCount = data[[1]]$WT.count
mutCount = data[[1]]$mut.count
totalCopyNumber = data[[1]]$subclonal.CN
copyNumberAdjustment = data[[1]]$no.chrs.bearing.mut
normalCopyNumber = data[[1]]$normalCopyNumber
for(s in 2:length(subsamplenames)){
  WTCount = cbind(WTCount,data[[s]]$WT.count)
  mutCount = cbind(mutCount,data[[s]]$mut.count)
  totalCopyNumber = cbind(totalCopyNumber,data[[s]]$subclonal.CN)
  copyNumberAdjustment = cbind(copyNumberAdjustment,data[[s]]$no.chrs.bearing.mut)
  normalCopyNumber = cbind(normalCopyNumber,data[[1]]$normalCopyNumber)
}
# changed
non.zero = which(rowSums(mutCount)>0 & !is.na(rowSums(totalCopyNumber)) & rowSums(copyNumberAdjustment!=0)>0)
# original
# non.zero = which(rowSums(mutCount)>0 & !is.na(rowSums(totalCopyNumber)) & rowSums(copyNumberAdjustment==0)==0)
mutCount = mutCount[non.zero,]
WTCount = WTCount[non.zero,]
totalCopyNumber = totalCopyNumber[non.zero,]
copyNumberAdjustment = copyNumberAdjustment[non.zero,]
normalCopyNumber = normalCopyNumber[non.zero,]
for(s in 1:length(subsamplenames)){
  data[[s]]=data[[s]][non.zero,]
  if (s == 1){
    #check dimension after filtering
    print(paste("Number of retained SNVs =", nrow(data[[1]])))
    # Chr and Pos of filtered SNVs for post-clustering mapping
    write.table(data[[s]][,1:2],paste0(output_folder,"/",samplename,"_union_filtered_SNVs.txt"),col.names=T,row.names=F,quote=F,sep="\t")
  }
}
mutation.copy.number = array(NA,dim(totalCopyNumber))
for(i in 1:length(subsamplenames)){
  mutation.copy.number[,i] = mutationBurdenToMutationCopyNumber(mutCount[,i] / (mutCount[,i]+WTCount[,i]) , totalCopyNumber[,i], cellularity[i], normalCopyNumber[,i])
}
#
#
no.iters=1250
# if(nrow(mutation.copy.number)>50000){
#   no.iters=250
# }


#This function does everything. It may be better to run separate function to perform Gibbs sampling and mutation assignment
DirichletProcessClustering(mutCount = mutCount, WTCount = WTCount, totalCopyNumber = totalCopyNumber, copyNumberAdjustment = copyNumberAdjustment,
                           mutation.copy.number = mutation.copy.number, cellularity = cellularity, output_folder = output_folder, no.iters = no.iters, no.iters.burn.in = ceiling(no.iters/5),
                           subsamplesrun = subsamplenames,samplename=samplename, conc_param = 1, cluster_conc = 5, mut.assignment.type = 1, most.similar.mut = NA, mutationTypes = NA,max.considered.clusters=20)
# min.frac.snvs.cluster = NA
