.libPaths("/gpfs0/users/wedge/naser/R/x86_64-pc-linux-gnu-library/3.3")
print(.libPaths())
library(DPClust)

#patch to cope with zero allele frequency in some samples
library(R.utils)
source("subclone_Dirichlet_Gibbs_sampler_nD_binomial_DCW.R")
source("plotnD.DCW.R")
reassignInPackage("plotnD",pkgName="DPClust",plotnD.DCW)
reassignInPackage("subclone.dirichlet.gibbs",pkgName="DPClust",subclone.dirichlet.gibbs.DCW)
reassignInPackage("Gibbs.subclone.density.est",pkgName="DPClust",Gibbs.subclone.density.est.DCW)
reassignInPackage("Gibbs.subclone.density.est.nd",pkgName="DPClust",Gibbs.subclone.density.est.nd.DCW)
####
####starting_point####


#samplename="tenX"

subsamplenames=read.table("samplenames")
subsamplenames=as.vector(subsamplenames$V1)

#get cellularities
ce=read.table("cellularity")
cellularity=as.vector(ce$V1)


#get DP info
DP.files = paste(subsamplenames,"_allDirichletProcessInfo_fromMultipleSamplesWithoutChromosomeLosses.txt",sep="")
print(paste(DP.files))

#naser: exclude unwanted variants using the "non.zero" filter (as opposed to the original "-exclude.indices") as applied in the naser_multidimclustparallel_updated_final.R
data=list()
for (s in 1:length(DP.files)){
data[[s]]=read.table(DP.files[s],sep="\t",header=T)
}
#check original dimension
print(paste(dim(data[[1]])))
#
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
non.zero = which(rowSums(mutCount)>0 & !is.na(rowSums(totalCopyNumber)) & rowSums(copyNumberAdjustment==0)==0)
mutCount = mutCount[non.zero,]
WTCount = WTCount[non.zero,]
totalCopyNumber = totalCopyNumber[non.zero,]
copyNumberAdjustment = copyNumberAdjustment[non.zero,]
normalCopyNumber = normalCopyNumber[non.zero,]
for(s in 1:length(subsamplenames)){
data[[s]]=data[[s]][non.zero,]
}
mutation.copy.number = array(NA,dim(totalCopyNumber))
for(i in 1:length(subsamplenames)){
mutation.copy.number[,i] = mutationBurdenToMutationCopyNumber(mutCount[,i] / (mutCount[,i]+WTCount[,i]) , totalCopyNumber[,i], cellularity[i], normalCopyNumber[,i])
}
#check dimension after filtering
print(paste(dim(data[[1]])))
#
#
no.iters=1000
if(nrow(mutation.copy.number)>50000){
	no.iters=250
}

output_folder = "DP_output"
if(!file.exists(output_folder)){
	dir.create(output_folder)
}

#This function does everything. It may be better to run separate function to perform Gibbs sampling and mutation assignment
DirichletProcessClustering(mutCount = mutCount, WTCount = WTCount, totalCopyNumber = totalCopyNumber, copyNumberAdjustment = copyNumberAdjustment,
mutation.copy.number = mutation.copy.number, cellularity = cellularity, output_folder = output_folder, no.iters = no.iters, no.iters.burn.in = ceiling(no.iters/5),
subsamplesrun = subsamplenames,samplename=samplename, conc_param = 1, cluster_conc = 5, mut.assignment.type = 1, most.similar.mut = NA, mutationTypes = NA, min.frac.snvs.cluster = NA)
