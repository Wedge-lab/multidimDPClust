# .libPaths(c("/gpfs0/users/wedge/naser/R/x86_64-pc-linux-gnu-library/3.3","/gpfs0/users/wedge/frangou/R/x86_64-pc-linux-gnu-library/3.3"))
#args=commandArgs(TRUE)
#run = as.integer(args[1])
#only ONE PATIENT: so run=1
source("../R/interconvertMutationBurdens.R")
run=1
print(paste(run))
##### miao change 1000 to 250 ,p8
noiters=1250
conc_param=1
cluster_conc = 5
density.smooth = 0.01
burn.in = 250


args <- commandArgs(trailingOnly = TRUE)
case = args[1]
outpath = args[2]

samplenames=c(case)

data_dir = paste0(outpath, "/", case, "/DPinput")

setwd(data_dir)

sample_info = read.delim2("./sample_info.txt")
samplenames_list = sample_info$samplenames
subsamplenames_list= sapply(strsplit(samplenames_list, "_"), function(x) paste(x[-1], collapse = "_"))



subsamples = list()
subsamples[[run]] = subsamplenames_list

samplename = samplenames[run]

cellularity = as.numeric(sample_info$cellularity)

print(cellularity)

no.subsamples = length(subsamples[[run]])
choose.number = 3
choose.from = no.subsamples
no.perms = choose(no.subsamples,3)



get.subsample.indices<-function(choose.index,choose.number,choose.from){
	subsample.indices = 1:choose.number
	temp.choose.index = choose.index
	for(i in 1:choose.number){
		last.subtotal = 0
		subtotal = choose(choose.from-subsample.indices[i],choose.number-i)
		while(temp.choose.index>subtotal){
			subsample.indices[i] = subsample.indices[i] + 1
			last.subtotal = subtotal
			subtotal = subtotal + choose(choose.from-subsample.indices[i],choose.number-i)			
		}
		if(i<choose.number){
			subsample.indices[i+1] = subsample.indices[i]+1
		}
		temp.choose.index = temp.choose.index - last.subtotal
	}
	return(subsample.indices)
}
print(paste(no.subsamples[[run]]))
#required for mutation.copy.number and copyNumberAdjustment


data=list()
for(s in 1:no.subsamples){
	data[[s]]=read.table(paste(samplenames_list[s],"_allDirichletProcessInfo_fromMultipleSamplesWithoutChromosomeLosses.txt",sep=""),sep="\t",header=T)
	#data[[s]]=read.table(paste(samplename,subsamples[[run]][subsample.indices[s]],"_muts_withAllSubclonalCNinfoAndWithoutMutsOnDeletedChrs_Dec2013.txt",sep=""),sep="\t",header=T)
}
#### miao change, this line is for input line 63
run_dir = paste0(outpath,"/",case, "/")
# set file path for reading data for input,Line 71
setwd(run_dir)

last_output_folder = "05_PAoutput" # read results from assign mutation step, accordingly changed line 120,147,168
new_output_folder = "06_Finaloutput"  # for files generated in this step
if(!file.exists(new_output_folder)){
  dir.create(new_output_folder)
}

#
#
# library(MASS)
# library(MCMCpack)

# library(mvtnorm)
#source("/lustre/scratch109/sanger/dw9/2D_Dirichlet_Process/multidimensionalDensityEstimator.R")
#source("/lustre/scratch109/sanger/dw9/2D_Dirichlet_Process/interconvertMutationBurdens.R")
# library(DPClust)
#
#
#required for plotting the PDF: mutation.copy.number/copyNumberAdjustment
WTCount = data[[1]]$WT.count
mutCount = data[[1]]$mut.count
totalCopyNumber = data[[1]]$subclonal.CN
copyNumberAdjustment = data[[1]]$no.chrs.bearing.mut
normalCopyNumber = data[[1]]$normalCopyNumber
for(s in 2:length(subsamples[[run]])){
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
for(s in 1:length(subsamples[[run]])){
	data[[s]]=data[[s]][non.zero,]
}
mutation.copy.number = array(NA,dim(totalCopyNumber))
for(i in 1:length(subsamples[[run]])){
	mutation.copy.number[,i] = mutationBurdenToMutationCopyNumber(mutCount[,i] / (mutCount[,i]+WTCount[,i]) , totalCopyNumber[,i], cellularity[i], normalCopyNumber[,i])
}

#naser: code above updated - 190918

no.muts = nrow(data[[1]]) #naser: nrow(info) until runDP filter has changed from -exclude.indices to non.zero ??
print(paste(no.muts))
node.assignments=NULL


############################# Error in no.perms : object 'no.perms' not found
for(choose.index in 1:no.perms){
	print(choose.index)
	subsample.indices = get.subsample.indices(choose.index,choose.number,choose.from)
	print(subsample.indices)
	perm.assignments = read.table(paste(last_output_folder,"/",samplename,"_",paste(subsample.indices,collapse="_"),"_DP_and_cluster_info_0.01.txt",sep=""),header=T,sep="\t")

	if(is.null(node.assignments)){
		node.assignments = array(data.matrix(perm.assignments[ncol(perm.assignments)]),c(no.muts,1))
	}else{
		node.assignments = cbind(node.assignments,data.matrix(perm.assignments[ncol(perm.assignments)]))
	}
}

#get different combinations of assignments from the different permutations
unique.assignments = unique(node.assignments)
write.table(unique.assignments,paste(new_output_folder,"/",samplename,"_matchedClustersInParallelRuns.txt",sep=""),sep="\t",row.names=F,quote=F)

no.consensus.nodes = nrow(unique.assignments)
print(paste("#consensus nodes = ",no.consensus.nodes,sep=""))
consensus.assignments = vector(mode="numeric",length = no.consensus.nodes)	
for(n in 1:no.consensus.nodes){
	print(n)
	consensus.assignments[sapply(1:no.muts,function(a,u,i){all(u==a[i,])},a=node.assignments,u=unique.assignments[n,])]=n
}

#get probabilities of assignment to each set of assignments
prob.consensus.assignments = array(1,c(no.muts,no.consensus.nodes)) 
#prob.consensus.assignments = array(0,c(no.muts,no.consensus.nodes)) #get average probability, rather than product
for(choose.index in 1:no.perms){
	print(choose.index)
	subsample.indices = get.subsample.indices(choose.index,choose.number,choose.from)
	perm.assignments = read.table(paste(last_output_folder,"/",samplename,"_",paste(subsample.indices,collapse="_"),"_DP_and_cluster_info_0.01.txt",sep=""),header=T,sep="\t")
	for(n in 1:no.consensus.nodes){
		prob.consensus.assignments[,n] = prob.consensus.assignments[,n] * perm.assignments[,unique.assignments[n,choose.index]+2]
		#prob.consensus.assignments[,n] = prob.consensus.assignments[,n] + perm.assignments[,unique.assignments[n,choose.index]+2]
	}
}
#prob.consensus.assignments = prob.consensus.assignments / no.perms

out = data[[1]][,1:2]
for(s in 1:no.subsamples){
out = cbind(out,data[[s]]$subclonal.fraction)
}
out = cbind(out,prob.consensus.assignments)
names(out) = c("chr","pos",paste(samplename,subsamples[[run]],"_subclonalFraction",sep=""),paste("prob.cluster",1:no.consensus.nodes,sep=""))
write.table(out,paste(new_output_folder,"/",samplename,"_allClusterProbabilitiesFromParallelRuns.txt",sep=""),sep="\t",row.names=F,quote=F)

##############################################################################################################################################
all.CIs = array(NA,c(no.consensus.nodes,no.subsamples,no.perms,2))
for(choose.index in 1:no.perms){
	subsample.indices = get.subsample.indices(choose.index,choose.number,choose.from)
	print(subsample.indices)
	CIs = data.matrix(read.table(paste(last_output_folder,"/",samplename,"_",paste(subsample.indices,collapse="_"),"_confInts_",density.smooth,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F))
	for(n in 1:no.consensus.nodes){
		all.CIs[n,subsample.indices,choose.index,1] = CIs[unique.assignments[n,choose.index],2*(1:choose.number)-1]
		all.CIs[n,subsample.indices,choose.index,2] = CIs[unique.assignments[n,choose.index],2*(1:choose.number)]
	}
}
median.CIs = array(NA,c(no.consensus.nodes,no.subsamples,2))
for(n in 1:no.consensus.nodes){
	for(s in 1:no.subsamples){
		for(c in 1:2){
			median.CIs[n,s,c] = median(all.CIs[n,s,,c],na.rm=T)
		}
	}
}
median.CIs.2D = t(sapply(1:no.consensus.nodes,function(m,i){as.vector(t(m[i,,]))},m=median.CIs))
write.table(cbind(1:no.consensus.nodes,median.CIs.2D,table(consensus.assignments)),sep="\t",row.names=F,quote=F,col.names=c("cluster.no",paste(rep(paste(samplename,subsamples[[run]],sep=""),each=2),rep(c("lowerCI","upperCI"),times=no.subsamples),sep="_"),"no.muts"),paste(new_output_folder,"/",samplename,"_consensusClustersByParallelNodeAssignment.txt",sep=""))
out = data[[1]][,1:2]
for(s in 1:no.subsamples){
	out = cbind(out,data[[s]]$subclonal.fraction)
}
out = cbind(out,consensus.assignments)
names(out) = c("chr","pos",paste(samplename,subsamples[[run]],"_subclonalFraction",sep=""),"cluster.no")
write.table(out,paste(new_output_folder,"/",samplename,"_allClusterassignmentsFromParallelRuns.txt",sep=""),sep="\t",row.names=F,quote=F)

pdf(paste(new_output_folder,"/","consensus_cluster_assignment_",samplename,"_combined_",density.smooth,".pdf",sep=""),height=4,width=4)
#its hard to distinguish more than 8 different colours
max.cols = 8
cols = rainbow(min(max.cols,no.consensus.nodes))
plot.data = mutation.copy.number/copyNumberAdjustment
plot.data[is.na(plot.data)]=0
for(i in 1:(no.subsamples-1)){
	for(j in (i+1):no.subsamples){
		plot(plot.data[,i],plot.data[,j],type = "n",xlab = paste(samplename,(subsamples[[run]][subsample.indices])[i]," subclonal fraction",sep=""), ylab = paste(samplename,subsamples[[run]][j]," subclonal fraction",sep=""),xlim = c(0,max(plot.data[,i])*1.25))
		for(n in 1:no.consensus.nodes){
			pch=20 + floor((n-1)/max.cols)
			#pch is not implmeneted above 25
			if(pch>25){
				pch=pch-20
			}
			points(plot.data[,i][consensus.assignments==n],plot.data[,j][consensus.assignments==n],col=cols[(n-1) %% max.cols + 1],pch=pch,cex=0.5)
		}
		pch=20 + floor((0:(no.consensus.nodes-1))/max.cols)
		pch[pch>25] = pch[pch>25]-20
		legend(max(plot.data[,i])*1.05,max(plot.data[,j]),legend = 1:no.consensus.nodes,col=cols[(0:(no.consensus.nodes-1)) %% max.cols + 1],pch=pch,cex=1)
	}
}	
dev.off()

q(save="no")
