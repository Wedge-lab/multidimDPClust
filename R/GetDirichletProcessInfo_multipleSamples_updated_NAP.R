source("interconvertMutationBurdens.R")
nucleotides = c("A","C","G","T")
#########input data: NASER#########
samplenames=read.table("samplenames")
samplenames=as.vector(samplenames$V1)
ce=read.table("cellularity")
cellularity=as.vector(ce$cellularity)
sub=read.table("subclone.files")
subclone.files=as.vector(sub$V1)
m=read.table("mutCount",header=T,stringsAsFactors=F)
library(tidyr)
info=m %>% separate(ID,into=c("chr","pos"),sep="\\_")
INFO=info[,1:2]
mutCount=array(m[,-1])
w=read.table("WTCount",header=T,stringsAsFactors=F)
WTCount=array(w[,-1])
######samplenames,cellularity,subclone.files,INFO, mutCount and WTCount are read ABOVE#####
#
#GetDirichletProcessInfo_multipleSamples<-function(samplenames, cellularity, mutCount, WTCount, subclone.files, is.male = F, out.files = NULL, phase.dir = NULL, info = NULL, SNP.phase.file = NULL, mut.phase.file = NULL){
#keep.muts.not.explained.by.CN added 020314
GetDirichletProcessInfo_multipleSamples <- function(samplenames, cellularity, mutCount, WTCount, subclone.files, is.male = T, out.files = NULL, phase.dir = NULL, info = INFO, SNP.phase.file = NULL, mut.phase.file = NULL, keep.muts.not.explained.by.CN=T){
	non.zero.indices=list()
	p.vals2 = list()
	
	if(is.null(out.files)){
		out.files = cbind(paste(samplenames,"_allDirichletProcessInfo_fromMultipleSamples.txt",sep=""),paste(samplenames,"_allDirichletProcessInfo_fromMultipleSamplesWithoutChromosomeLosses.txt",sep=""))
	}
	
	data = list()
	no.subsamples = length(samplenames)
	for(s in 1:no.subsamples){
		data[[s]] = data.frame(cbind(info,WTCount[,s],mutCount[,s]))
		names(data[[s]])[(ncol(data[[s]])-1):ncol(data[[s]])] =c("WT.count","mut.count")
		full.samplename = samplenames[s]
                data[[s]]$pos=as.numeric(data[[s]]$pos) #naser: read as character and not used in mathematical equations below -> thus NA output		
		subclone.data = read.table(subclone.files[s],header=T,stringsAsFactors=F)	
		
		data[[s]]$subclonal.CN = NA
		data[[s]]$nMaj1 = NA
		data[[s]]$nMin1 = NA
		data[[s]]$frac1 = NA
		data[[s]]$nMaj2 = NA
		data[[s]]$nMin2 = NA
		data[[s]]$frac2 = NA

		print("subclones file OK")
		
		for(r in 1:nrow(subclone.data)){
			CN = (subclone.data$nMaj1_A[r] + subclone.data$nMin1_A[r]) * subclone.data$frac1_A[r]
			if(subclone.data$frac1_A[r] != 1){
				CN = CN + (subclone.data$nMaj2_A[r] + subclone.data$nMin2_A[r]) * subclone.data$frac2_A[r]
			}
			data[[s]]$subclonal.CN[data[[s]][,1]==subclone.data$chr[r] & data[[s]][,2]>=subclone.data$startpos[r] & data[[s]][,2]<=subclone.data$endpos[r]] = CN
			data[[s]]$nMaj1[data[[s]][,1]==subclone.data$chr[r] & data[[s]][,2]>=subclone.data$startpos[r] & data[[s]][,2]<=subclone.data$endpos[r]] = subclone.data$nMaj1_A[r]
			data[[s]]$nMin1[data[[s]][,1]==subclone.data$chr[r] & data[[s]][,2]>=subclone.data$startpos[r] & data[[s]][,2]<=subclone.data$endpos[r]] = subclone.data$nMin1_A[r]
			data[[s]]$frac1[data[[s]][,1]==subclone.data$chr[r] & data[[s]][,2]>=subclone.data$startpos[r] & data[[s]][,2]<=subclone.data$endpos[r]] = subclone.data$frac1_A[r]
			data[[s]]$nMaj2[data[[s]][,1]==subclone.data$chr[r] & data[[s]][,2]>=subclone.data$startpos[r] & data[[s]][,2]<=subclone.data$endpos[r]] = subclone.data$nMaj2_A[r]
			data[[s]]$nMin2[data[[s]][,1]==subclone.data$chr[r] & data[[s]][,2]>=subclone.data$startpos[r] & data[[s]][,2]<=subclone.data$endpos[r]] = subclone.data$nMin2_A[r]
			data[[s]]$frac2[data[[s]][,1]==subclone.data$chr[r] & data[[s]][,2]>=subclone.data$startpos[r] & data[[s]][,2]<=subclone.data$endpos[r]] = subclone.data$frac2_A[r]	
		}
		
		data[[s]]$phase="unphased"
		#for(chr in unique(data[[s]][,1])){
		#	chr.inds = which(data[[s]][,1]==chr)
		#	ph.file = paste("phasing_data/",samplename,"_muts_linkedMuts_segmented_chr",chr,".txt",sep="")
		#	if(file.exists(ph.file)){
		#		phase.data[[s]] = read.table(ph.file,sep="\t",header=F,row.names=NULL,stringsAsFactors=F)
		#		indices = match(data[[s]][chr.inds,2],phase.data[[s]][,2])
		#		data[[s]]$phase[chr.inds[!is.na(indices)]] = phase.data[[s]][indices[!is.na(indices)],14]
		#	}
		#}
		#data[[s]]$phase[is.na(data[[s]]$phase)]="unphased"
		
		data[[s]]$normalCopyNumber = 2
		#assume that col1 contains chromosome
		data[[s]]$normalCopyNumber[data[[s]][,1]=="X"] = 1
		data[[s]]$mutation.copy.number = mutationBurdenToMutationCopyNumber(data[[s]]$mut.count/ (data[[s]]$mut.count + data[[s]]$WT.count) , data[[s]]$subclonal.CN, cellularity[s],data[[s]]$normalCopyNumber)
				
		data[[s]]$mutation.copy.number[is.nan(data[[s]]$subclonal.CN)]=NA
		data[[s]]$subclonal.CN[is.nan(data[[s]]$subclonal.CN)] = NA

		#convert MCN to subclonal fraction - tricky for amplified mutations
		data[[s]]$subclonal.fraction = data[[s]]$mutation.copy.number
		expected.burden.for.MCN = mutationCopyNumberToMutationBurden(rep(1,nrow(data[[s]])),data[[s]]$subclonal.CN,cellularity[s],data[[s]]$normalCopyNumber)
		non.zero.indices[[s]] = which(data[[s]]$mut.count>0 & !is.na(expected.burden.for.MCN))
		#test for mutations in more than 1 copy
		if(length(non.zero.indices[[s]])>0){
			p.vals = sapply(1:length(non.zero.indices[[s]]),function(v,e,i){prop.test(v$mut.count[i],v$mut.count[i] + v$WT.count[i],e[i],alternative="greater")$p.value},v=data[[s]][non.zero.indices[[s]],], e= expected.burden.for.MCN[non.zero.indices[[s]]])
			amplified.muts = non.zero.indices[[s]][p.vals<=0.05]
		}else{
			amplified.muts = integer(0)
		}

		data[[s]]$no.chrs.bearing.mut = 1	

		#copy numbers of subclones can only differ by 1 or 0 (as assumed when calling subclones)
		if(length(amplified.muts)>0){		
			for(a in 1:length(amplified.muts)){
				max.CN2=0
				#use phasing info - if on 'deleted' (lower CN chromosome), use the minor copy number
				if(data[[s]]$phase[amplified.muts[a]]=="MUT_ON_DELETED"){
					print("mut on minor chromosome")
					max.CN1 = data[[s]]$nMin1[amplified.muts[a]]
					frac1 = data[[s]]$frac1[amplified.muts[a]]
					frac2=0
					if(!is.na(data[[s]]$nMin2[amplified.muts[a]])){
						#swap subclones, so that the one with the higher CN is first
						if(data[[s]]$nMin2[amplified.muts[a]]>max.CN1){
							max.CN2 = max.CN1
							max.CN1 = data[[s]]$nMin2[amplified.muts[a]]
							frac2 = frac1
							frac1 = data[[s]]$frac2[amplified.muts[a]]
						}else{
							max.CN2 = data[[s]]$nMin2[amplified.muts[a]]
							frac2 = data[[s]]$frac2[amplified.muts[a]]
						}
					}					
				}else{
					max.CN1 = data[[s]]$nMaj1[amplified.muts[a]]
					frac1 = data[[s]]$frac1[amplified.muts[a]]
					frac2=0
					if(!is.na(data[[s]]$nMaj2[amplified.muts[a]])){
						#swap subclones, so that the one with the higher CN is first
						if(data[[s]]$nMaj2[amplified.muts[a]]>max.CN1){
							max.CN2 = max.CN1
							max.CN1 = data[[s]]$nMaj2[amplified.muts[a]]
							frac2 = frac1
							frac1 = data[[s]]$frac2[amplified.muts[a]]						
						}else{
							max.CN2 = data[[s]]$nMaj2[amplified.muts[a]]
							frac2 = data[[s]]$frac2[amplified.muts[a]]
						}
					}	
				}
				best.err = data[[s]]$mutation.copy.number[amplified.muts[a]] - 1
				best.CN=1
				for(j in 1:max.CN1){
					for(k in (j-1):min(j,max.CN2)){
						potential.CN = j * frac1 + k * frac2
						err = abs(data[[s]]$mutation.copy.number[amplified.muts[a]]/potential.CN-1)
						if(err<best.err){
							data[[s]]$no.chrs.bearing.mut[amplified.muts[a]] = potential.CN
							best.err=err
							best.CN = potential.CN
						}
					}
				}
				data[[s]]$subclonal.fraction[amplified.muts[a]] = data[[s]]$mutation.copy.number[amplified.muts[a]] / best.CN
			}
		}

	##########################################################################

		#test for subclonal mutations
		#test whether mut burden is less than expected value for MCN = 1
		if(length(non.zero.indices[[s]])>0){
			p.vals1 = sapply(1:length(non.zero.indices[[s]]),function(v,e,i){prop.test(v$mut.count[i],v$mut.count[i] + v$WT.count[i],e[i], alternative="less")$p.value},v=data[[s]][non.zero.indices[[s]],], e= expected.burden.for.MCN[non.zero.indices[[s]]])
			#test whether mut burden is above error rate (assumed to be 1 in 200)
			p.vals2[[s]] = sapply(1:length(non.zero.indices[[s]]),function(v,i){prop.test(v$mut.count[i],v$mut.count[i] + v$WT.count[i],0.005,alternative="greater")$p.value},v=data[[s]][non.zero.indices[[s]],])

			#test for non-zero muts causes problems for muts that have been lost in most cells
			#subclonal.muts = non.zero.indices[[s]][p.vals1<=0.05 & p.vals2[[s]]<=0.05]
			subclonal.muts = non.zero.indices[[s]][p.vals1<=0.05]
		}else{
			subclonal.muts = integer(0)
		}
		
    # use subclonal CN that minimises the difference in subclonal fraction from 1
    # NAP: get best.CN for subclonal mutation to test whether assigning to a CNA subclone better matches the observed VAF with the expected VAF (obtained from MCN * best.CN) #
    if(length(subclonal.muts)>0){
      for(a in 1:length(subclonal.muts)){
        #if there are no subclonal CNVs, don't adjust subclonal fraction
        if(is.na(data[[s]]$frac2[subclonal.muts[a]])){next}
        #assume subclonal muts are on one chromosome copy, therefore mutation copy number must be subclonal fraction of the higher CN subclone (i.e. lost in lower CN subclone) or 1 (i.e. present in both subclones)
        if(data[[s]]$nMaj1[subclonal.muts[a]]+data[[s]]$nMin1[subclonal.muts[a]] > data[[s]]$nMaj2[subclonal.muts[a]]+data[[s]]$nMin2[subclonal.muts[a]]){
          # NAP: additional loop for special case when subclonal fractions have 2:1 and 1:1 copy number states - most likely scenario is gain of 1 copy in smaller subclone vs subclonal copy loss
          if (data[[s]]$nMaj1[subclonal.muts[a]]==2 & data[[s]]$nMin1[subclonal.muts[a]]==1 & data[[s]]$nMaj2[subclonal.muts[a]]==1 & data[[s]]$nMin2[subclonal.muts[a]]==1){
            possible.subclonal.fractions = c(1+data[[s]]$frac1[subclonal.muts[a]],1) # NAP: 1+frac1 = assuming subclonal gain
          } else {
            possible.subclonal.fractions = c(data[[s]]$frac1[subclonal.muts[a]],1)
            }
          } else {
          # NAP: additional loop for special case when subclonal fractions have 1:1 and 2:1 copy number states - most likely scenario is gain of 1 copy in larger subclone vs subclonal copy loss
          if (data[[s]]$nMaj1[subclonal.muts[a]]==1 & data[[s]]$nMin1[subclonal.muts[a]]==1 & data[[s]]$nMaj2[subclonal.muts[a]]==2 & data[[s]]$nMin2[subclonal.muts[a]]==1){
            possible.subclonal.fractions = c(1+data[[s]]$frac2[subclonal.muts[a]],1) # NAP: 1+frac2 = assuming subclonal gain
          } else {
          possible.subclonal.fractions = c(data[[s]]$frac2[subclonal.muts[a]],1)
          }
        }
        best.CN = possible.subclonal.fractions[which.min(abs(data[[s]]$mutation.copy.number[subclonal.muts[a]]/possible.subclonal.fractions - 1))] # NAP: which fraction is closer to the MCN of the mutation ; min(|proportion-1|)
  
        #NAP: if you apply best.CN (which is either frac1/1+frac1 or frac2/1+frac2 and NOT 1), is the deviation of the observed VAF non-significant from expectedVAF based on MCN? I.e. are we close enough to the observed VAF?
        
        #extra test 200313 - check whether subclonal CN results in clonal mutation, otherwise subclonal CN doesn't explain subclonal MCN
        if(best.CN != 1 & prop.test(data[[s]]$mut.count[subclonal.muts[a]],data[[s]]$mut.count[subclonal.muts[a]]+data[[s]]$WT.count[subclonal.muts[a]],expected.burden.for.MCN[subclonal.muts[a]] * best.CN)$p.value > 0.05){
          data[[s]]$subclonal.fraction[subclonal.muts[a]] = data[[s]]$mutation.copy.number[subclonal.muts[a]] / best.CN
          data[[s]]$no.chrs.bearing.mut[subclonal.muts[a]] = best.CN
        }
      }
    }

		write.table(data[[s]], out.files[s,1],sep="\t",row.names=F,quote=F)
	}
		
	combined.nMaj = pmax(data[[1]]$nMaj1,data[[1]]$nMaj2,na.rm=T)
	combined.nMin = pmax(data[[1]]$nMin1,data[[1]]$nMin2,na.rm=T)
	combined.frac = data[[1]]$frac1
	for(s in 2:no.subsamples){
		combined.nMaj = cbind(combined.nMaj,pmax(data[[s]]$nMaj1,data[[s]]$nMaj2,na.rm=T))
		combined.nMin = cbind(combined.nMin,pmax(data[[s]]$nMin1,data[[s]]$nMin2,na.rm=T))
		combined.frac = cbind(combined.frac,data[[s]]$frac1)
	}
	max.maj.CNs = sapply(1:nrow(combined.nMaj),function(v,i){max(v[i,],na.rm=T)},v=combined.nMaj)
	max.min.CNs = sapply(1:nrow(combined.nMin),function(v,i){max(v[i,],na.rm=T)},v=combined.nMin)
	
	
	#June 2014 check other samples
	possible.zero.muts = list()
	non.mutant = list()
	for(s in 1:no.subsamples){
		possible.zero.muts[[s]] = intersect((1:nrow(data[[s]]))[-non.zero.indices[[s]]],which(!is.na(data[[s]]$nMin1)))
		if(length(non.zero.indices[[s]])>0){
		  possible.zero.muts[[s]] = c(possible.zero.muts[[s]],non.zero.indices[[s]][p.vals2[[s]]>0.05])
		}
		non.mutant[[s]] = possible.zero.muts[[s]]
		if(length(possible.zero.muts[[s]])>0){
			#del.indices = which(data[possible.zero.muts,paste("nMin1.",subsamplenames[s],sep="")]==0)
			#060113 - check for any reduction in CN
			del.indices = which(combined.frac[possible.zero.muts[[s]],s] == 1 & (combined.nMin[possible.zero.muts[[s]],s]<max.min.CNs[possible.zero.muts[[s]]] | combined.nMaj[possible.zero.muts[[s]],s]<max.maj.CNs[possible.zero.muts[[s]]]))
			possible.zero.muts[[s]] = possible.zero.muts[[s]][del.indices]
		}
	}
	for(s in 1:no.subsamples){
		non.del.indices = integer(0)
		if(keep.muts.not.explained.by.CN){
			for(t in (1:no.subsamples)[-s]){
				non.del.inds = which(combined.nMin[possible.zero.muts[[s]],s]>=combined.nMin[possible.zero.muts[[s]],t] & combined.nMaj[possible.zero.muts[[s]],s]>=combined.nMaj[possible.zero.muts[[s]],t])
				#non.del.inds = non.del.inds[!(non.del.inds %in% possible.zero.muts[[t]])]
				#260514
				non.del.inds = non.del.inds[!(possible.zero.muts[[s]][non.del.inds] %in% non.mutant[[t]])]
				non.del.indices = union(non.del.indices,non.del.inds)
			}
		}
		print(paste(samplenames[s],": ",length(non.del.indices)," of ",length(possible.zero.muts[[s]])," wrongly called as lost through CN loss",sep=""))
		if(length(non.del.indices)>0){
			data[[s]][possible.zero.muts[[s]][-non.del.indices],"subclonal.fraction"] = NA
			data[[s]][possible.zero.muts[[s]][-non.del.indices],"no.chrs.bearing.mut"] = 0
		}else{
			data[[s]][possible.zero.muts[[s]],"subclonal.fraction"] = NA
			data[[s]][possible.zero.muts[[s]],"no.chrs.bearing.mut"] = 0			
		}

	}
	
	for(s in 1:no.subsamples){
		write.table(data[[s]], out.files[s,2],sep="\t",row.names=F,quote=F)
	}
}

GetDirichletProcessInfo_multipleSamples(samplenames, cellularity, mutCount, WTCount, subclone.files, is.male = T, out.files = NULL, phase.dir = NULL, info = INFO, SNP.phase.file = NULL, mut.phase.file = NULL, keep.muts.not.explained.by.CN=T)
