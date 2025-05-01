# library(dpclust3p)
# library(DPClust)

create_loci_file=function(vcf_file,loci_file_name, pass_only=TRUE, exclude_chrY = TRUE){
  read_data = tryCatch(read.delim(vcf_file, comment.char="#", header=F, stringsAsFactor=F, nrows=1),
                       error = function(e) NULL)
  if (!is.null(read_data)) {
    vcf.cols = ncol(read_data)
    vcf.cols.default = 10 # vcf file standard contains 10 columns
    vcf.colClasses = c(NA, NA, "NULL", NA, NA, "NULL", NA, rep("NULL", 3+(vcf.cols-vcf.cols.default)))
    vcf.loci = read.delim(vcf_file, comment.char="#", header=F, stringsAsFactor=F, colClasses=vcf.colClasses)
    if(pass_only)
      vcf.loci = vcf.loci[vcf.loci[,5]=="PASS",-5]
    else
      vcf.loci = vcf.loci[,-5]
    colnames(vcf.loci) = c("chromosome", "pos", "ref","alt")
    if(exclude_chrY){
      vcf.loci = vcf.loci[vcf.loci[,1]!="Y" & vcf.loci[,1]!="chrY",]
    }
    write.table(vcf.loci,loci_file_name,quote=F,row.names=F,col.names=F,sep="\t")
  }
}


args <- commandArgs(trailingOnly = TRUE)
case = args[1]
snv_path = args[2]
run_info_path = args[3]
outpath = args[4]


snvInput=snv_path
outputDir = paste0(outpath, "/", case, "/DPinput/")

print(outputDir)
if (!dir.exists(outputDir)) {
  dir.create(outputDir, recursive = TRUE)
} 

sample_info_all = read.delim2(paste0(run_info_path, "/run_info.txt"))

sample_info <- sample_info_all[sample_info_all$case == case, ]
samplenames = sample_info$Tumour

df <- data.frame(samplenames)

print(df)

write.table(df, file = paste0(outputDir,"/sample_name.txt"), sep = "\t", quote = FALSE, row.names = FALSE)


for (samplename in samplenames){
  print(samplename)
  vcf_file=paste0(snvInput, samplename, "/0001.vcf")
  # AF_file = paste0(outputDir,"/",samplename,"_allelefrequencies.txt")
  loci_file = paste0(outputDir,"/",samplename,"_loci.txt")
  # dumpCounts.Sanger(vcf_file,AF_file)
  # AFs = read.table(paste0(outputDir,"/",samplename,"_allelefrequencies.txt"),comment.char="",stringsAsFactors = F,header=T,sep="")
  create_loci_file(vcf_file,loci_file)
  # loci = read.table(paste0(outputDir,"/",samplename,"_loci.txt"),comment.char="",stringsAsFactors = F,sep="")
  
}


combined_loci = NULL

for(samplename in samplenames){
  loci_file = paste0(outputDir,"/",samplename,"_loci.txt")
  print(loci_file)
  loci = read.table(loci_file,header=F,sep="\t",row.names=NULL,stringsAsFactors = F)
  print(nrow(loci))
  if(is.null(combined_loci)){
    combined_loci = loci
  }else{
    combined_loci = rbind(combined_loci,loci)
  }
}
combined_loci = combined_loci[!duplicated(combined_loci),]
combined_loci[,1] = factor(combined_loci[,1],levels = paste0("chr",c(1:22,"X")))
combined_loci = combined_loci[order(combined_loci[,1],combined_loci[,2]),]
combined_loci_file = paste0(outputDir,"/",case,"_loci_combined.txt")
write.table(combined_loci,combined_loci_file,col.names = F, row.names = F, quote=F, sep="\t")



