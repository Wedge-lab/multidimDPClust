# multidimDPClust

This set of code allows the DPClust analysis to be undertaken at a multidimensional scale (i.e. multi-region sampling)

## STEP ONE

The first step is to create the multidimensional input files for DPClust.

The R code for this step is **GetDirichletProcessInfo_multipleSamples_updated_NAP.R** which can be run like

Rscript GetDirichletProcessInfo_multipleSamples_updated_NAP.R

This code requires **interconvertMutationBurdens.R** as input which is sourced at the beginning.

For this to run you need to produce a number of files as input:

### 1) samplenames

This is a file with the names of the samples for a particular patient:

Patient1a\
Patient1b\
Patient1c

### 2) cellularity

This is a file with the cellularity/purity values for the samples:

0.76\
0.83\
0.86

### 3) subclone.files

This file contains the name (including path if required) of the subclones.txt files obtained from Battenberg

Patient1a.subclones.txt\
Patient1b.subclones.txt\
Patient1c.subclones.txt

### 4) mutCount and WTCount

These files contain the raw counts of ref (WTCount) and alt (mutCount) alleles at each SNV to be used for mutation clustering with DPClust.

ID	Patient1a	Patient1b	Patient1c\
1_14000	7		8		0\
2_23400 0		10		9\
14_1230	12		0		0\
17_9876	11		9		13

where ID is CHROM_POSITION of SNVs

*IMPORTANT NOTE: remove any SNV which is absent in all samples (i.e. 0 in all samples)*

To generate these files the following steps need to be taken:

i) get union list of SNVs across the multiple samples from their respective VCF files

ii) Run **alleleCounter** on the BAM files of the multiple samples for the SNVs in the union loci list. The BAM index file (.bai) is required - if not present, run **samtools index BAMfile** to get it.

iii) match counts from alleleCounter with the SNV alleles of the union loci list and generate the WTCount and mutCount files.

*NOTE: please contact me if you have any difficulty with steps i)-iii).*

## STEP TWO

The second step is to run the DP clustering on the samples.

**submitDP_NAP.sh** submits the bash code **runDP_new_nonzero_NAP.sh** which in turn calls and executes **runDP_new_nonzero_NAP.R**

The runDP R code requires **subclone_Dirichlet_Gibbs_sampler_nD_binomial_DCW.R** and **plotnD.DCW.R** as input to reassign functions in the DPClust package for the purposes of running multidimensional DP and are sourced at the beginning.

## STEP THREE

The third step is to run the assignment of mutations to clusters per triplet combination of samples (when no. of samples > 3).

The R code for this step is **multidimclustparallel_NAP.R** which is called and executed by the bash code **multidimclustparallel_NAP.sh**. To submit the bash code to SGE, the number of parallel jobs in the jobarray must be calculated.
For example in R: choose(N,3) will give the number of triplet combinations where N is the number of samples obtained from a patient. If we have e.g. 11 samples, we will have choose(11,3)=165 triplet combinations.

The submission command would look like this:

**qsub -q QUEUE_NAME -cwd -t 1-165 -o /dev/null -e /dev/null multidimclustparallel_NAP.sh**


## STEP FOUR

The fourth and final step is to combine the clustering assignments across all triplets.

The R code for this is **combineClusteringParallelByAssignment_NAP.R** which can be run like

Rscript combineClusteringParallelByAssignment_NAP.R
