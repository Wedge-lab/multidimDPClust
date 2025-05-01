#!/bin/bash


# The value of run_in_parelle is the answer of this question:
# How many ways can we choose 3 samples out of X samples without considering the order? (Combinatorial calculation)
# for example: if you have 8 samples, the value of run_in_parelle is 56

# for ((i=4;i<=12;i++));do sh pipeline.sh $i;done

patient_number=$1



########################################################################################
#########Please modify this part accroding to your own run_info.txt file################
########################################################################################

# This script could be run as : sh pipeline.sh 1
#  1 is the patient_number here

if [ "$patient_number" -eq 1 ]; then
	case=DW111230_original
	sample_num=3
	run_in_parallel=1

elif [ "$patient_number" -eq 2 ]; then
	case=multi_case_reseq_refit_WGD
	sample_num=8
	run_in_parallel=56
else
    echo "Invalid patient number. Please provide 1 or 2."
    exit 1
fi

echo "Case: $case"
echo "Sample number: $sample_num"
echo "Run in Parallel: $run_in_parallel"

########################################################################################
################## Please change the file path into your own file path##################
########################################################################################
Bam_file_path=/mnt/bmh01-rds/UoOxford_David_W/r44399mg/CRPM/Bam_New/no_hardcut_data_remove_DUP/
snv_path=/mnt/bmh01-rds/UoOxford_David_W/r44399mg/CRPM/SNVs/Consensus/
BB_file_path=/mnt/bmh01-rds/UoOxford_David_W/r44399mg/CRPM/BB2210/
run_info_path=/mnt/bmh01-rds/UoOxford_David_W/r44399mg/CRPM/Code/msDPClust_seperate_12/BASH/
outpath=/mnt/bmh01-rds/UoOxford_David_W/r44399mg/CRPM/msDPCtest

ab=`qsub -terse 00_run.txt ${case} ${snv_path} ${run_info_path} ${outpath} `
ab=${ab%%.*}
echo $ab


cd=`qsub -hold_jid ${ab} -terse -t 1-${sample_num} 01_allele_count.txt ${case} ${Bam_file_path} ${outpath} `
cd=${cd%%.*}
echo $cd



ef=`qsub -hold_jid ${cd} -terse 02_run.txt ${case} ${run_info_path} ${outpath} `
ef=${ef%%.*}
echo $ef


hj=`qsub -hold_jid ${ef} -terse 03_Get_DP_Input.txt ${case} ${run_info_path} ${BB_file_path} ${outpath} `
hj=${hj%%.*}
echo $hj


hjg=`qsub -hold_jid ${hj} -terse 04_runDP_GS.txt ${case} ${outpath} `
hjg=${hjg%%.*}
echo $hjg


hjl=`qsub -hold_jid ${hjg} -terse -t 1-${run_in_parallel} 05_parallel_assign_mutations.txt ${case} ${sample_num} ${outpath} `
hjl=${hjl%%.*}
echo $hjl

hjk=`qsub -hold_jid ${hjl} -terse 06_combine.txt ${case} ${outpath} `
echo $hjk
echo submit done
