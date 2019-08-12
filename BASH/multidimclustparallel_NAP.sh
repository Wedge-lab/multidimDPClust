#!/bin/bash
module load R/3.3.0
R CMD BATCH '--no-restore --no-save --args 1 11 3 '$SGE_TASK_ID /well/wedge/naser/melanoma/exomes/multdim/multidimclustparallel_NAP.R naser_multidimclustparallel.$SGE_TASK_ID.Rout # 1 patient, 11 samples and run in triplets (3)
exit $?
