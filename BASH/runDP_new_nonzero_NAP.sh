#!/bin/bash
module load R/3.3.0
R CMD BATCH --no-save --no-restore runDP_new_nonzero_NAP.R runDP_new_nonzero_NAP.Rout
