#!/bin/bash

#$ -cwd -V
#$ -q long.qc
#$ -N runDP -j y 

./runDP_new_nonzero_NAP.sh
