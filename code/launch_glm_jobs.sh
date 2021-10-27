#!/usr/bin/env bash

cd /t1-data/project/psychdementia/mfernand/PROJECTS/STRS-WILL
mkdir results
cd ./results

module add R-cbrg/current

export RSCRIPT=/t1-data/project/psychdementia/mfernand/PROJECTS/STRS-WILL/code
export REPEATS_PATH=/t1-data/project/psychdementia/mfernand/PROJECTS/STRS-WILL/chr.list
export COVARIATES=/t1-data/project/psychdementia/shared/Repeats/Variables/Covariates_MRI_HIP_1stAttempt

for chr in $(cat $REPEATS_PATH)
do
sbatch -p batch -c 8 --job-name glm_${chr} -o %j.out -e %j.err --mail-user=mfernand --mail-type=ALL \
--wrap="Rscript --vanilla $RSCRIPT/f_001_glm_gaussian.R \"${chr}\" ${COVARIATES} 8"
done