#!/bin/bash

#PBS -l nodes=1:ppn=10
#PBS -l gres=localhd:20
#PBS -l mem=16gb
#PBS -l vmem=16gb
#PBS -l walltime=05:00:00

module load samtools/1.5
module load star/2.6.1c
module load java/1.8.0_161
module load fastqc/0.11.5
module load rsem/1.2.22
module load R/4.0.3
module load python/3.7.6
module load tabix

cd /hpf/largeprojects/ccmbio/yliang/test_place/test_snakemake

snakefile=main_test.snake
config_file=config.yaml
sample_file=test_samples.txt

snakemake -s $snakefile --jobs 1 --config sample_file=$sample_file --configfile=$config_file -k -r






