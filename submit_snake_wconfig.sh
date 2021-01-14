#!/bin/bash

#PBS -l nodes=1:ppn=10
#PBS -l gres=localhd:20
#PBS -l mem=100gb
#PBS -l vmem=100gb
#PBS -l walltime=80:00:00

module load samtools/1.5
module load star/2.6.1c
module load java/1.8.0_161
module load fastqc/0.11.5
module load rsem/1.2.22
module load bedtools/2.21.0
module load R/3.5.1
module load python/3.5.6_vcf2db
module load bcftools/1.6
module load vep/97
module load vt/fb0288b
module load bcftools
module load tabix

cd /hpf/largeprojects/ccmbio/yliang/test_place/test_snakemake

snakefile=main_test.snake
config_file=config.yaml
sample_file=test_samples.txt

snakemake -s $snakefile --jobs 1 --config sample_file=$sample_file --configfile=$config_file -k -r






