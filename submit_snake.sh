#!/bin/bash

#PBS -l nodes=1:ppn=11
#PBS -l gres=localhd:20
#PBS -l mem=80gb
#PBS -l vmem=80gb
#PBS -l walltime=80:00:00

module load samtools/1.5
module load star/2.6.1c
module load java/1.8.0_161
module load fastqc/0.11.5
module load rsem/1.2.22
module load R

cd $outdir

snakemake -s Snakefile --cores 11 --jobs 10 --configfile=config.yaml -k -r

