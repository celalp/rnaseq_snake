#!/bin/bash
# properties = {properties}

unset R_LIBS
module purge
module load samtools/1.5
module load star/2.6.1c
module load java/1.8.0_161
module load fastqc/0.11.5
module load rsem/1.2.22
module load R/3.5.1
module load python/3.5.6_vcf2db


{exec_job}