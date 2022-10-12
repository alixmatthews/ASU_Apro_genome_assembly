#!/bin/bash

#SBATCH --job-name=Index_PROW981_NCBI
#SBATCH --partition=comp01
#SBATCH --output=Index_PROW981_NCBI_%j.txt
#SBATCH --error=Index_PROW981_NCBI_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:10:00

module load bwa/0.7.17
module load samtools/1.15.1
module load java/sunjdk_1.8.0 # for gatk and picard
module load gatk/4.2.6.1
module load fastqc/0.11.5 # workaround for picard-tools
module load picard-tools/2.17.10



cd /local_scratch/$SLURM_JOB_ID/


##----------------------------------------------------------------------------------------
## MAKE AND ASSIGN DIRECTORIES

## MAIN DIR
PROJECT_DIR=/scrfs/storage/amatthews/20210816_projects/20221012_finalgenome

## SNP DIR (mkdir beforehand)
SNP_DIR=$PROJECT_DIR/IndexRef

## REF DIR (mkdir beforehand)
REF_DIR=$SNP_DIR/PROW981R1_NCBI

## REFERENCE GENOME
REF_FULL=scaffolds_reduced_contigs_kept_NCBI-200000000


##----------------------------------------------------------------------------------------
## INDEX, .FAI, AND .DICT THE REFERENCES
## only have to do this once per reference
## REF_FULL
bwa index ${REF_DIR}/${REF_FULL}.fsa

## ensure .fai file is created
samtools faidx ${REF_DIR}/${REF_FULL}.fsa

## create sequence dictionary for reference
java -jar /share/apps/bioinformatics/picard/picard-tools-2.17.10/picard.jar CreateSequenceDictionary REFERENCE=${REF_DIR}/${REF_FULL}.fsa OUTPUT=${REF_DIR}/${REF_FULL}.dict


