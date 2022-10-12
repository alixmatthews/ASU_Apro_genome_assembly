#!/bin/bash

#SBATCH --job-name=02_IndexRef_reducedref
#SBATCH --partition=comp01
#SBATCH --output=02_IndexRef_reducedref_%j.txt
#SBATCH --error=02_IndexRef_reducedref_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:10:00

module load bwa/0.7.17
module load samtools/1.10
module load java/sunjdk_1.8.0 # for gatk and picard
module load gatk/4.0.0
module load fastqc/0.11.5 # workaround for picard-tools
module load picard-tools/2.17.10

cd /local_scratch/$SLURM_JOB_ID/


##----------------------------------------------------------------------------------------
## MAKE AND ASSIGN DIRECTORIES

## MAIN DIR
PROJECT_DIR=/scrfs/storage/amatthews/20210816_projects/20210816_snp

## SNP DIR 
SNP_DIR=$PROJECT_DIR/02_IndexRef

## REF DIR (new subdir for the full genome)
REF_DIR=$SNP_DIR/ref_full

## REFERENCE GENOME
REF_FULL=scaffolds


##----------------------------------------------------------------------------------------
## INDEX, .FAI, AND .DICT THE REFERENCES
## only have to do this once per reference
## REF_FULL
bwa index ${REF_DIR}/${REF_FULL}_reduced.fasta

## ensure .fai file is created
samtools faidx ${REF_DIR}/${REF_FULL}_reduced.fasta

## create sequence dictionary for reference
java -jar /share/apps/bioinformatics/picard/picard-tools-2.17.10/picard.jar CreateSequenceDictionary REFERENCE=${REF_DIR}/${REF_FULL}_reduced.fasta OUTPUT=${REF_DIR}/${REF_FULL}_reduced.dict


