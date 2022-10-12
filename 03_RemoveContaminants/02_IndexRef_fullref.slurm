#!/bin/bash

#SBATCH --job-name=02_IndexRef_fullref
#SBATCH --partition=comp01
#SBATCH --output=02_IndexRef_%j.txt
#SBATCH --error=02_IndexRef_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00

##----------------------------------------------------------------------------------------
## LOAD MODULES

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

## FASTQ LOCATIONS
FASTQ=$PROJECT_DIR/01_DataCombo_20210903/Adapter_Removed_bb

## SNP DIR 
SNP_DIR=$PROJECT_DIR/02_IndexRef

## REF DIR (new subdir for the full genome)
REF_DIR=$SNP_DIR/ref_full
mkdir $REF_DIR

## COPY FULL REF OVER TO REF_DIR
cp /scrfs/storage/amatthews/20210816_projects/20210816_exp/01_ASSEMBLY/spades_20211014/PROW_981_R1_TGCTGTGA-GAATCGTG/scaffolds.fasta $REF_DIR


## REFERENCE GENOME
REF_FULL=scaffolds


##----------------------------------------------------------------------------------------
## INDEX, .FAI, AND .DICT THE REFERENCES
## only have to do this once per reference
## REF_FULL
bwa index $REF_DIR/$REF_FULL.fasta

## ensure .fai file is created
samtools faidx $REF_DIR/$REF_FULL.fasta

## create sequence dictionary for reference
java -jar /share/apps/bioinformatics/picard/picard-tools-2.17.10/picard.jar CreateSequenceDictionary REFERENCE=$REF_DIR/$REF_FULL.fasta OUTPUT=$REF_DIR/$REF_FULL.dict



