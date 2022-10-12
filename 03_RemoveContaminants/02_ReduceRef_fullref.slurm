#!/bin/bash

#SBATCH --job-name=02_ReduceRef_fullref
#SBATCH --partition=comp01
#SBATCH --output=02_ReduceRef_fullref_%j.txt
#SBATCH --error=02_ReduceRef_fullref_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=01:00:00

module load java/sunjdk_1.8.0
module load bbmap/38.82

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
## REDUCE THE REFERENCE TO CONTIGS >= 5000BP

reformat.sh in=${REF_DIR}/${REF_FULL}.fasta out=${REF_DIR}/${REF_FULL}_reduced.fasta minlength=5000
