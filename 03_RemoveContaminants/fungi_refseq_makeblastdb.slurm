#!/bin/bash

#SBATCH --job-name=fungi_refseq_makeblastdb
#SBATCH --partition=comp01
#SBATCH --output=fungi_refseq_makeblastdb_%j.txt
#SBATCH --error=fungi_refseq_makeblastdb_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=01:00:00

module load java/sunjdk_1.8.0
module load python/3.7.3-anaconda
source /share/apps/bin/conda-3.7.3.sh
conda activate bioconda3-el7
module load blast/2.11.0 # I think this is the only one I really need, but loading the others first just in case

cd /local_scratch/$SLURM_JOB_ID/


############ LOCATIONS AND DIR NAMES ############
# location of the fungi.fasta from RefSeq 
# /scrfs/storage/amatthews/20220505_fungi_refseq/fna_files


## PROJECT DIRECTORY 
PROJECT_DIR=/scrfs/storage/amatthews/20220505_fungi_refseq

## BACTERIA DIR FOR BACTERIAL REFERENCES FROM REFSEQ NCBI DOWNLOAD
FUNGI_DIR=${PROJECT_DIR}/fna_files

mkdir -p ${FUNGI_DIR}/databases
mkdir ${FUNGI_DIR}/databases/fungi_refseq_reference

############ MAKE BLASTDB ############

makeblastdb -in ${FUNGI_DIR}/fungi_refseq.fa -title reference -dbtype nucl -out ${FUNGI_DIR}/databases/fungi_refseq_reference



