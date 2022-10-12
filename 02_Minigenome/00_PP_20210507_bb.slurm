#!/bin/bash

#SBATCH --job-name=bbduk_trim
#SBATCH --partition=comp06
#SBATCH --output=bbduk_trim_%j.txt
#SBATCH --error=bbduk_trim_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=06:00:00

module load java/sunjdk_1.8.0
module load fastqc/0.11.5
module load bbmap/38.82

module load python/3.7.3-anaconda
source /share/apps/bin/conda-3.7.3.sh
conda activate bioconda3-el7

cd /local_scratch/$SLURM_JOB_ID/

############ CONFIGURATION INFORMATION ############

############ LOCATIONS AND SAMPLE NAMES ############

## SAMPLE FILE NAMES 
## value1 is sample name
filenames="/scrfs/storage/amatthews/20210412_snp/00_PP_20210427/filenames.txt"

## PROJECT DIRECTORY 
PROJECT_SAMPLE_DIR=/scrfs/storage/amatthews/20210412_snp

## WHERE THE RAW DATA ARE
RAW_DATA=${PROJECT_SAMPLE_DIR}/raw_data


## ASSIGN 00_PP FOLDER TO QTRIM_DATA
QTRIM_DATA=${PROJECT_SAMPLE_DIR}/00_PP_20210427 


## MAKE NEW DIRECTORIES
## ADAPTER REMOVED
ADP_REM_BB=$QTRIM_DATA/Adapter_Removed_bb 
mkdir $ADP_REM_BB 

## FASTQC AFTER QUALITY TRIM
A_Q_T_BB=$QTRIM_DATA/FASTQC/FASTQ_After_Q_Trim_bb
mkdir $A_Q_T_BB



############ SOFTWARE SETTINGS ############

## QUALITY STUFFS
FASTQ_TYPE="sanger"
QUALIY_THRESHOLD=30
MIN_LENGTH=75 ## can optimize this later if desired


## ADAPTER INFO THAT NEEDS TO BE REMOVED
ADAPTER_SEQ=/scrfs/storage/amatthews/20210412_snp/00_PP_20210427/adapters.fa


############ START THE LOOP THROUGH FILES ############

while IFS=" " read -r value1
do {


############ ASSIGN FILES TO OBJECTS ############

FIRST_SAMPLE_LOC=${QTRIM_DATA}/Adapter_Removed/${value1}/${value1}_trim_R1.fastq
SECOND_SAMPLE_LOC=${QTRIM_DATA}/Adapter_Removed/${value1}/${value1}_trim_R2.fastq


############ ADAPTER REMOVING ############

mkdir $ADP_REM_BB/${value1}
cd $ADP_REM_BB/${value1}

bbduk.sh in1=$FIRST_SAMPLE_LOC in2=$SECOND_SAMPLE_LOC out1=${value1}_bb_trim_R1.fastq out2=${value1}_bb_trim_R2.fastq ref=$ADAPTER_SEQ ktrim=r k=23 mink=11 hdist=1 ftr=135 tpe tbo qtrim=rl trimq=$QUALIY_THRESHOLD minlen=$MIN_LENGTH


############ FASTQC AFTER TRIMMING ############

mkdir $A_Q_T_BB/${value1}
cd $A_Q_T_BB/${value1}

fastqc $ADP_REM_BB/${value1}/${value1}_bb_trim_R1.fastq $ADP_REM_BB/${value1}/${value1}_bb_trim_R2.fastq -o $PWD


} done <"$filenames"

############ MULTIQC ############

## run multiqc - just give it the path to fastqc files
multiqc /scrfs/storage/amatthews/20210412_snp/00_PP_20210427/FASTQC/FASTQ_After_Q_Trim_bb
