#!/bin/bash

#SBATCH --job-name=00_PP_20210906_b_exp
#SBATCH --partition=comp06
#SBATCH --output=00_PP_20210906_b_exp_%j.txt
#SBATCH --error=00_PP_20210906_b_exp_%j.err
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
filenames="/scrfs/storage/amatthews/20210816_projects/20210816_exp/20210816_exp_filenames.txt"

## PROJECT DIRECTORY 
PROJECT_SAMPLE_DIR=/scrfs/storage/amatthews/20210816_projects

## WHERE THE RAW DATA ARE
RAW_DATA=${PROJECT_SAMPLE_DIR}/raw_data

## ASSIGN 00_PP FOLDER TO QTRIM_DATA
QTRIM_DATA=${PROJECT_SAMPLE_DIR}/20210816_exp/00_PP_20210903

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
ADAPTER_SEQ=/scrfs/storage/amatthews/20210816_projects/20210816_exp/00_PP_20210903/adapters.fa


############ START THE LOOP THROUGH FILES ############

while IFS=" " read -r value1
do {


############ ASSIGN FILES TO OBJECTS ############

FIRST_SAMPLE_LOC=${QTRIM_DATA}/Adapter_Removed/${value1}/${value1}_trim_R1.fastq
SECOND_SAMPLE_LOC=${QTRIM_DATA}/Adapter_Removed/${value1}/${value1}_trim_R2.fastq


############ ADAPTER REMOVING ############

mkdir $ADP_REM_BB/${value1}
cd $ADP_REM_BB/${value1}

bbduk.sh in1=$FIRST_SAMPLE_LOC in2=$SECOND_SAMPLE_LOC out1=${value1}_bb_trim_R1.fastq out2=${value1}_bb_trim_R2.fastq ref=$ADAPTER_SEQ ktrim=r k=23 mink=11 hdist=1 ftr=140 tpe tbo qtrim=rl trimq=$QUALIY_THRESHOLD minlen=$MIN_LENGTH

# "ktrim=r” is for right-trimming (3′ adapters)
# “k” specifies the kmer size to use (23)
# “mink” allows it to use shorter kmers at the ends of the read (for example, k=11 for the last 11 bases)
# “hdist” means “hamming distance”; this allows 1 mismatch
# "frt" means trim the right end, retaining the value and lower (ftr=140 means keep positions 140 and lower)
# For normal paired-end fragment libraries, I (bbduk folks) recommend adding the flags “tbo”, which specifies to also trim adapters based on pair overlap detection using BBMerge (which does not require known adapter sequences), and “tpe”, which specifies to trim both reads to the same length (in the event that an adapter kmer was only detected in one of them).
# “qtrim=rl” quality-trim to Q10 using the Phred algorithm on both sides




############ DEDUPLICATING ############

dedupe.sh in1=$ADP_REM_BB/${value1}/${value1}_bb_trim_R1.fastq in2=$ADP_REM_BB/${value1}/${value1}_bb_trim_R2.fastq out=$ADP_REM_BB/${value1}/${value1}_bb_trim_dd.fastq ac=false

# ac=false sets absorbContainment=false (which is necessary for paired reads)
# naming this output *_trim_dd so that my reformatted output can be *_trim_R1/R2 (to match endings with the other outputs that were not deduplicated)

reformat.sh in=$ADP_REM_BB/${value1}/${value1}_bb_trim_dd.fastq out1=$ADP_REM_BB/${value1}/${value1}_dd_bb_trim_R1.fastq out2=$ADP_REM_BB/${value1}/${value1}_dd_bb_trim_R2.fastq

# so now when calling deduplicated samples, will need to add "_dd" to the end of their file names 




############ FASTQC AFTER TRIMMING AND DD'ING ############

mkdir $A_Q_T_BB/${value1}
cd $A_Q_T_BB/${value1}

# run fastqc on these dd'd samples
fastqc $ADP_REM_BB/${value1}/${value1}_dd_bb_trim_R1.fastq  $ADP_REM_BB/${value1}/${value1}_dd_bb_trim_R2.fastq -o $PWD



} done <"$filenames"

############ MULTIQC ############

## run multiqc - just give it the path to fastqc files
multiqc /scrfs/storage/amatthews/20210816_projects/20210816_exp/00_PP_20210903/FASTQC/FASTQ_After_Q_Trim_bb



