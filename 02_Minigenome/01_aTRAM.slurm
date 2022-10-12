#!/bin/bash

#SBATCH --job-name=aTRAM
#SBATCH --partition=comp72
#SBATCH --output=aTRAM_%j.txt
#SBATCH --error=aTRAM_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00

module load exonerate/2.4.0
module load python/3.7.3-anaconda
source /share/apps/bin/conda-3.7.3.sh
conda activate aTRAM

cd /local_scratch/$SLURM_JOB_ID/

############################################################
############ CONFIGURATION INFORMATION ############
############ LOCATIONS AND SAMPLE NAMES ############
############################################################

## SAMPLE FILE NAMES 
filenames="/scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/atram_spp.txt" #this file will only need to include one individual sample from each mite species (one reference per species)

## MAIN DIR
PROJECT_DIR=/scrfs/storage/amatthews/20210412_snp

## ATRAM DIR
## contains mt reference (Proctophyllodes_miliariae) from Dona et al. in amino acids; 13 PCGs
## contains nu reference (scabies mite) from OrthoDB Acari in amino acids; 1692 genes
ATRAM_DIR=$PROJECT_DIR/01_aTRAM_20210507
mkdir ATRAM_DIR

## ATRAM DATABASES
ATRAM_DB=$ATRAM_DIR/aTRAM_db
mkdir ATRAM_DB


############ SET UP PER SAMPLE DIRS ############

while IFS=" " read -r value1
do {


## FASTQ FILE DIRECTORY
SAMPLE_FASTQ=$PROJECT_DIR/00_PP_20210427/Adapter_Removed_bb/${value1}

# ATRAM_DB DIR
# note to self: may need another /{value1} see --blast-db argument for atram_preprocessor.py and atram.py commands
mkdir $ATRAM_DB/{value1}

# SAMPLE OUTPUT / (NU, MT) / (ALL, FILTERED, STITCH)
SAMPLE_OUT=$ATRAM_DIR/${value1}_out
mkdir $SAMPLE_OUT
SAMPLE_OUT_MT=$SAMPLE_OUT/${value1}_mt
mkdir -p $SAMPLE_OUT_MT
mkdir $SAMPLE_OUT_MT/all
mkdir $SAMPLE_OUT_MT/filtered
mkdir $SAMPLE_OUT_MT/stitch
SAMPLE_OUT_NU=$SAMPLE_OUT/${value1}_nu
mkdir -p $SAMPLE_OUT_NU
mkdir $SAMPLE_OUT_NU/all
mkdir $SAMPLE_OUT_NU/filtered
mkdir $SAMPLE_OUT_NU/stitch


############################################################
############ RUN aTRAM ############
############################################################

############ PREPROCESSOR ############

atram_preprocessor.py --blast-db=$ATRAM_DB/${value1}/${value1} --end-1=$SAMPLE_FASTQ/${value1}_bb_trim_R1.fastq  --end-2=$SAMPLE_FASTQ/${value1}_bb_trim_R2.fastq
## second /${value1} on --blast-db line is referring to the prefix to use for the file names it creates in this step


############ ASSEMBLER ############
############ MITOCHONDRIAL ############

atram.py --blast-db=$ATRAM_DB/${value1}/${value1} --query-split=$ATRAM_DIR/mt_ref.fasta --assembler=spades --protein --db-gencode=5 --output-prefix=$SAMPLE_OUT_MT/all/
## second /${value1} on --blast-db line is referring to the prefix to use/search for the file names it created in the previous step

cd $SAMPLE_OUT_MT
find . -name "*filtered_contigs.fasta" -exec mv -i {} -t "./filtered/" ";"
## move filtered_contigs files to a new folder called "filtered" (separate them from non-filtered)

############ NUCLEAR ############
atram.py --blast-db=$ATRAM_DB/${value1}/${value1} --query-split=$ATRAM_DIR/nu_ref.fasta --assembler=spades --protein --db-gencode=5 --output-prefix=$SAMPLE_OUT_NU/all/
## second /${value1} on --blast-db line is referring to the prefix to use/search for the file names it created in the previous step

cd $SAMPLE_OUT_NU
find . -name "*filtered_contigs.fasta" -exec mv -i {} -t "./filtered/" ";"
## move filtered_contigs files to a new folder called "filtered" (separate them from non-filtered)


############ STITCHER ############
############ MITOCHONDRIAL ############

atram_stitcher.py --assemblies-dir=$SAMPLE_OUT_MT/filtered --reference-genes=$ATRAM_DIR/mt_ref.fasta --taxa=/scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/atram_spp.txt --output-prefix=$SAMPLE_OUT_MT/stitch/${value1}

############ NUCLEAR ############

atram_stitcher.py --assemblies-dir=$SAMPLE_OUT_NU/filtered --reference-genes=$ATRAM_DIR/nu_ref.fasta --taxa=/scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/atram_spp.txt --output-prefix=$SAMPLE_OUT_NU/stitch/${value1}


############ FINISH aTRAM ############

} done <"$filenames"

