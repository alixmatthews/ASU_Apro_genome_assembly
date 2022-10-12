#!/bin/bash

#SBATCH --job-name=02_BUSCO_PROW981reduced_arachnida_exp
#SBATCH --partition=comp01
#SBATCH --output=02_BUSCO_PROW981reduced_arachnida_exp_%j.txt
#SBATCH --error=02_BUSCO_PROW981reduced_arachnida_exp_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:30:00

##----------------------------------------------------------------------------------------
## LOAD MODULES

module load python/anaconda-3.8
source /share/apps/bin/conda-3.8.sh
conda activate BUSCO

# I did this (below) to create a BUSCO environment with 5.2.2 in the home login node beforehand
# module load python/anaconda-3.8
# source /share/apps/bin/conda-3.8.sh
# mamba create -n BUSCO -c conda-forge -c bioconda busco=5.2.2


cd /local_scratch/$SLURM_JOB_ID/


##----------------------------------------------------------------------------------------
## "ADJUSTABLES"/DIRECTORIES

## PROJECT DIRECTORY 
PROJECT_DIR=/scrfs/storage/amatthews/20210816_projects


### ASSIGN OUTPUT DIR ###
OUTPUT_DIR=${PROJECT_DIR}/20210816_exp/02_BUSCO



############ SAMPLE(S) ARRAY ############

array=(

PROW_981_R1_TGCTGTGA-GAATCGTG_reduced
                      
)



############ RUN BUSCO ############



for file_name in "${array[@]}"

do {
	
	mkdir -p ${OUTPUT_DIR}
	mkdir ${OUTPUT_DIR}/PROW_981_R1_reduced
	mkdir ${OUTPUT_DIR}/PROW_981_R1_reduced/${file_name}
	mkdir ${OUTPUT_DIR}/PROW_981_R1_reduced/${file_name}/busco_output
	
	BUSCO_DIR=${OUTPUT_DIR}/PROW_981_R1_reduced/${file_name}/busco_output
	
	cd ${BUSCO_DIR}
	
	busco --mode genome --cpu 32 --force --in ${PROJECT_DIR}/20210816_snp/02_IndexRef/ref_full/scaffolds_reduced.fasta --out ${file_name} --lineage_dataset arachnida_odb10

}


done












