#!/bin/bash

#SBATCH --job-name=BUSCO_PROW981R1_NCBI_arachnida
#SBATCH --partition=comp01
#SBATCH --output=BUSCO_PROW981R1_NCBI_arachnida_%j.txt
#SBATCH --error=BUSCO_PROW981R1_NCBI_arachnida_%j.err
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


### ASSIGN OUTPUT DIR ### (mkdir beforehand)
OUTPUT_DIR=${PROJECT_DIR}/20221012_finalgenome/BUSCO



############ SAMPLE(S) ARRAY ############

array=(

scaffolds_reduced_contigs_kept_NCBI-200000000
                      
)



############ RUN BUSCO ############



for file_name in "${array[@]}"

do {
	
	mkdir -p ${OUTPUT_DIR}
	mkdir ${OUTPUT_DIR}/PROW981R1_NCBI
	mkdir ${OUTPUT_DIR}/PROW981R1_NCBI/${file_name}
	mkdir ${OUTPUT_DIR}/PROW981R1_NCBI/${file_name}/busco_output_arachnida
	
	BUSCO_DIR=${OUTPUT_DIR}/PROW981R1_NCBI/${file_name}/busco_output_arachnida
	
	cd ${BUSCO_DIR}
	
	busco --mode genome --cpu 32 --force --in ${PROJECT_DIR}/20221012_finalgenome/IndexRef/PROW981R1_NCBI/scaffolds_reduced_contigs_kept_NCBI-200000000.fsa --out ${file_name} --lineage_dataset arachnida_odb10

}


done












