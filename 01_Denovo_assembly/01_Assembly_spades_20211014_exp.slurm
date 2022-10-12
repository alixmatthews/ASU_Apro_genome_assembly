#!/bin/bash

#SBATCH --job-name=01_Assembly_spades_20211014_exp
#SBATCH --partition=himem72
#SBATCH --output=01_Assembly_spades_20211014_exp_%j.txt
#SBATCH --error=01_Assembly_spades_20211014_exp_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=12:00:00

module load spades/3.12.0 # I think this is the only one I really need, but loading the others just in case
module load java/sunjdk_1.8.0
module load python/3.7.3-anaconda
source /share/apps/bin/conda-3.7.3.sh
conda activate bioconda3-el7


cd /local_scratch/$SLURM_JOB_ID/


############ CONFIGURATION INFORMATION ############

############ LOCATIONS AND DIR NAMES ############

## PROJECT DIRECTORY 
PROJECT_DIR=/scrfs/storage/amatthews/20210816_projects/20210816_exp


## FASTQ FILE DIRECTORY
SAMPLE_FASTQ=${PROJECT_DIR}/00_PP_20210903/Adapter_Removed_bb


## ASSIGN AND MAKE ASSEMBLY DIRECTORY (mkdir beforehand)
ASSEM_DIR=${PROJECT_DIR}/01_ASSEMBLY


## MAKE ASSEMBLY DIRECTORY FOR TODAY'S ATTEMPT (ASSEMBLYPROGRAM_DATE)
SPADES_DIR=${ASSEM_DIR}/spades_20211014
mkdir ${SPADES_DIR}

OUTPUT_DIR=${SPADES_DIR}




############ SAMPLE(S) ARRAY ############

array=(

PROW_953_R1_TGACAACC-CTGTTGAC
PROW_954_R1_CTGACTAC-ACATTGCG
PROW_981_R1_TGCTGTGA-GAATCGTG
                      
)



############ RUN SPADES ############

# removed --careful flag, the other assembly was better without the --careful flag!

for file_name in "${array[@]}"

do {
	mkdir ${OUTPUT_DIR}/${file_name}
	cd ${PROJECT_DIR}/

	paired_end1=${SAMPLE_FASTQ}/${file_name}/${file_name}_dd_bb_trim_R1.fastq
	paired_end2=${SAMPLE_FASTQ}/${file_name}/${file_name}_dd_bb_trim_R2.fastq

spades.py -t 15 -k 21, 33, 55, 77, 99, 111 --only-assembler --pe1-1 $paired_end1 --pe1-2 $paired_end2 -o $OUTPUT_DIR/${file_name}  

}


done







