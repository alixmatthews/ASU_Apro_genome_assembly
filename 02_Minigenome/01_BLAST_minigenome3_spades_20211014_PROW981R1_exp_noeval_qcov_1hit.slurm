#!/bin/bash

#SBATCH --job-name=01_BLAST_minigenome3_spades_20211014_PROW981R1_exp_noeval_qcov_1hit
#SBATCH --partition=comp01
#SBATCH --output=01_BLAST_minigenome3_spades_20211014_PROW981R1_exp_noeval_qcov_1hit_%j.txt
#SBATCH --error=01_BLAST_minigenome3_spades_20211014_PROW981R1_exp_noeval_qcov_1hit_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=00:05:00

module load java/sunjdk_1.8.0
module load python/3.7.3-anaconda
source /share/apps/bin/conda-3.7.3.sh
conda activate bioconda3-el7
module load blast/2.11.0 # I think this is the only one I really need, but loading the others first just in case

cd /local_scratch/$SLURM_JOB_ID/


############ LOCATIONS AND DIR NAMES ############

## PROJECT DIRECTORY 
PROJECT_DIR=/scrfs/storage/amatthews/20210816_projects/20210816_exp

## ASSEMBLY DIR FOR SPECIFIC ASSEMBLY METHOD/DAY/SAMPLE
ASSEMBLY_DIR=${PROJECT_DIR}/01_ASSEMBLY/spades_20211014/PROW_981_R1_TGCTGTGA-GAATCGTG/


############ MAKE BLASTDB - already done ############

# makeblastdb -in ${ASSEMBLY_DIR}/scaffolds.fasta -title reference -dbtype nucl -out ${ASSEMBLY_DIR}/databases/PROW_981_R1_reference


############ RUN BLASTN AGAINST MITOGENOME ############

blastn -db ${ASSEMBLY_DIR}/databases/PROW_981_R1_reference -query ${PROJECT_DIR}/01_ASSEMBLY/PROW_831_R_GATCTATC-AGCCTCAT_exons_mt95_best_nu95_ref_3.fasta -max_target_seqs 1 -max_hsps 1 -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" > ${ASSEMBLY_DIR}/PROW_981_R1_minigenome3_noeval_qcov_1hit_blastout.txt



