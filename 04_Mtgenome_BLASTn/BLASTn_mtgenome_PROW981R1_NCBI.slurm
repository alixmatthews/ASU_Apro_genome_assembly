#!/bin/bash

#SBATCH --job-name=BLASTn_mtgenome_PROW981R1_NCBI
#SBATCH --partition=comp01
#SBATCH --output=BLASTn_mtgenome_PROW981R1_NCBI_%j.txt
#SBATCH --error=BLASTn_mtgenome_PROW981R1_NCBI_%j.err
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

## PROJECT DIRECTORY 
PROJECT_DIR=/scrfs/storage/amatthews/20210816_projects/20221012_finalgenome

## BLASTn dir (mkdir beforehand)
ASSEMBLY_DIR=${PROJECT_DIR}/BLASTn


############ MAKE BLASTDB ############

makeblastdb -in /scrfs/storage/amatthews/20210816_projects/20221012_finalgenome/IndexRef/PROW981R1_NCBI/scaffolds_reduced_contigs_kept_NCBI-200000000.fsa -title reference -dbtype nucl -out ${ASSEMBLY_DIR}/databases/PROW981R1_NCBI


############ RUN BLASTN AGAINST MITOGENOME ############

blastn -db ${ASSEMBLY_DIR}/databases/PROW981R1_NCBI -query ${ASSEMBLY_DIR}/Trubecula_mtgenome.fasta -evalue 1e-3  -outfmt "7 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore" > ${ASSEMBLY_DIR}/PROW981R1_NCBI_blastout.txt



