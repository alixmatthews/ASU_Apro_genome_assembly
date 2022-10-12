#!/bin/bash

#SBATCH --job-name=fungi_refseq_reduced_genome_blastn
#SBATCH --partition=comp01
#SBATCH --output=fungi_refseq_reduced_genome_blastn_%j.txt
#SBATCH --error=fungi_refseq_reduced_genome_blastn_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:30:00

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

## FUNGI DIR FOR FUNGI REFERENCES FROM REFSEQ NCBI DOWNLOAD
FUNGI_DIR=${PROJECT_DIR}/fna_files




############ MAKE BLASTDB - already done ############

# makeblastdb -in ${FUNGI_DIR}/fungi_refseq.fa -title reference -dbtype nucl -out ${FUNGI_DIR}/databases/fungi_refseq_reference



############ RUN BLASTN AGAINST FUNGI_REFSEQ_BLASTDB ############

for file in /scrfs/storage/amatthews/20210816_projects/20210816_snp/02_IndexRef/ref_full/scaffolds_reduced_split/*.fa;
do
   blastn -db ${FUNGI_DIR}/databases/fungi_refseq_reference -query ${file} -evalue 1e-3  -num_threads 31 -max_target_seqs 10 -max_hsps 10 -outfmt "7 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore" > ${file%.*}.txt;
done


