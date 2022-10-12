# Final *Amerodectes protonotaria* genome pipeline

- Steps needed from preprocessed reads obtained from sequencing facility to NCBI-deposited genome

---
---

### 00. Preprocessing reads

- This was done on all samples from this batch sent off for sequencing (1 mite, 5 mites, 20 mites). I call them the 'experimental' mites

#### Slurms (run in order):
- Need this file for both slurms: `20210816_exp_filenames.txt`
- `00_PP_20210903_a_exp.slurm`
- `00_PP_20210903_b_exp.slurm`
  - Need `adapters.fa` 

#### Output
- `00_PP_20210903_a_exp_multiqc_report.html`
- `00_PP_20210906_b_exp_multiqc_report.html`

---
---

### 01. *De novo* assembly

- This was done on a few samples, and only one lucky sample made it through to the next steps!

#### Slurm
- `01_Assembly_spades_20211014_exp.slurm`

---
---

### 02. Minigenome assembly

- We used PROW_831_R_GATCTATC-AGCCTCAT as the pooled minigenome assembly, but these files also include all the other samples that were sequenced at this same time from the `20210412_snp` sequencing run.
- This was before I had more cleanly named files! Bear with me!

---
#### 02.1. Preprocessing reads

##### Slurms (run in order):
- Need this file for both slurms: `filenames.txt`
- `00_PP_20210427.slurm`
- `00_PP_20210507_bb`
  - Need `adapters.fa` 

##### Output
- `00_PP_20210427_multiqc_report.html`
- `00_PP_20200507_multiqc_report.html`

---
#### 02.2. aTRAM

##### 02.2.1. OrthoDB
- First need to get the *Sarcoptes scabiei* mites orthologous sequences from OrthoDB

```
curl 'https://www.orthodb.org/fasta?universal=1&singlecopy=1&limit=5000&level=6933&species=52283_0' -o Acari_scabies_nuclear_orthologs.fasta
```

- Then also pulled the mitochondrial genome reference sequences for *Proctophyllodes miliariae* from published database
- Translate both of these into amino acid sequences

##### 02.2.2. Run aTRAM 

###### Slurms: 
- Need this `atram_spp.txt` (PROW 831 is the sample of interest)
- `01_aTRAM.slurm`

##### 02.2.3. Merge nu and mt data
- Then I ran a filter-merge-concatenate pipeline by hand. This is to retain genes that were >95% complete. It is outlined here: `02_FMC_manual.sh` 
- These files are needed: `PROW_831_R_GATCTATC-AGCCTCAT_mt95_ref.txt` and `PROW_831_R_GATCTATC-AGCCTCAT_nu95_ref.txt`

---
#### 02.3. BLASTn

- Had to split the minigenome into 3 pieces because of blast limits

##### Slurms
- `01_BLAST_minigenome1_spades_20211014_PROW981R1_exp_noeval_qcov_1hit.slurm`
- `01_BLAST_minigenome2_spades_20211014_PROW981R1_exp_noeval_qcov_1hit.slurm`
- `01_BLAST_minigenome3_spades_20211014_PROW981R1_exp_noeval_qcov_1hit.slurm`

##### Outputs
- `PROW_981_R1_minigenome1_noeval_qcov_1hit_blastout.txt`
- `PROW_981_R1_minigenome2_noeval_qcov_1hit_blastout.txt`
- `PROW_981_R1_minigenome3_noeval_qcov_1hit_blastout.txt`
- `PROW_981_R1_minigenome_combo_noeval_qcov_1hit_blastout.xlsx`

---
---

### 03. Identify and remove potential contaminants from SPAdes *de novo* assembly

#### 03.1. Index the SPAdes assembly

##### Slurm: 
- `02_IndexRef_fullref.slurm`
  - Take a look at scaffold lengths from the assembly:
  
  ```
  grep ">" scaffolds.fasta | cut -f4 -d"_" > PROW_981_scaffold_len.txt
  # output: PROW_981_scaffold_len.txt
  ```

- Switch to R and run the code below
- Figured out there were tons and tons of contigs (26,013), which was clogging up many of the HPC resources. Checked the number of contigs <1000 and <5000 in R, and then plotted how many reads mapped to the remaining contigs under each condition. We did not lose that many mapped reads when removing contigs <1000 or <5000 (and the diffs between species was not major).

```
library(ggplot2)
scaffolds<-read.csv("PROW_981_scaffold_len.txt", sep = ",", header = FALSE)
head(scaffolds)
tail(scaffolds)
total <- length(scaffolds$V1) # 26013 scaffolds

lt1000 <- length(which(scaffolds$V1 < 1000)) # 19597 scaffolds less than 1000bp
lt5000 <- length(which(scaffolds$V1 < 5000)) # 23433 scaffolds less than 5000bp

total-lt5000 # 2580 = number of contigs >5000bp
total-lt1000 # 6416 = number of contigs >1000bp
```

---

#### 03.2. Remove contigs <5000bp from the original reference

##### Slurms: 
- `02_ReduceRef_fullref.slurm` 
- `02_IndexRef_reducedref.slurm`

```
grep ">" scaffolds_reduced.fasta | wc -l
2580
```

---

#### 03.3. BUSCO on the reduced (i.e., without contigs <5000bp) genome

##### Slurms: 
- `02_BUSCO_PROW981reduced_arachnida_exp.slurm`
- `02_BUSCO_PROW981reduced_arthropoda_exp.slurm` 

##### Outputs: 
- `02_BUSCO_PROW981reduced_arachnida_exp_short_summary.txt`
- `02_BUSCO_PROW981reduced_arthropoda_exp_short_summary.txt`

---

#### 03.4. Remove microbes
- Now it's time to remove any and all microbes from the reduced genome. This is a long step that doesn't have many slurms associated, so writing all of the fun steps out here instead.

- I did the first three steps on my laptop/Linux machine to get the assembly_summary.txt file and parse the ftpdirpaths and ftpfilepaths from that .txt file. Then I put the outputs on a dir in my `scrfs/storage/amatthews` 

##### 03.4.1. Obtain Bacteria, Fungi, and Archaea RefSeq sequences from NCBI
- Download the `assembly_summary.txt` from RefSeq for bacteria, fungi, and archaea. 
- NOTE: Bacteria is the only one without a specialized prefix from here onward.

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt # this is bacteria, but is not renamed
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt # renamed to fungi_assembly_summary.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt # renamed to archaea_assembly_summary.txt
```


##### 03.4.2. List the FTP path (column 20) for the assemblies of interest
- In this case those that have "Complete Genome" assembly_level (column 12) and "latest" version_status (column 11)

```
awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $20}' assembly_summary.txt > ftpdirpaths
awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $20}' fungi_assembly_summary.txt > fungi_ftpdirpaths
awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $20}' archaea_assembly_summary.txt > archaea_ftpdirpaths
```

##### 03.4.3. Append the filename of interest
- In this case "*_genomic.fna.gz" to the FTP directory names
    
```
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths > ftpfilepaths
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' fungi_ftpdirpaths > fungi_ftpfilepaths
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' archaea_ftpdirpaths > archaea_ftpfilepaths
```


##### 03.4.4. More file organization

###### 03.4.4.1. Make separate directories for all the microbes on scrfs/storage/amatthews for easy accessibility

```
cd scrfs/storage/amatthews
mkdir 20220418_bacteria_refseq
mkdir 20220505_fungi_refseq
mkdir 20220505_archaea_refseq
```

###### 03.4.4.2. Put the ftpfilepaths (from 03.4.3.), ftpdirpaths (from 03.4.2.), and X_assembly_summary.txt (from 03.4.1.) in the appropriate directories using FileZilla

```
cd 20220418_bacteria_refseq
mkdir fna_files
cd fna_files
# do the same for archaea and fungi
```


##### 03.4.5. Start an srun
- There were some firewall issues on the HPC with just doing wget on its own so... have to open an srun on a cloud72 server

```
srun --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --partition cloud72 --qos cloud --time=72:00:00 --pty /bin/bash
# then the shell is open (bash) for me to run away with the code! 
# Don't forget to cd into the fna_files directory, though
```


##### 03.4.6. Download the datafile for each FTP path in the list
- This takes a while for bacteria because there are 25,999 files

```
wget -i ../ftpfilepaths
```

##### 03.4.7. Cat the files 
- This takes a little bit for bacteria (not as long as the wget!)
```
zcat *.gz > bac_refseq.fa
```

##### 03.4.8. Repeat 03.4.6. and 03.4.7. for Fungi and Archaea

``` 
## first for the Fungi 
cd ../..
cd 20220505_fungi_refseq
mkdir fna_files
cd fna_files
wget -i ../fungi_ftpfilepaths
zcat *.gz > fungi_refseq.fa # there are only 22 files so it doesn't take long

## then for the Archaea
cd ../..
cd 20220505_archaea_refseq
mkdir fna_files
cd fna_files
wget -i ../archaea_ftpfilepaths
zcat *.gz > archaea_refseq.fa # there are only 429 files so also doesn't take very long
```

##### 03.4.9. Make the local blast database 
###### Slurms:
- `bac_refseq_makeblastdb.slurm`
- `fungi_refseq_makeblastdb.slurm` 
- `archaea_refseq_makeblastdb.slurm`


##### 03.4.10. Split the .fasta file 
This is necessary for the BLASTn to work (because there are limits; https://ncbiinsights.ncbi.nlm.nih.gov/2020/06/18/new-blast-settings/)

```
# on AHPCC...
ssh pinnacle-l6
module load python/anaconda-3.8
source /share/apps/bin/conda-3.8.sh
mamba create -n faSplit ucsc-fasplit
conda activate faSplit

# pwd: /scrfs/storage/amatthews/20210816_projects/20210816_snp/02_IndexRef/ref_full
mkdir scaffolds_reduced_split
faSplit byname scaffolds_reduced.fasta ./scaffolds_reduced_split/

# can check to see if the number of nucleotides matches with the contig names for a few random ones
grep -v ">" NODE_987_length_15793_cov_1.813863.fa | wc | awk '{print $3-$1}'
grep -v ">" NODE_2153_length_6366_cov_1.965787.fa | wc | awk '{print $3-$1}'
grep -v ">" NODE_14_length_210216_cov_2.410195.fa | wc | awk '{print $3-$1}'

# all contigs are showing the appropriate amount of nucleotides! Good to go...
```

##### 03.4.11. Run BLASTn
- Be sure to remove all the NODE_X.txt files from the previous slurm before running the next slurm (`rm NODE_*txt`); 
  - otherwise it will overwrite the files. So also have to run these slurms one at a time. Could potentially make a new directory for each slurm's output to avoid having to wait in sequence, but waiting is fine and gives some time to do other things... :) 

###### Slurms: 

**NOTE: DON'T FORGET TO `rm NODE_*txt` IN BETWEEN EACH SLURM, DUE TO THE WAY THIS IS SET UP!**
- `bac_refseq_reduced_genome_blastn.slurm`
   - original, before seeing to change the number of threads, https://voorloopnul.com/blog/how-to-correctly-speed-up-blast-using-num_threads/)
- `bac_refseq_reduced_genome_blastn_numthreadstest.slurm`
   - (with num_threads changed and is much faster), 
- `fungi_refseq_reduced_genome_blastn.slurm`
- `archaea_refseq_reduced_genome_blastn.slurm`




##### 03.4.12. Simplify your life by cat'ing all the results together (and deleting empty files)
- change the cat line to whatever blast ran most recently...

Blastouts: 
- `blastn_output_simple_combo.txt` (this is for bacteria)
- `archaea_blastn_output_simple_combo.txt`
- `fungi_blastn_output_simple_combo.txt`

```
# remove lines 1-5 and last line of all the blastn output files, then cat these together and remove the files`
`for file in NODE_*txt; 
do 
sed '1,5d;$d' ${file} > ${file}simple; 
done 

cat *simple > archaea_blastn_output_simple_combo.txt

rm *simple
```

---

#### 03.5. QC on the BLASTn results
- Determine which contigs should go and which should stay in our reference. 
- Also used a guided blastn with our minigenome (blast with minigenome x assembly) to find contigs that contained orthologous genes and should stay in the reference.
- QC flow is found in the .xlsx files (`blastn_output_simple_combo.xlsx` (bacteria), `fungi_blastn_output_simple_combo.xlsx`, `archaea_blastn_output_simple_combo.xlsx`, and an overview of just the list of contigs here: `contigs_lists.xlsx`

More info about the final 2400 contigs kept below, found in list form here `contigs_to_keep.txt`: 

```
260 contigs: matches btwn the assembly x minigenome that do not match with any microbes (keep them all for the reference)
183 contigs: matches btwn the assembly x minigenome that DO match with some microbes (I went through these and qc'd them, we should keep them all for the reference)
1236 contigs: assembly contigs that did not have a blast hit to minigenome or to microbes
721 contigs: assemblyxmicrobe matches (but did not match with minigenomexassembly), these are the ones that passed our filters
```

---

#### 03.6. Extract contigs to keep from the initial assembly
- On Linux machine at school

```
conda install -c bioconda seqtk

seqtk subseq scaffolds_reduced.fasta contigs_to_keep.txt > scaffolds_reduced_contigs_kept.fasta

grep ">" scaffolds_reduced_contigs_kept.fasta | wc -l
2400
```

---

#### 03.7. Index this new 2400 contig reference

##### Slurm: 
- `02_IndexRef_reducedrefnomicrobes.slurm`

---

#### 03.8. Submit to NCBI
Do what they say to do on their README files! 




---
---

### 04. Mitochondrial genome BLASTn

#### Slurms (run in order):
- `Index_PROW981_NCBI.slurm` # first need to index the .fsa file
- `BLASTn_mtgenome_PROW981R1_NCBI`

#### Output:
- `PROW981R1_NCBI_blastout.txt`

---
---


### 05. QUAST

#### Slurm:
- `QUAST_PROW981R1_NCBI.slurm`

#### Output:
- `report.txt`

---
---

### 06. BUSCO

#### Slurms: 
- `BUSCO_PROW981R1_NCBI_arachnida.slurm`
- `BUSCO_PROW981R1_NCBI_arthropoda.slurm`

#### Outputs:
- `short_summary.specific.arachnida_odb10.scaffolds_reduced_contigs_kept_NCBI-200000000.txt`
- `short_summary.specific.arthropoda_odb10.scaffolds_reduced_contigs_kept_NCBI-200000000.txt`


