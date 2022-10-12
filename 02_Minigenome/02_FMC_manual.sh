#### Exons only ####
# PWD: /scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/CERW_717_R_ACAGGCGC-AGGCAGAG_out/CERW_717_R_ACAGGCGC-AGGCAGAG_mt/stitch

file_mt95="CERW_717_R_ACAGGCGC-AGGCAGAG_mt95_ref.txt"

while IFS= read -r line
do {
gene=/scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/CERW_717_R_ACAGGCGC-AGGCAGAG_out/CERW_717_R_ACAGGCGC-AGGCAGAG_mt/stitch/*${line}*stitched_exons.fasta
awk '/^>/{if(N)exit;++N;} {print;}' $gene > "${line}_exons_mt95_ref.fasta" 
} done <"$file_mt95"

awk 'NR == FNR { o[n++] = $0; next } /^>/ && i < n { $0 = ">" o[i++] } 1' CERW_717_R_ACAGGCGC-AGGCAGAG_mt95_ref.txt *_exons_mt95_ref.fasta > CERW_717_R_ACAGGCGC-AGGCAGAG_merged_exons_mt95_ref.fasta



# PWD: /scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/CERW_717_R_ACAGGCGC-AGGCAGAG_out/CERW_717_R_ACAGGCGC-AGGCAGAG_nu/stitch

file_nu95="CERW_717_R_ACAGGCGC-AGGCAGAG_nu95_ref.txt"

while IFS= read -r line
do {
gene=/scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/CERW_717_R_ACAGGCGC-AGGCAGAG_out/CERW_717_R_ACAGGCGC-AGGCAGAG_nu/stitch/*${line}*stitched_exons.fasta
awk '/^>/{if(N)exit;++N;} {print;}' $gene > "${line}_exons_nu95_ref.fasta" 
} done <"$file_nu95"

awk 'NR == FNR { o[n++] = $0; next } /^>/ && i < n { $0 = ">" o[i++] } 1' CERW_717_R_ACAGGCGC-AGGCAGAG_nu95_ref.txt *_exons_nu95_ref.fasta > CERW_717_R_ACAGGCGC-AGGCAGAG_merged_exons_nu95_ref.fasta




# PWD: /scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/PROW_831_R_GATCTATC-AGCCTCAT_out/PROW_831_R_GATCTATC-AGCCTCAT_mt/stitch

file_mt95="PROW_831_R_GATCTATC-AGCCTCAT_mt95_ref.txt"

while IFS= read -r line
do {
gene=/scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/PROW_831_R_GATCTATC-AGCCTCAT_out/PROW_831_R_GATCTATC-AGCCTCAT_mt/stitch/*${line}*stitched_exons.fasta
awk '/^>/{if(N)exit;++N;} {print;}' $gene > "${line}_exons_mt95_ref.fasta" 
} done <"$file_mt95"

awk 'NR == FNR { o[n++] = $0; next } /^>/ && i < n { $0 = ">" o[i++] } 1' PROW_831_R_GATCTATC-AGCCTCAT_mt95_ref.txt *_exons_mt95_ref.fasta > PROW_831_R_GATCTATC-AGCCTCAT_merged_exons_mt95_ref.fasta




# PWD: /scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/PROW_831_R_GATCTATC-AGCCTCAT_out/PROW_831_R_GATCTATC-AGCCTCAT_nu/stitch

file_nu95="PROW_831_R_GATCTATC-AGCCTCAT_nu95_ref.txt"

while IFS= read -r line
do {
gene=/scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/PROW_831_R_GATCTATC-AGCCTCAT_out/PROW_831_R_GATCTATC-AGCCTCAT_nu/stitch/*${line}*stitched_exons.fasta
awk '/^>/{if(N)exit;++N;} {print;}' $gene > "${line}_exons_nu95_ref.fasta" 
} done <"$file_nu95"

awk 'NR == FNR { o[n++] = $0; next } /^>/ && i < n { $0 = ">" o[i++] } 1' PROW_831_R_GATCTATC-AGCCTCAT_nu95_ref.txt *_exons_nu95_ref.fasta > PROW_831_R_GATCTATC-AGCCTCAT_merged_exons_nu95_ref.fasta




#### Best filtered only ####

# PWD: /scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/CERW_717_R_ACAGGCGC-AGGCAGAG_out/CERW_717_R_ACAGGCGC-AGGCAGAG_mt/filtered
### NOTE: do not forget the beginning / before PATH on gene= line :) 

file_mt95="CERW_717_R_ACAGGCGC-AGGCAGAG_mt95_ref.txt"
while IFS= read -r line
do {
gene=/scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/CERW_717_R_ACAGGCGC-AGGCAGAG_out/CERW_717_R_ACAGGCGC-AGGCAGAG_mt/filtered/*${line}*filtered_contigs.fasta
awk '/^>/{if(N)exit;++N;} {print;}' $gene > "${line}_best_mt95_ref.fasta" 
} done <"$file_mt95"
awk 'NR == FNR { o[n++] = $0; next } /^>/ && i < n { $0 = ">" o[i++] } 1' CERW_717_R_ACAGGCGC-AGGCAGAG_mt95_ref.txt *_best_mt95_ref.fasta > CERW_717_R_ACAGGCGC-AGGCAGAG_merged_best_mt95_ref.fasta


# PWD: /scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/CERW_717_R_ACAGGCGC-AGGCAGAG_out/CERW_717_R_ACAGGCGC-AGGCAGAG_nu/filtered
file_nu95="CERW_717_R_ACAGGCGC-AGGCAGAG_nu95_ref.txt"
while IFS= read -r line
do {
gene=/scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/CERW_717_R_ACAGGCGC-AGGCAGAG_out/CERW_717_R_ACAGGCGC-AGGCAGAG_nu/filtered/*${line}*filtered_contigs.fasta
awk '/^>/{if(N)exit;++N;} {print;}' $gene > "${line}_best_nu95_ref.fasta" 
} done <"$file_nu95"

awk 'NR == FNR { o[n++] = $0; next } /^>/ && i < n { $0 = ">" o[i++] } 1' CERW_717_R_ACAGGCGC-AGGCAGAG_nu95_ref.txt *_best_nu95_ref.fasta > CERW_717_R_ACAGGCGC-AGGCAGAG_merged_best_nu95_ref.fasta


# PWD: /scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/PROW_831_R_GATCTATC-AGCCTCAT_out/PROW_831_R_GATCTATC-AGCCTCAT_mt/filtered
file_mt95="PROW_831_R_GATCTATC-AGCCTCAT_mt95_ref.txt"
while IFS= read -r line
do {
gene=/scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/PROW_831_R_GATCTATC-AGCCTCAT_out/PROW_831_R_GATCTATC-AGCCTCAT_mt/filtered/*${line}*filtered_contigs.fasta
awk '/^>/{if(N)exit;++N;} {print;}' $gene > "${line}_best_mt95_ref.fasta" 
} done <"$file_mt95"

awk 'NR == FNR { o[n++] = $0; next } /^>/ && i < n { $0 = ">" o[i++] } 1' PROW_831_R_GATCTATC-AGCCTCAT_mt95_ref.txt *_best_mt95_ref.fasta > PROW_831_R_GATCTATC-AGCCTCAT_merged_best_mt95_ref.fasta


# PWD: /scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/PROW_831_R_GATCTATC-AGCCTCAT_out/PROW_831_R_GATCTATC-AGCCTCAT_nu/filtered
file_nu95="PROW_831_R_GATCTATC-AGCCTCAT_nu95_ref.txt"
while IFS= read -r line
do {
gene=/scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/PROW_831_R_GATCTATC-AGCCTCAT_out/PROW_831_R_GATCTATC-AGCCTCAT_nu/filtered/*${line}*filtered_contigs.fasta
awk '/^>/{if(N)exit;++N;} {print;}' $gene > "${line}_best_nu95_ref.fasta" 
} done <"$file_nu95"

awk 'NR == FNR { o[n++] = $0; next } /^>/ && i < n { $0 = ">" o[i++] } 1' PROW_831_R_GATCTATC-AGCCTCAT_nu95_ref.txt *_best_nu95_ref.fasta > PROW_831_R_GATCTATC-AGCCTCAT_merged_best_nu95_ref.fasta




#### merge mt95 exons and nu95 best


cat /scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/CERW_717_R_ACAGGCGC-AGGCAGAG_out/CERW_717_R_ACAGGCGC-AGGCAGAG_mt/stitch/CERW_717_R_ACAGGCGC-AGGCAGAG_merged_exons_mt95_ref.fasta /scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/CERW_717_R_ACAGGCGC-AGGCAGAG_out/CERW_717_R_ACAGGCGC-AGGCAGAG_nu/filtered/CERW_717_R_ACAGGCGC-AGGCAGAG_merged_best_nu95_ref.fasta > /scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/CERW_717_R_ACAGGCGC-AGGCAGAG_out/CERW_717_R_ACAGGCGC-AGGCAGAG_exons_mt95_best_nu95_ref.fasta



cat /scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/PROW_831_R_GATCTATC-AGCCTCAT_out/PROW_831_R_GATCTATC-AGCCTCAT_mt/stitch/PROW_831_R_GATCTATC-AGCCTCAT_merged_exons_mt95_ref.fasta /scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/PROW_831_R_GATCTATC-AGCCTCAT_out/PROW_831_R_GATCTATC-AGCCTCAT_nu/filtered/PROW_831_R_GATCTATC-AGCCTCAT_merged_best_nu95_ref.fasta > /scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/PROW_831_R_GATCTATC-AGCCTCAT_out/PROW_831_R_GATCTATC-AGCCTCAT_exons_mt95_best_nu95_ref.fasta





#### merge mt95 exons and nu95 exons

cat /scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/CERW_717_R_ACAGGCGC-AGGCAGAG_out/CERW_717_R_ACAGGCGC-AGGCAGAG_mt/stitch/CERW_717_R_ACAGGCGC-AGGCAGAG_merged_exons_mt95_ref.fasta /scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/CERW_717_R_ACAGGCGC-AGGCAGAG_out/CERW_717_R_ACAGGCGC-AGGCAGAG_nu/stitch/CERW_717_R_ACAGGCGC-AGGCAGAG_merged_exons_nu95_ref.fasta > /scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/CERW_717_R_ACAGGCGC-AGGCAGAG_out/CERW_717_R_ACAGGCGC-AGGCAGAG_exons_mt95_exons_nu95_ref.fasta



cat /scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/PROW_831_R_GATCTATC-AGCCTCAT_out/PROW_831_R_GATCTATC-AGCCTCAT_mt/stitch/PROW_831_R_GATCTATC-AGCCTCAT_merged_exons_mt95_ref.fasta /scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/PROW_831_R_GATCTATC-AGCCTCAT_out/PROW_831_R_GATCTATC-AGCCTCAT_nu/stitch/PROW_831_R_GATCTATC-AGCCTCAT_merged_exons_nu95_ref.fasta > /scrfs/storage/amatthews/20210412_snp/01_aTRAM_20210507/PROW_831_R_GATCTATC-AGCCTCAT_out/PROW_831_R_GATCTATC-AGCCTCAT_exons_mt95_exons_nu95_ref.fasta


















