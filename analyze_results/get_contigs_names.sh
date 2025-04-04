#!/bin/bash


path_to_contigs="${1:-/shared/projects/microbiome_translocation/results/Douek_Cleveland/Contigs/}"
ID="${2:-/shared/projects/microbiome_translocation/data/Douek_Cleveland/sra_list_RNA.txt}"

cd $path_to_contigs

# Storing folder names in an array 
declare -a SRA_IDs  
for folder in $(cat $ID) ;do
    SRA_IDs=("${SRA_IDs[@]}" "$folder")
done



for srr in ${SRA_IDs[@]}; do

    awk 'sub(/^>/, "")' ${path_to_contigs}${srr}/trinity_out_dir.Trinity.fasta > ${path_to_contigs}${srr}/Contigs_name_length.txt

done