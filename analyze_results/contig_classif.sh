#!/bin/bash

#############################
#SBATCH -J Blast
#SBATCH --array=0-54
#############################

path_to_contigs="${1:-/shared/projects/microbiome_translocation/results/Douek_Cleveland/Contigs/}"
ID="${2:-/shared/projects/microbiome_translocation/data/Douek_Cleveland/sra_list_RNA.txt}"


declare -a SRA_IDs  
for folder in $(cat $ID) ;do
    SRA_IDs=("${SRA_IDs[@]}" "$folder")
done


# Number of job in the job Array
length=${#SRA_IDs[@]}
# length divided by n, where n is the number of task asked (= number of job)
n=55
#n_operation=$(( (length / n) + 1 ))
n_operation=$(( (length / n) )) 
index_iter=$((n_operation * SLURM_ARRAY_TASK_ID))

echo ${SRA_IDs[@]:$index_iter:$n_operation}




for SRA_ID in ${SRA_IDs[@]:$index_iter:$n_operation} ; do              

    current_path=${path_to_contigs}${SRA_ID}/
    cd $current_path

    echo $SRA_ID

    awk -F"," 'NR>1 {gsub(/"/, "", $1); gsub(/"/, "", $7); gsub(/"/, "", $9); print $1, "\t", $7, "\t", $9}' Contig_classification.csv > ContigsToReads/contigs_classif.txt
    awk -F"\t" '!seen[$2]++ {print $2,"\t","9606","\t","Homo_sapiens"}' ContigsToReads/human_readsToContigs.txt  >> ContigsToReads/contigs_classif.txt
    awk -F"\t" '!seen[$2]++ {print $2,"\t","0","\t","unclassified"}' ContigsToReads/unclassified_readsToContigs.txt >> ContigsToReads/contigs_classif.txt


done
