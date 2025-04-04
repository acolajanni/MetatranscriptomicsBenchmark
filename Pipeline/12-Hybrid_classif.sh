#!/bin/bash

#############################
#SBATCH -J make_hybrid
#SBATCH --array=0-54
#############################


################################################################################
# > March 2025
# > Script : Hybrid_classif
# > Function : Merge classification to create "hybrid" / "Combined" strategy
# @ COLAJANNI Antonin
################################################################################

MAIN_PATH="${1:-/shared/projects/microbiome_translocation/results/Douek_Cleveland/}"
CLASSIF_1="${2:-/shared/projects/microbiome_translocation/results/Douek_Cleveland/Contigs/}"
CLASSIF_2="${3:-/shared/projects/microbiome_translocation/results/Douek_Cleveland/kraken/microbialDB/kuniq/}"
ID="${4:-/shared/projects/microbiome_translocation/data/Douek_Cleveland/sra_list_RNA.txt}"
method_name="${5:-hybrid_Blast-kuniq}"



declare -a SRA_IDs  
for folder in $(cat $ID) ;do
    SRA_IDs=("${SRA_IDs[@]}" "$folder")
done


SRA_ID=${SRA_IDs[$SLURM_ARRAY_TASK_ID]}
echo $SRA_ID



if [[ $CLASSIF_1 == *"/Contigs"* ]]; then
    CLASSIF_1_path=${CLASSIF_1}${SRA_ID}/ContigsToReads/ReadsClassif.txt
else
    CLASSIF_1_path=${CLASSIF_1}${SRA_ID}/ReadsClassif.txt
fi


CLASSIF_2_path=${CLASSIF_2}${SRA_ID}/ReadsClassif.txt



mkdir ${MAIN_PATH}hybrid/
mkdir ${MAIN_PATH}hybrid/${method_name}/
mkdir ${MAIN_PATH}hybrid/${method_name}/${SRA_ID}
outdir=${MAIN_PATH}hybrid/${method_name}/${SRA_ID}/


### List reads to remove:
if [[ $CLASSIF_1 == *"/Contigs"* ]]; then
    
    awk 'NR>1 && $2 != "unclassified" && $2 != "unassembled"' $CLASSIF_1_path | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8}' > ${outdir}ReadsClassif.txt
    cut -f1 ${outdir}ReadsClassif.txt > ${outdir}classified_readnames.txt

else # else = (kraken)

    ### For kraken we have taxids , we can filter with that (or filtering by "superkingdom" containing "unclassified")
    awk -F'\t' -v OFS='\t' '
        $2 ~ /^(0|1|198431|155900|131567|28384)$/ || $3 ~ /unclassified/ { next } 
        { $2 = ""; print $1, substr($0, index($0, $3)) }
        ' $CLASSIF_1_path > ${outdir}ReadsClassif.txt

    cut -f1 ${outdir}ReadsClassif.txt > ${outdir}classified_readnames.txt
    
fi

awk 'NR==FNR {exclude[$1]; next} !($1 in exclude)' ${outdir}classified_readnames.txt $CLASSIF_2_path > ${outdir}remaining_reads.txt

### Second classification always kraken
awk -F'\t' -v OFS='\t' '
{
    if ($2 == "0") {
        for (i=3; i<=NF; i++) $i = "unclassified";
    }
    print $1, substr($0, index($0, $3));
}' ${outdir}remaining_reads.txt > ${outdir}ReadsClassif_2.txt

cut -f1 ${outdir}ReadsClassif_2.txt > ${outdir}classified_readnames_second.txt

### merging files
cat ${outdir}ReadsClassif_2.txt >> ${outdir}ReadsClassif.txt 


rm ${outdir}remaining_reads.txt ${outdir}ReadsClassif_2.txt
