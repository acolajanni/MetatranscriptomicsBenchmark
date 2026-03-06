#!/bin/bash

################################################################################
# > December 2025
# > Script : Cleanifier.sh
# > Function filter reads with cleanifier -to be used after STAR
# @ COLAJANNI Antonin
################################################################################

PATH_RES="${1:-/Douek_cell2021/}"
SEQ_TYPE="${2:-PE}"
ID="${3:-/shared/projects/microbiome_translocation/data/Douek_cell2021/sra_subset_in_manuscript.txt}"
GENOME_input="${4:-hg19}"
PATH_INPUT="${5:-/shared/projects/microbiome_translocation/data/Douek_cell2021/STAR_mapping_hg19/}"
Index="${6:exact}"


if [ $Index = "exact" ] ; then 
    index=/shared/projects/microbiome_translocation/doc/cleanifier_index/test/t2t_pangenome_hla_variants_cdna_cleanifier_exact_k29_w33 
else 
    index=/shared/projects/microbiome_translocation/doc/cleanifier_index/test/t2t_pangenome_hla_variants_cdna_cleanifier_index_k29_w33
fi


FileNamePrefix_in=_STAR_unmapped


# Align reads to genome : STAR
FileNamePrefix_out=_cleaned

# Loop through reads
FOLDERS=${PATH_INPUT}${ID}*/

mkdir -p ${PATH_RES}cleaned/
cd ${PATH_RES}cleaned/

ID_list=$(awk '{ print $1 }' ${ID})
echo $ID_list
for SRA_ID in $ID_list ; do
    if [ $SRA_ID = "aaa" ] ; then 
        continue
    fi
    
    mkdir -p ${PATH_RES}cleaned/${SRA_ID}
    output_path=${PATH_RES}cleaned/${SRA_ID}/
      
    input_read=${PATH_INPUT}${SRA_ID}/${SRA_ID}${FileNamePrefix_in}.fastq
    
    ### Map against reference 
    # cleanifier filter --index $index --fastq $input_read \
    #   --prefix ${output_path}${SRA_ID}${FileNamePrefix_out} \
    #   --sensitive --keep-host --threads 8 
      
    # ### 
    # awk 'BEGIN{
    #     tax="9606\tEukaryota\tChordata\tMammalia\tPrimates\tHominidae\tHomo\tHomo sapiens"
    #   }
    #   (/^@SRR/){
    #     gsub(/^@/,"",$1);      # remove leading @
    #     print $1 "\t" tax
    #   }' ${output_path}${SRA_ID}${FileNamePrefix_out}_filter.fastq > ${output_path}ReadsClassif_cleanifier.tsv

    
    ### Counting reads in each file
    for f in ${output_path}${SRA_ID}${FileNamePrefix_out}*.fastq ; do
      echo -e "$f\t$(grep -c '^@' "$f")"
    done > ${output_path}Counting_reads.tsv


done
