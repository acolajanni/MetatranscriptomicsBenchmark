#!/bin/bash
# need taxonkit
# conda activate kmcp

module load r/4.2.3


MAIN_PATH="${1:-/shared/projects/microbiome_translocation/data/Douek_Cleveland/}"
PATH_RES="${2:-/shared/projects/microbiome_translocation/results/Douek_Cleveland/krakenuniq/microbialDB/}"
ID="${3:-/shared/projects/microbiome_translocation/data/Douek_Cleveland/sra_list_RNA.txt}"





cd $PATH_RES

filename=${SRA_ID}.output.txt


ID_list=$(awk '{ print $1 }' ${ID})
for SRA_ID in $ID_list ; do
    echo $SRA_ID

    OUTPUT=${SRA_ID}.output.txt

    current_PATH=${PATH_RES}${SRA_ID}
    cd $current_PATH

    # print read + taxid
    awk '{print $2 "\t" $3}' $OUTPUT | taxonkit reformat -I 2 -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" > ReadsClassif.txt

done 