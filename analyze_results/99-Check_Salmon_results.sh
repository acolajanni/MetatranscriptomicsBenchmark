#!/bin/bash
# need taxonkit
# conda activate kmcp

################################################################################
# > september 2023
# > Script : SalmonQuantification
# > Function : Using Salmon for the quantification
# @ COLAJANNI Antonin
################################################################################

PATH_RES="${1:-/shared/projects/microbiome_translocation/results/Douek_cell2021/Contigs/}"
PATH_DATA="${2:-/shared/projects/microbiome_translocation/data/Douek_cell2021/STAR_mapping/}"
SEQ_TYPE="${3:-PE}"
ID="${4:-/shared/projects/microbiome_translocation/data/Douek_cell2021/sra_subset_in_manuscript.txt}"
compressed="${5:-FALSE}"
extension="${6:-.fastq}"

module purge
module load seqkit/2.1.0 
module load salmon/1.10.2


if [ $compressed = "TRUE" ] ; then 
    extension=${extension}.gz
fi



ID_list=$(awk '{ print $1 }' ${ID})
for SRA_ID in $ID_list ; do     

    echo $SRA_ID
        
    cd ${PATH_RES}${SRA_ID}
        
    current_folder=$(pwd)
    QUANT_folder=${current_folder}/Quantification/
        
    mkdir -p $QUANT_folder
    cd $QUANT_folder

    if [ $SEQ_TYPE = "PE" ] ; then 
            #path_to_read_1=${PATH_DATA}${SRA_ID}/${SRA_ID}_STAR_unmapped_1.fastq
            #path_to_read_2=${PATH_DATA}${SRA_ID}/${SRA_ID}_STAR_unmapped_2.fastq

            path_to_read_1=${PATH_DATA}${SRA_ID}/${SRA_ID}*_unmapped_1${extension} 
            path_to_read_2=${PATH_DATA}${SRA_ID}/${SRA_ID}*_unmapped_2${extension} 
    else
            #path_to_read=${PATH_DATA}${SRA_ID}/${SRA_ID}_STAR_unmapped.fastq
            path_to_read=${PATH_DATA}${SRA_ID}/${SRA_ID}*_unmapped${extension} 

    fi


    mkdir -p ${PATH_RES}${SRA_ID}/testing_salmon/

    salmon index -t ${PATH_RES}${SRA_ID}/trinity_out_dir.Trinity.fasta -i ${PATH_RES}${SRA_ID}/testing_salmon/test

    salmon quant -i ${PATH_RES}${SRA_ID}/testing_salmon/test \
        -l A \
        -1 $path_to_read_1 \
        -2 $path_to_read_2 \
        -p 12 \
        --validateMappings \
        --writeMappings=${PATH_RES}${SRA_ID}/testing_salmon/mapping.sam \
        -o ${PATH_RES}${SRA_ID}/testing_salmon/


done