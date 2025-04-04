#! /bin/sh

################################################################################
# > september 2023
# > Script : get filtered contigs.sh
# > Function : retrieve list of reads and extract them from the reads
# @ COLAJANNI Antonin
################################################################################


PATH_RES="${1:-/shared/projects/microbiome_translocation/data/Douek_cell2021/Contigs}"
ID="${2:-/shared/projects/microbiome_translocation/data/Douek_cell2021/sra_subset_in_manuscript.txt}"


ID_list=$(awk '{ print $1 }' ${ID})
for SRA_ID in $ID_list ; do
    echo $SRA_ID

    # Taking the list of filtered ID, and taking the contigs that corresponds to those IDs
    seqkit grep --pattern-file ${PATH_RES}${SRA_ID}/${SRA_ID}_contigs_filtered_id.txt ${PATH_RES}${SRA_ID}/trinity_out_dir.Trinity.fasta > ${PATH_RES}${SRA_ID}/${SRA_ID}_filtered_contigs.Trinity.fasta

    # Taking the list of Plantae ids 
    if [ $(cat ./plantea_contigs_id.txt | wc -l) != 0 ] ; then 
        seqkit grep --pattern-file ${PATH_RES}${SRA_ID}/plantea_contigs_id.txt ${PATH_RES}${SRA_ID}/trinity_out_dir.Trinity.fasta > ${PATH_RES}${SRA_ID}/${SRA_ID}_plantae_contigs.Trinity.fasta
    fi


done