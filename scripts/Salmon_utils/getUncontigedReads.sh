#!/bin/bash

### Get un-assembled reads by trinity to add them in the file 'ReadsClassif.txt'

PATH_DATA="${1:-/shared/projects/microbiome_translocation/data/Simulation/with_replacement/raw_reads_45/}"
PATH_RES="${2:-/shared/projects/microbiome_translocation/results/Simulation/with_replacement/Contigs/}"
SRA_ID="${3:-A4}"
prefix="${4:-Simulated}"


# 1- Retrieve read names

if [[ $prefix == "SRR" ]]; then
    motif='^@SRR'
elif [[ $prefix == "Simulated" ]]; then
    # Else it means its simulated. Simulated identifiants starts with R00000001.....
    motif='^@R0'

fi


shopt -s nullglob  # Ensures globbing does not return the pattern itself if no file matches

uncompressed_file=(${PATH_DATA}${SRA_ID}/${SRA_ID}*_unmapped_1.fastq)  # Expands to matching file(s)
compressed_file=(${PATH_DATA}${SRA_ID}/${SRA_ID}*_unmapped_1.fastq.gz)  # Expands to matching file(s)

if [[ -f "${uncompressed_file[0]}" ]]; then
    echo "Uncompressed file"
    grep "$motif" "${uncompressed_file[0]}" | sed 's/^@//; s/\/1$//' > "${PATH_DATA}${SRA_ID}/${SRA_ID}_ReadNames.txt"

elif [[ -f "${compressed_file[0]}" ]]; then
    echo "Compressed file"
    zcat "${compressed_file[0]}" | grep "$motif" | sed 's/^@//; s/\/1$//' > "${PATH_DATA}${SRA_ID}/${SRA_ID}_ReadNames.txt"

else
    echo "nothing"
fi




# 2- Get the missing IDs
PATH_CLASSIF=${PATH_RES}${SRA_ID}/ContigsToReads/

awk '{print $1}' ${PATH_CLASSIF}ReadsClassif.txt | tail -n+2 > ${PATH_CLASSIF}classified_Reads.txt

# Step 2: Find IDs in ${SRA_ID}_ReadNames.txt that are NOT in classif_ids.txt
grep -v -F -x -f ${PATH_CLASSIF}classified_Reads.txt ${PATH_DATA}${SRA_ID}/${SRA_ID}_ReadNames.txt > ${PATH_CLASSIF}missing_ids.txt

# Step 3: Add those IDs at the end of Reads Classif
awk '{print $0 "\t'unassembled'\t'unassembled'\t'unassembled'\t'unassembled'\t'unassembled'\t'unassembled'\t'unassembled'\t'unassembled'"}' ${PATH_CLASSIF}missing_ids.txt >> ${PATH_CLASSIF}ReadsClassif.txt




