#!/bin/bash

PATH_DATA="${1:-/shared/projects/microbiome_translocation/data/Douek_Cleveland/raw_reads/}"
prefix="${2:-SRR}"
compressed="${3:-FALSE}"
raw_reads="${4:-FALSE}"
SRA_ID="${5:-A0}"

# 1- Retrieve read names

if [[ $raw_reads == "FALSE" ]] ; then
    read_path=${PATH_DATA}${SRA_ID}/${SRA_ID}*_1*.fastq
    echo $read_path
else
    read_path=${PATH_DATA}/${SRA_ID}*_1*.fastq
fi

if [[ $compressed == "TRUE" ]] ; then
    read_path=${read_path}.gz
fi

output_dir=${PATH_DATA}/ReadNames/
mkdir -p $output_dir


if [[ $prefix == "SRR" ]]; then

    if [[ $compressed == "TRUE" ]]; then
        zcat $read_path | grep '^@SRR' | sed 's/^@//; s/\/1$//' | sort > ${output_dir}${SRA_ID}_ReadNames.txt
    else
        grep '^@SRR' $read_path | sed 's/^@//; s/\/1$//' | sort > ${output_dir}${SRA_ID}_ReadNames.txt
    fi

elif [[ $prefix == "Simulated" ]]; then
    # Else it means its simulated. Simulated identifiants starts with R00000001.....

    if [[ $compressed == "TRUE" ]]; then
        zcat $read_path |  grep '^@R0' | sed 's/^@//; s/\/1$//' | sort > ${output_dir}${SRA_ID}_ReadNames.txt
    else
        grep '^@R0' $read_path | sed 's/^@//; s/\/1$//' | sort > ${output_dir}${SRA_ID}_ReadNames.txt
    fi

fi