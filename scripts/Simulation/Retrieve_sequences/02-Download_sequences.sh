#!/bin/bash

#############################
#SBATCH -J DL
#SBATCH --array=1-20
#############################



DATABASE_PATH=/shared/projects/microbiome_translocation/database_clean/

cd $DATABASE_PATH
mkdir -p sequences

cd ${DATABASE_PATH}sequences/
current_path=${DATABASE_PATH}sequences/



for ftp_path in $(awk -v FS='\t' '{print $9}' ${DATABASE_PATH}/refseq_selected_genomes_p${SLURM_ARRAY_TASK_ID}.tsv) ; do


    assembly_id="${ftp_path%*/}"
    assembly_id="${assembly_id##*/}"

    mkdir -p ${current_path}${assembly_id}
    cd ${current_path}${assembly_id}
    

    if ! find . -type f -name "*.fna.gz" -size +0c | grep -q '.'; then
        
        echo $ftp_path
        echo "No non-empty files matching *.fna.gz found"
        # If it exists, take rna file, else get the cds file
        wget -r -nd --no-parent -A '*_rna_from_genomic.fna.gz' $ftp_path || wget -r -nd -A '*_cds_from_genomic.fna.gz' $ftp_path
    else
        echo "There is a non-empty file matching *.fna.gz"
        continue
    fi






done

