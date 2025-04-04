#!/bin/bash

#############################
#SBATCH -J simu
#SBATCH --array=0-9
#############################

module purge
module load seqkit


PATH_MAIN=/shared/projects/microbiome_translocation/
PATH_DATABASE=${PATH_MAIN}database_clean/
SEQUENCE_DATABASE=${PATH_DATABASE}sequences/

PATH_SELECTED_GENOMES=${PATH_DATABASE}${SLURM_ARRAY_TASK_ID}/

genomeID_array=$(awk -v FS='\t' '{print $24}' ${PATH_SELECTED_GENOMES}selected_genomes.tsv)
taxid_array=$(awk -v FS='\t' '{print $1}' ${PATH_SELECTED_GENOMES}selected_genomes.tsv)
phylum_array=$(awk -v FS='\t' '{print $17}' ${PATH_SELECTED_GENOMES}selected_genomes.tsv)

   
# change into a bash array
genomeID_array=(`echo ${genomeID_array}`)
taxid_array=(`echo ${taxid_array}`)
phylum_array=(`echo ${phylum_array}`)

mkdir -p ${PATH_DATABASE}transcript_database/${SLURM_ARRAY_TASK_ID}/
rm -r ${PATH_DATABASE}transcript_database/${SLURM_ARRAY_TASK_ID}/*

PATH_TRANSCRIPT_DATABASE=${PATH_DATABASE}transcript_database/${SLURM_ARRAY_TASK_ID}/

#${!array[@]} ==> returns the indices
for i in ${!genomeID_array[@]} ; do 

    genomeID=${genomeID_array[$i]}
    taxid=${taxid_array[$i]}
    phylum=${phylum_array[$i]}


    if [ ! -f "${PATH_TRANSCRIPT_DATABASE}/${phylum}_transcripts.fasta" ]; then
        touch "${PATH_TRANSCRIPT_DATABASE}/${phylum}_transcripts.fasta"
    fi


    ### CAREFULL: Phyla must not have blanc space in names ==> Unclassified Viruses Phylum 

    #Taking the fna.gz files and append them in the phylum transcripts file
    #Add taxid to the header of sequences in the fasta file +
    #Shortening FASTA headers  
    seqkit replace --threads 12 -p '(.+)' -r taxid:${taxid}'_${1}' ${SEQUENCE_DATABASE}${genomeID}/*.fna.gz \
        | sed '/^>/s/^>\([^ ]*\) .*/>\1 /' \
        | sed  's/lcl|//g' \
        | sed 's/ *$//g' >> ${PATH_TRANSCRIPT_DATABASE}/${phylum}_transcripts.fasta


done



seqkit replace --threads 12 -p '(.+)' -r taxid:9606'_${1}' /shared/projects/microbiome_translocation/database/Grch38p14/GCF_000001405.40_GRCh38.p14_cds_cleaned.fna.gz \
    | sed '/^>/s/^>\([^ ]*\) .*/>\1 /' \
    | sed  's/lcl|//g' \
    | sed 's/ *$//g' > ${PATH_TRANSCRIPT_DATABASE}Homo_sapiens_transcripts.fasta


seqkit replace --threads 12 -p '(.+)' -r taxid:9606'_${1}' /shared/projects/microbiome_translocation/database/Grch38p14/GCF_000001405.40_GRCh38.p14_rna_from_genomic.fna.gz \
    | sed '/^>/s/^>\([^ ]*\) .*/>\1 /' \
    | sed  's/lcl|//g' \
    | sed 's/ *$//g' > ${PATH_TRANSCRIPT_DATABASE}Homo_sapiens_rna_transcripts.fasta


seqkit replace --threads 12 -p '(.+)' -r taxid:9606'_${1}' /shared/projects/microbiome_translocation/database/T2T_CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_rna_from_genomic.fna.gz \
    | sed '/^>/s/^>\([^ ]*\) .*/>\1 /' \
    | sed  's/lcl|//g' \
    | sed 's/ *$//g' \
    | head -n -1 > ${PATH_TRANSCRIPT_DATABASE}Homo_sapiens_t2t_transcripts.fasta




mkdir -p ${PATH_TRANSCRIPT_DATABASE}/dedup_transcripts_ID/
for f in ${PATH_TRANSCRIPT_DATABASE}/*.fasta  ; do 

    filename=$(basename $f)

    echo $filename

    awk '/^>/ {
        header=$0; 
        count[header]++; 
        if (count[header] > 1) {
            print header "_" count[header] - 1;
        } else { 
            print header;
        }
    } !/^>/ {print}' $f > ${PATH_TRANSCRIPT_DATABASE}/dedup_transcripts_ID/${filename}


done


    # # awk '/^>/ {
    # #     header=$0; 
    # #     count[header]++; 
    # #     if (count[header] > 1) {
    # #         print header "_" count[header] - 1;
    # #     } else { 
    # #         print header;
    # #     }
    # # } !/^>/ {print}' ${PATH_TRANSCRIPT_DATABASE}Homo_sapiens_rna_transcripts.fasta > ${PATH_TRANSCRIPT_DATABASE}/dedup_transcripts_ID/Homo_sapiens_rna_transcripts.fasta


    #     awk '/^>/ {
    #     header=$0; 
    #     count[header]++; 
    #     if (count[header] > 1) {
    #         print header "_" count[header] - 1;
    #     } else { 
    #         print header;
    #     }
    # } !/^>/ {print}' ${PATH_TRANSCRIPT_DATABASE}Homo_sapiens_t2t_transcripts.fasta > ${PATH_TRANSCRIPT_DATABASE}/dedup_transcripts_ID/Homo_sapiens_t2t_transcripts.fasta