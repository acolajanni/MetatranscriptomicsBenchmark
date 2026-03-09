#!/bin/bash

#############################
#SBATCH -J simu
#SBATCH --array=0-49
#############################

module purge
module load seqkit

PATH_MAIN=/shared/projects/microbiome_translocation/
PATH_DATABASE=${PATH_MAIN}database_clean/
SEQUENCE_DATABASE=${PATH_DATABASE}sequences/

PATH_SELECTED_GENOMES=${PATH_DATABASE}/selected_genomes/realist/${SLURM_ARRAY_TASK_ID}/

file=selected_genomes_realist.tsv

genomeID_array=$(awk -v FS='\t' '{print $24}' ${PATH_SELECTED_GENOMES}${file})
taxid_array=$(awk -v FS='\t' '{print $1}' ${PATH_SELECTED_GENOMES}${file})
phylum_array=$(awk -v FS='\t' '{print $17}' ${PATH_SELECTED_GENOMES}${file})
family_array=$(awk -v FS='\t' '{print $20}' ${PATH_SELECTED_GENOMES}${file})
genus_array=$(awk -v FS='\t' '{print $21}' ${PATH_SELECTED_GENOMES}${file})

   
# change into a bash array
genomeID_array=(`echo ${genomeID_array}`)
taxid_array=(`echo ${taxid_array}`)
phylum_array=(`echo ${phylum_array}`)
family_array=(`echo ${family_array}`)
genus_array=(`echo ${genus_array}`)


mkdir -p ${PATH_DATABASE}transcript_database/realist/${SLURM_ARRAY_TASK_ID}/
rm -r ${PATH_DATABASE}transcript_database/realist/${SLURM_ARRAY_TASK_ID}/*

PATH_TRANSCRIPT_DATABASE=${PATH_DATABASE}transcript_database/realist/${SLURM_ARRAY_TASK_ID}/

#${!array[@]} ==> returns the indices
for i in ${!genomeID_array[@]} ; do 

    genomeID=${genomeID_array[$i]}
    taxid=${taxid_array[$i]}
    phylum=${phylum_array[$i]}
    family=${family_array[$i]}
    genus=${genus_array[$i]}


    if [ ! -f "${PATH_TRANSCRIPT_DATABASE}/${genus}_transcripts.fasta" ]; then
        touch "${PATH_TRANSCRIPT_DATABASE}/${genus}_transcripts.fasta"
    fi


    ### CAREFULL: Phyla must not have blanc space in names ==> Unclassified Viruses Phylum 

    #Taking the fna.gz files and append them in the phylum transcripts file
    #Add taxid to the header of sequences in the fasta file +
    #Shortening FASTA headers  
    seqkit replace -p '(.+)' -r taxid:${taxid}'_${1}' ${SEQUENCE_DATABASE}${genomeID}/*.fna.gz \
        | sed '/^>/s/^>\([^ ]*\) .*/>\1 /' \
        | sed  's/lcl|//g' \
        | sed 's/ *$//g' >> ${PATH_TRANSCRIPT_DATABASE}/${genus}_transcripts.fasta


done


# seqkit replace --threads 10 -p '(.+)' -r taxid:9606'_${1}' /shared/projects/microbiome_translocation/database/Grch38p14/GCF_000001405.40_GRCh38.p14_rna_from_genomic.fna.gz \
#     | sed '/^>/s/^>\([^ ]*\) .*/>\1 /' \
#     | sed  's/lcl|//g' \
#     | sed 's/ *$//g' > ${PATH_TRANSCRIPT_DATABASE}Homo_sapiens_rna_transcripts.fasta




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
