#!/bin/bash

#############################
#SBATCH -J simu
#SBATCH --array=0-9
#############################

module purge
module load seqkit
module load r/4.2.3


condition_list="${1:-~/data/Simulation/with_replacement/condition_list_Simu2.tsv}"
PATH_DATA="${2:-~/data/Simulation/with_replacement/}"
suffix_reads="${3:-raw}"
restrict_group="${4:-None}"


eukaryota=(
    "Apicomplexa" "Ascomycota" "Basidiomycota" "Chytridiomycota" "Discosea"
    "Euglenozoa" "Evosea" "Fornicata" "Heterolobosea" "Microsporidia" 
    "Mucoromycota" "Oomycota" "Parabasalia")

bacteria=(
  "Actinomycetota"  "Bacillota"  "Bacteroidota"  "Campylobacterota"  "Fusobacteriota"
  "Lentisphaerota"  "Mycoplasmatota" "Pseudomonadota"  "Spirochaetota" "Synergistota"
  "Thermodesulfobacteriota" "Verrucomicrobiota")

viruses=(
  "Artverviricota"  "Cossaviricota"  "Cressdnaviricota"  "Hofneiviricota"  "Peploviricota"
  "Phixviricota"  "Pisuviricota" "Preplasmiviricota"  "Uroviricota"  "unclassified-Viruses-phylum")

archaea=("Euryarchaeota"  "Candidatus-Thermoplasmatota")


# Function to test if a value is in an array
contains() {
  local match="$1"
  shift
  local array=("$@")
  for item in "${array[@]}"; do
    if [[ "$item" == "$match" ]]; then
      return 0  # found
    fi
  done
  return 1  # not found
}


PATH_MAIN=~/

# Retrieve conditions through config file "condition_list.tsv" and other variables
index_row=$(( SLURM_ARRAY_TASK_ID + 1))
echo $index_row

condition=$(tail -n +2 $condition_list | awk -v OFS="\t" "NR==$index_row")

reads_per_phylum=$(echo $condition | cut -d ' ' -f 1 )
read_per_transcript=$(echo $condition | cut -d ' ' -f 2 )
human_read=$(echo $condition | cut -d ' ' -f 3 )
SRA_ID=$(echo $condition | cut -d ' ' -f 4 )
serie_ID=$(echo $condition | cut -d ' ' -f 6 )
human_read_rna=$(echo $condition | cut -d ' ' -f 7 )
human_read_t2t=$(echo $condition | cut -d ' ' -f 8 )
human_dir_OR=$(echo $condition | cut -d ' ' -f 9 )
additional_seq=$(echo $condition | cut -d ' ' -f 10 )


n_transcript=$(( reads_per_phylum / read_per_transcript))


compression=.gz
extension=.fastq


echo $SRA_ID


SIMULATION_DIR=${PATH_DATA}raw_reads_${n_transcript}/
mkdir -p ${PATH_DATA}raw_reads_${n_transcript}/${SRA_ID}
OUTPUT_DIR=${PATH_DATA}raw_reads_${n_transcript}/${SRA_ID}/

R1_out=${OUTPUT_DIR}${SRA_ID}_${suffix_reads}_1${extension}
R2_out=${OUTPUT_DIR}${SRA_ID}_${suffix_reads}_2${extension}

rm $R1_out
touch $R1_out

rm $R2_out
touch $R2_out

rm $R2_out.gz
rm $R1_out.gz

# For each phylum, gather the previously producted reads except for human if we don't want them
for phylum_path in ${SIMULATION_DIR}${serie_ID}/* ; do
    phylum="${phylum_path##*/}"
    echo $phylum


    if [[ "$phylum" == "Homo_sapiens" && "$human_read" == "FALSE" ]]; then
        echo skip
        continue
    fi

    if [[ "$phylum" == "Homo_sapiens_rna" && "$human_read_rna" == "FALSE" ]]; then
        echo skip
        continue
    fi

    if [[ "$phylum" == "Homo_sapiens_t2t" && "$human_read_t2t" == "FALSE" ]]; then
        echo skip
        continue
    fi


    # Only proceed if restrict_group is active and phylum doesn't belong
    if [[ "$restrict_group" != "None" ]]; then
        if [[ "$restrict_group" == "eukaryota" ]] && ! contains "$phylum" "${eukaryota[@]}"; then
            echo skip
            continue
        elif [[ "$restrict_group" == "bacteria" ]] && ! contains "$phylum" "${bacteria[@]}"; then
            echo skip
            continue
        elif [[ "$restrict_group" == "viruses" ]] && ! contains "$phylum" "${viruses[@]}"; then
            echo skip
            continue
        elif [[ "$restrict_group" == "archaea" ]] && ! contains "$phylum" "${archaea[@]}"; then
            echo skip
            continue
        fi
    fi


    echo no skip

    echo " "

    R1=${phylum_path}/${phylum}_${read_per_transcript}rt_R1${extension}${compression}
    R2=${phylum_path}/${phylum}_${read_per_transcript}rt_R2${extension}${compression}

    echo $R1

    # carefull Seqkit removes the compression

    # ADD /1 and /2 at the end of read header to prevent potential problems downstream
    seqkit replace --threads 4 -p '(.+)' -r '${1}/1' $R1 >> $R1_out
    seqkit replace --threads 4 -p '(.+)' -r '${1}/2' $R2 >> $R2_out     

done


if [[ "$human_dir_OR" != "$n_transcript" ]]; then


    phylum_path=${PATH_DATA}raw_reads_${human_dir_OR}/${serie_ID}/Homo_sapiens_rna/
    phylum=Homo_sapiens_rna

    R1=${phylum_path}/${phylum}_${read_per_transcript}rt_R1${extension}${compression}
    R2=${phylum_path}/${phylum}_${read_per_transcript}rt_R2${extension}${compression}

    seqkit replace --threads 4 -p '(.+)' -r '${1}/1' $R1 >> $R1_out
    seqkit replace --threads 4 -p '(.+)' -r '${1}/2' $R2 >> $R2_out     



fi

if [[ "$additional_seq" != "None" ]]; then


    SIMULATION_DIR=${PATH_DATA}raw_reads_${additional_seq}/
    for additionnal_path in ${SIMULATION_DIR}${serie_ID}/* ; do

        additional="${additionnal_path##*/}"
        echo $additional

        R1=${additionnal_path}/${additional}_${read_per_transcript}rt_R1${extension}${compression}
        R2=${additionnal_path}/${additional}_${read_per_transcript}rt_R2${extension}${compression}

        echo $R1

        # carefull Seqkit removes the compression

        # ADD /1 and /2 at the end of read header to prevent potential problems downstream
        seqkit replace --threads 4 -p '(.+)' -r '${1}/1' $R1 >> $R1_out
        seqkit replace --threads 4 -p '(.+)' -r '${1}/2' $R2 >> $R2_out     

done

fi







gzip $R1_out
gzip $R2_out


# rm ${OUTPUT_DIR}${SRA_ID}_unmapped_1${extension}
# rm ${OUTPUT_DIR}${SRA_ID}_unmapped_2${extension}