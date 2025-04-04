#!/bin/bash
#############################
#SBATCH -J retrieve_reads
#SBATCH --array=0-39
#############################

echo "#############################" 
echo "User:" $USER
echo "Date:" `date`
echo "Host:" `hostname`
echo "Directory:" `pwd`
echo "SLURM_JOBID:" $SLURM_JOBID
echo "SLURM_SUBMIT_DIR:" $SLURM_SUBMIT_DIR
echo "SLURM_JOB_NODELIST:" $SLURM_JOB_NODELIST
echo "#############################" 

################################################################################
# > february 2024
# > Script : Extract kraken reads
# > Function : based on kraken classification, extract reads like human, unclassified, etc
# @ COLAJANNI Antonin
################################################################################

PATH_unmapped="${1:-/shared/projects/microbiome_translocation/results/Douek_Cleveland/STAR_mapping_hg19/}"
PATH_RES="${2:-/shared/projects/microbiome_translocation/results/Douek_Cleveland/}"
PATH_output="${3:-/shared/projects/microbiome_translocation/results/Douek_Cleveland/hybrid/}"
wanted_classif="${4-unclassified}"
ID="${5:-/shared/projects/microbiome_translocation/data/Douek_Cleveland/sra_subset_in_manuscript.txt}"
compressed="${6:-FALSE}"

module load seqkit

declare -a SRA_IDs  
for folder in $(cat $ID) ;do
    SRA_IDs=("${SRA_IDs[@]}" "$folder")
done

SRA_ID=${SRA_IDs[$SLURM_ARRAY_TASK_ID]}
echo $SRA_ID

if [ "$SRA_ID" = "aaa" ] ; then
    exit
fi

extension=.fastq

if [ $compressed = "TRUE" ] ; then 
    extension=.fastq.gz
fi


SRA_DIR=${PATH_RES}/${SRA_ID}/



### 1: get wanted read IDs
ReadsClassif=${SRA_DIR}ReadsClassif.txt

if [ $wanted_classif = "human" ] ; then 
    awk '!/Mammalia/ {print $1}' $ReadsClassif > ${SRA_DIR}/human_ID.txt
elif [ $wanted_classif = "unclassified" ] ; then
    awk '$2 == 0 || $2 == 1 || $2 == 131567 || $2 == 155900 || $2 == 198431 {print $1}' "$ReadsClassif" > "${SRA_DIR}/unclassified_ID.txt"
fi

### 2: Fetch them in fastq files
mkdir -p ${PATH_output}/${SRA_ID}/ 
outdir=${PATH_output}/${SRA_ID}/ 


#Retrieve fasta sequence associated with these contig ID -- add /1 and /2 to read IDs 
seqkit grep --pattern-file <(sed 's/$/\/1/' ${SRA_DIR}/${wanted_classif}_ID.txt) \
    ${PATH_unmapped}/${SRA_ID}/${SRA_ID}_*_1${extension} > \
    ${outdir}/${SRA_ID}_unmapped_1.fastq

seqkit grep --pattern-file <(sed 's/$/\/2/' ${SRA_DIR}/${wanted_classif}_ID.txt) \
    ${PATH_unmapped}/${SRA_ID}/${SRA_ID}_*_2${extension} > \
    ${outdir}/${SRA_ID}_unmapped_2.fastq