#!/bin/bash

#############################
#SBATCH -J Align
#SBATCH --array=0-54%10
#############################

# useful informations to print
echo "#############################" 
echo "User:" $USER
echo "Date:" `date`
echo "Host:" `hostname`
echo "Directory:" `pwd`
echo "SLURM_JOBID:" $SLURM_JOBID
echo "SLURM_SUBMIT_DIR:" $SLURM_SUBMIT_DIR
echo "SLURM_JOB_NODELIST:" $SLURM_JOB_NODELIST
echo "#############################" 

#############################
module purge
module load sra-tools/2.11.0
module load trimmomatic/0.39
module load bowtie2/2.5.1
module load star/2.7.5a
module load samtools/1.15.1
module load bedtools/2.30.0 
module load seqkit/2.1.0 


MAIN_PATH="${1:-/shared/projects/microbiome_translocation/data/Douek_Cleveland/}"
ID="${2:-/shared/projects/microbiome_translocation/data/Douek_Cleveland/sra_list_RNA.txt}"
RANDOM_NAME="${3-full}"

SCRIPT_DIR=/shared/projects/microbiome_translocation/fastq_scripts/Alignment/

cd $MAIN_PATH


for folder in $(cat $ID) ;do
    SRA_IDs=("${SRA_IDs[@]}" "$folder")
done



echo ${SRA_IDs[SLURM_ARRAY_TASK_ID]}
SRA=${SRA_IDs[SLURM_ARRAY_TASK_ID]}

printf "%s\n" ${SRA} aaa > ${MAIN_PATH}sra_tmp${SLURM_ARRAY_TASK_ID}_${RANDOM_NAME}.txt

id=${MAIN_PATH}sra_tmp${SLURM_ARRAY_TASK_ID}_${RANDOM_NAME}.txt


#######################################################################################################
########################## - F I L T E R I N G   H U M A N   G E N O M E S - ##########################
#######################################################################################################

cd $SCRIPT_DIR




################################### BOWTIE2 + STAR (hg38) ON TRIMMED READS ###################################

echo "----------------"
echo GRCh38 
echo "----------------"

# @1: Path to project
# @2: Sequencing type (PE/SE)
# @3: ID ==> GSM... / SRR...
# @4: version of human genome to align against
# @5: Path of input reads
# @6: logical: if all the reads are in the same folder or FALSE if reads in ../SRR.../read
# @7: logical: if nothing align returns input reads 
# @8: logical: if input reads are compressed or not
bash ./03-Bowtie_align.sh \
   ${MAIN_PATH} \
   PE $id hg38 \
   ${MAIN_PATH}/Trimmed_reads/ \
   TRUE FALSE TRUE


# # @1: Path to project
# # @2: Sequencing type (PE/SE)
# # @3: ID ==> GSM... / SRR...
# # @4: version of human genome to align against
# # @5: Path of input reads
# # @6: Genomic mat : RNA/DNA
# # @7: index to map with (hg38, chm13, ...)
bash ./03bis-STAR_align.sh \
   ${MAIN_PATH} \
   PE $id hg38 \
   ${MAIN_PATH}Bowtie2_mapping_hg38/ \
   $INPUT_SAMPLE \
   RNA hg38

echo "----------------"
echo CHM13 
echo "----------------"

################################### BOWTIE2 + STAR (CHM13/T2T) ON TRIMMED READS ###################################
bash ./03-Bowtie_align.sh \
   ${MAIN_PATH} \
   PE $id chm13 \
   ${MAIN_PATH}STAR_mapping_hg38/ \
   FALSE FALSE FALSE

bash ./03bis-STAR_align.sh \
   ${MAIN_PATH} \
   PE $id chm13 \
   ${MAIN_PATH}Bowtie2_mapping_chm13/ \
   $INPUT_SAMPLE \
   RNA chm13

   

echo "----------------"
echo GRCh37 
echo "----------------"
################################### BOWTIE2 + STAR (hg19) ON TRIMMED READS ###################################

bash ./03-Bowtie_align.sh \
   ${MAIN_PATH} \
   PE $id hg19 \
   ${MAIN_PATH}STAR_mapping_chm13/ \
   FALSE FALSE FALSE

bash ./03bis-STAR_align.sh \
   ${MAIN_PATH} \
   PE $id hg19 \
   ${MAIN_PATH}Bowtie2_mapping_hg19/ \
   $INPUT_SAMPLE \
   RNA hg19
######################################################################################################

echo "-----"




echo "Date end:" `date`
