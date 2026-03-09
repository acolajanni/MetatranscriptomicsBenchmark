#! /bin/bash

#############################
#SBATCH -J Realspades
#SBATCH --array=0-54
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



################################################################################
# > january 2026
# > Script : Contig_creation.sh
# > Function : Build contig with Spades
# @ COLAJANNI Antonin
################################################################################

module purge
module load spades/4.2.0
module load use.own


PATH_RES="${1:-~/results/}"
SEQ_TYPE="${2:-PE}"
ID="${3:-~/data/sra_list_RNA.txt}"
PATH_unmapped="${4:-~/data/STAR_mapping_hg19/}"
compressed="${5:-FALSE}"
extension="${6:-.fastq}"


# Storing folder names in an array 
declare -a SRA_IDs  
for folder in $(cat $ID) ;do
    SRA_IDs=("${SRA_IDs[@]}" "$folder")
done


SRA_ID=${SRA_IDs[$SLURM_ARRAY_TASK_ID]}


# if [ $extension = ".fasta" ] ; then 
#     seqtype=fa
# else
#     seqtype=fq
# fi

if [ $compressed = "TRUE" ] ; then 
    extension=${extension}.gz
fi

mkdir -p $PATH_RES"Contigs_rnaSpades/"
Contig_PATH=$PATH_RES"Contigs_rnaSpades/"
cd $Contig_PATH

# Using the SLURM array task ID to slice the array of SRA ids
# for SRA_ID in ${SRA_IDs[@]:$index_iter:$n_operation} ; do
echo $SRA_ID

if [ $SRA_ID = "aaa" ] ; then
    exit
fi

mkdir -p $SRA_ID
cd $SRA_ID
        
echo "$PWD"

current_dir=${Contig_PATH}/${SRA_ID}/

spades.py --rna -1 ${PATH_unmapped}${SRA_ID}/${SRA_ID}*_unmapped_1${extension} \
                -2 ${PATH_unmapped}${SRA_ID}/${SRA_ID}*_unmapped_2${extension} \
                -o ${current_dir}rnaSpades \
                -t 10 -m 64 



mv ${current_dir}rnaSpades/transcripts.fasta ${current_dir}/transcripts.rnaSpades.fasta
























