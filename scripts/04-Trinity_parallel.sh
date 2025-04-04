#! /bin/bash

#############################
#SBATCH -J Trinity
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



################################################################################
# > september 2023
# > Script : Contig_creation.sh
# > Function : Build contig with Trinity
# @ COLAJANNI Antonin
################################################################################

module purge
module load trinity/2.13.2
module load python/3.9 
module load use.own


PATH_DATA="${1:-/shared/projects/microbiome_translocation/data/Douek_Cleveland/}"
PATH_RES="${2:-/shared/projects/microbiome_translocation/results/Douek_Cleveland/}"
SEQ_TYPE="${3:-PE}"
ID="${4:-/shared/projects/microbiome_translocation/data/Douek_Cleveland/sra_list_RNA.txt}"
PATH_unmapped="${5:-/shared/projects/microbiome_translocation/data/Douek_Cleveland/Bowtie2_mapping_chm13/}"
compressed="${6:-FALSE}"
extension="${7:-.fastq}"


## By default parameter
MIN_CONTIG_LENGTH=200

# Storing folder names in an array 
declare -a SRA_IDs  
for folder in $(cat $ID) ;do
    SRA_IDs=("${SRA_IDs[@]}" "$folder")
done


SRA_ID=${SRA_IDs[SLURM_ARRAY_TASK_ID]}


if [ $extension = ".fasta" ] ; then 
    seqtype=fa
else
    seqtype=fq
fi

if [ $compressed = "TRUE" ] ; then 
    extension=${extension}.gz
fi

mkdir -p $PATH_RES"Contigs/"
Contig_PATH=$PATH_RES"Contigs/"
cd $Contig_PATH


echo $SRA_ID

if [ $SRA_ID = "aaa" ] ; then
    exit
fi

mkdir -p $SRA_ID
cd $SRA_ID
        
echo "$PWD"
    
if [ $SEQ_TYPE = "PE" ] ; then 
    echo - Trinity Paired-end -
        
    rm -r ${Contig_PATH}${SRA_ID}/trinity_out_dir

    Trinity --seqType $seqtype --full_cleanup \
            --left ${PATH_unmapped}${SRA_ID}/${SRA_ID}*_unmapped_1${extension} \
            --right ${PATH_unmapped}${SRA_ID}/${SRA_ID}*_unmapped_2${extension} \
            --min_kmer_cov 1 \
            --min_contig_length $MIN_CONTIG_LENGTH \
            --CPU 20 --max_memory 120G --no_version_check 1>trinity_out_${MIN_CONTIG_LENGTH}_${SRA_ID}.txt 2>trinity_err_${SRA_ID}${FileNamePrefix}.txt

  

else
    echo - Single-end -

    Trinity --seqType $seqtype --full_cleanup \
        	--min_kmer_cov 1 \
            --min_contig_length $MIN_CONTIG_LENGTH \
            --single ${PATH_unmapped}${SRA_ID}/${SRA_ID}${FileNamePrefix}_unmapped.fastq \
            --CPU 20 --max_memory 120G --no_version_check 1>"trinity_out_"$MIN_CONTIG_LENGTH"_"${SRA_ID}${FileNamePrefix}".txt" 2>"trinity_err_"${SRA_ID}${FileNamePrefix}".txt"
        
fi
    
    
cd $Contig_PATH
echo " --- "
