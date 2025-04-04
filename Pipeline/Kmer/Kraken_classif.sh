#!/bin/bash
#############################
#SBATCH -J KU
#SBATCH --array=0-54%10
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
# > Script : Kraken.sh
# > Function : Classification with kraken
# @ COLAJANNI Antonin
################################################################################

PATH_unmapped="${1:-/shared/projects/microbiome_translocation/results/Douek_Cleveland/STAR_mapping_hg19/}"
PATH_RES="${2:-/shared/projects/microbiome_translocation/results/Douek_Cleveland/}"
SEQ_TYPE="${3:-PE}"
ALGO="${4-Kuniq}"
DB="${5-kuniq_microbialDB}"
ID="${6:-/shared/projects/microbiome_translocation/data/Douek_Cleveland/sra_subset_in_manuscript.txt}"
compressed="${7:-FALSE}"
extension="${8:-.fastq}"


module purge
module load krakenuniq/1.0.3
module load kraken2



if [ $DB = "kuniq_microbialDB" ] ; then 
    KRAKENDBDIR=/shared/bank/krakenuniq/microbialdb/2020-08-16/
    mkdir -p ${PATH_RES}kraken/microbialDB/
    CURRENT_DIR=${PATH_RES}kraken/microbialDB/

elif [ $DB = "k2_nt" ] ; then
    KRAKENDBDIR=/shared/bank/nt/kraken2/2024-05-30/
    mkdir -p ${PATH_RES}kraken/nt/
    CURRENT_DIR=${PATH_RES}/kraken/nt/

fi


cd $CURRENT_DIR


echo $ALGO

echo $(pwd)

declare -a SRA_IDs  
for folder in $(cat $ID) ;do
    SRA_IDs=("${SRA_IDs[@]}" "$folder")
done

# Number of job in the job Array
length=${#SRA_IDs[@]}
# length divided by n, where n is the number of task asked (= number of job)
n=54
#n_operation=$(( (length / n) + 1 ))
n_operation=$(( (length / n) )) 
index_iter=$((n_operation * SLURM_ARRAY_TASK_ID))




if [ $compressed = "TRUE" ] ; then 
    extension=${extension}.gz
fi


SRA_ID=${SRA_IDs[$SLURM_ARRAY_TASK_ID]}
echo $SRA_ID

if [ "$SRA_ID" = "aaa" ] ; then
    exit
fi


mkdir -p ${CURRENT_DIR}${ALGO}/${SRA_ID}
SRA_DIR=${CURRENT_DIR}${ALGO}/${SRA_ID}/



if [ $ALGO = "kuniq" ] ; then 
    KRAKENDBDIR=/shared/bank/krakenuniq/microbialdb/2020-08-16/

    krakenuniq -t 12 \
        --db ${KRAKENDBDIR} \
        --preload-size 92G \
        --output ${SRA_DIR}${SRA_ID}.output.txt \
        --classified-out ${SRA_DIR}${SRA_ID}_classified#.fastq \
        --unclassified-out ${SRA_DIR}${SRA_ID}_unclassified#.fastq \
        --report ${SRA_DIR}kuniq_Report_${SRA_ID}.txt \
        --exact \
        --paired ${PATH_unmapped}${SRA_ID}/${SRA_ID}*_unmapped_1${extension} ${PATH_unmapped}${SRA_ID}/${SRA_ID}*_unmapped_2${extension} 
        

    report=${SRA_DIR}kuniq_Report_${SRA_ID}.txt \
    # Getting new taxnomy
    tail -n +4 $report | awk '{print $7}' | taxonkit lineage | taxonkit reformat -I 1 -f "{k}|{p}|{c}|{o}|{f}|{g}|{s}|{t}" -F -t > ${SRA_DIR}${SRA_ID}_taxid_lineage_from_old_taxid.txt


elif [ $ALGO = "kraken2" ] ; then

    kraken2 -t 12 \
        --memory-mapping \
        --db ${KRAKENDBDIR} \
        --output ${SRA_DIR}${SRA_ID}.output.txt \
        --classified-out ${SRA_DIR}${SRA_ID}_classified#.fastq \
        --unclassified-out ${SRA_DIR}${SRA_ID}_unclassified#.fastq \
        --report ${SRA_DIR}k2Report_${SRA_ID}.txt \
        --paired ${PATH_unmapped}${SRA_ID}/${SRA_ID}*_unmapped_1${extension} ${PATH_unmapped}${SRA_ID}/${SRA_ID}*_unmapped_2${extension} 



    report=${SRA_DIR}k2Report_${SRA_ID}.txt     
    # Getting new taxnomy
    awk '{print $7}' $report | taxonkit lineage | taxonkit reformat -I 1 -f "{k}|{p}|{c}|{o}|{f}|{g}|{s}|{t}" -F -t > ${SRA_DIR}${SRA_ID}_taxid_lineage_from_old_taxid.txt


elif [ $ALGO = "k2uniq" ] ; then

    ## Minimizer data ==> kuniq
    kraken2 -t 12 \
        --memory-mapping \
        --db ${KRAKENDBDIR} \
        --output ${SRA_DIR}${SRA_ID}.output.txt \
        --classified-out ${SRA_DIR}${SRA_ID}_classified#.fastq \
        --unclassified-out ${SRA_DIR}${SRA_ID}_unclassified#.fastq \
        --report-minimizer-data \
        --report ${SRA_DIR}k2UniqReport_${SRA_ID}.txt \
        --paired ${PATH_unmapped}${SRA_ID}/${SRA_ID}*_unmapped_1${extension} ${PATH_unmapped}${SRA_ID}/${SRA_ID}*_unmapped_2${extension} 
    
    report=${SRA_DIR}k2UniqReport_${SRA_ID}.txt
    # Getting new taxnomy
    awk '{print $7}' $report | taxonkit lineage | taxonkit reformat -I 1 -f "{k}|{p}|{c}|{o}|{f}|{g}|{s}|{t}" -F -t > ${SRA_DIR}${SRA_ID}_taxid_lineage_from_old_taxid.txt

fi



awk -F'\t' '{print $2 "\t" $3}' ${SRA_DIR}${SRA_ID}.output.txt | taxonkit reformat -I 2 -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" -F > ${SRA_DIR}ReadsClassif.txt