#!/bin/bash

#############################
#SBATCH -J SalmonSimu
#SBATCH --array=0-49
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


################################################################################
# > september 2023
# > Script : SalmonQuantification
# > Function : Using Salmon for the quantification
# @ COLAJANNI Antonin
################################################################################

#############################

PATH_RES="${1:-~/results/Contigs/}"
PATH_DATA="${2:-~/data/STAR_mapping/}"
SEQ_TYPE="${3:-PE}"
ID="${4:-~/data/sra_list_RNA.txt}"
compressed="${5:-FALSE}"
extension="${6:-.fastq}"
data_type="${7:-Simulated}"
contigs_type="${8:-trinity}"

module purge
module load seqkit/2.1.0 
module load salmon/1.10.2
module load r/4.2.3


# Storing folder names in an array 
declare -a SRA_IDs  
for folder in $(cat $ID) ;do
    SRA_IDs=("${SRA_IDs[@]}" "$folder")
done


SCRIPT_DIR=~/fastq_scripts/


cd $PATH_RES

TAXONOMY=( strain species genus family order class phylum superkingdom unclassified human )

#TAXONOMY=( strain unclassified human )


if [ $compressed = "TRUE" ] ; then 
    extension=${extension}.gz
fi


# Number of job in the job Array
length=${#SRA_IDs[@]}
# length divided by n, where n is the number of task asked (= number of job)
n=50
#n_operation=$(( (length / n) + 1 ))
n_operation=$(( (length / n) )) 
index_iter=$((n_operation * SLURM_ARRAY_TASK_ID))

echo ${SRA_IDs[@]:$index_iter:$n_operation}
# Using the SLURM array task ID to slice the array of SRA ids
for SRA_ID in ${SRA_IDs[@]:$index_iter:$n_operation} ; do   

    echo $SRA_ID
    
    cd ${PATH_RES}${SRA_ID}
    
    current_folder=$(pwd)
    QUANT_folder=${current_folder}/Quantification/
    
    mkdir -p $QUANT_folder
    cd $QUANT_folder
        
    if [ $SEQ_TYPE = "PE" ] ; then 
        #path_to_read_1=${PATH_DATA}${SRA_ID}/${SRA_ID}_STAR_unmapped_1.fastq
        #path_to_read_2=${PATH_DATA}${SRA_ID}/${SRA_ID}_STAR_unmapped_2.fastq

        path_to_read_1=${PATH_DATA}${SRA_ID}/${SRA_ID}*_unmapped_1${extension} 
        path_to_read_2=${PATH_DATA}${SRA_ID}/${SRA_ID}*_unmapped_2${extension} 
    else
        #path_to_read=${PATH_DATA}${SRA_ID}/${SRA_ID}_STAR_unmapped.fastq
        path_to_read=${PATH_DATA}${SRA_ID}/${SRA_ID}*_unmapped${extension} 

    fi

    ### Handle unclassified and Human contigs

    if [ "$contigs_type" = "trinity" ] ; then

        contig_file=${PATH_RES}${SRA_ID}/Bowtie2_mapping_hg38/trinity_unmapped.Trinity.fasta
        unfiltered_contig_file=${PATH_RES}${SRA_ID}/trinity_out_dir.Trinity.fasta

        # Retrieve contig name
        awk '{print $2}' ${PATH_RES}${SRA_ID}/trinity_out_dir.Trinity.fasta.gene_trans_map > ${PATH_RES}${SRA_ID}/Contig_ID_all.txt

    elif [ "$contigs_type" = "rnaSpades" ] ; then
        
        contig_file=${PATH_RES}${SRA_ID}/Bowtie2_mapping_hg38/rnaSpades_unmapped.rnaSpades.fasta
        unfiltered_contig_file=${PATH_RES}${SRA_ID}/transcripts.rnaSpades.fasta

        # Retrieve contig name
        grep "^>" $unfiltered_contig_file | sed 's/^>//' > ${PATH_RES}${SRA_ID}/Contig_ID_all.txt

    fi

    # Get contig name from fasta (unmapped contig against hg38)
    awk 'sub(/^>/, "")' $contig_file > ${PATH_RES}${SRA_ID}/Bowtie2_mapping_hg38/contig_ID_nonHuman.txt

    # Get difference between all contigs and human to get the mapping to human contigs
    sort ${PATH_RES}${SRA_ID}/Bowtie2_mapping_hg38/contig_ID_nonHuman.txt ${PATH_RES}${SRA_ID}/Contig_ID_all.txt \
            | uniq -u > ${QUANT_folder}/contig_ID_human.txt
    

    # getting unclassified contig id by making the difference between all the nonhuman contigs and the classified ones
    sort ${PATH_RES}${SRA_ID}/Bowtie2_mapping_hg38/contig_ID_nonHuman.txt ${QUANT_folder}/contig_ID_superkingdom.txt \
            | uniq -u > ${QUANT_folder}/contig_ID_unclassified.txt


    for lvl in ${TAXONOMY[@]} ; do  
        mkdir -p $lvl"/"    
        cd $lvl"/"
        echo $lvl

        #Retrieve fasta sequence associated with these contig ID
        seqkit grep --pattern-file ${QUANT_folder}/contig_ID_${lvl}.txt \
            $unfiltered_contig_file > \
            ${QUANT_folder}${lvl}/${lvl}_contig.fasta


        salmon_output_file_name_and_path=${QUANT_folder}${lvl}/${lvl}_quant_${SRA_ID}

        # Create index for it
        # salmon index -p 6 -t ${lvl}_contig.fasta -i ${lvl}_salmon 
        
        if [ $SEQ_TYPE = "PE" ] ; then 

            echo "salmon"
            salmon quant -i ${lvl}_salmon \
                -l A \
                -1 $path_to_read_1 \
                -2 $path_to_read_2 \
                -p 6 \
                --validateMappings \
                --writeMappings=${QUANT_folder}${lvl}/${lvl}_quant_${SRA_ID}.sam \
                -o $salmon_output_file_name_and_path

        else
            salmon quant -i ${lvl}_salmon \
                -l A \
                -r $path_to_read \
                -p 6 \
                --validateMappings \
                -o $salmon_output_file_name_and_path
        fi
        

        cd $QUANT_folder

    done

    # Retrieve sam files - reads ==> Contigs ==> Classification
    bash ${SCRIPT_DIR}/Salmon_utils/get_reads_classif_from_salmon.sh \
        ${PATH_RES} ${SRA_ID}


    # Build file: read_name: classif (D P C O F G S T)
    Rscript --no-save --no-restore ${SCRIPT_DIR}/Salmon_utils/Get_salmon_readsClassif.r \
            ${PATH_RES} ${SRA_ID}

    # MISSING: Reads not in contigs ==> add them at the end of ReadsClassif.txt
    bash ${SCRIPT_DIR}/Salmon_utils/getUncontigedReads.sh ${PATH_DATA} ${PATH_RES} ${SRA_ID} $data_type

done


echo "Date end:" `date`
