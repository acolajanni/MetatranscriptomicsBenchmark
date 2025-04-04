#!/bin/bash

#############################
#SBATCH -J Blast
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
module load r/4.2.3
module load blast/2.14.0
module load python/3.9 
module load seqkit/2.1.0 
module load use.own


path_to_contigs="${1:-/shared/projects/microbiome_translocation/results/Douek_Cleveland/Contigs/}"
ID="${2:-/shared/projects/microbiome_translocation/data/Douek_Cleveland/sra_list_RNA.txt}"
use_unmapped_contigs="${3:-FALSE}"
restart_Blast="${4:-FALSE}"

cd $path_to_contigs

# Storing folder names in an array 
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

echo ${SRA_IDs[@]:$index_iter:$n_operation}

# Using the SLURM array task ID to slice the array of SRA ids
for srr in ${SRA_IDs[@]:$index_iter:$n_operation} ; do

    if [ "$srr" = "aaa" ] ; then
        exit
    fi

    mkdir -p ${path_to_contigs}${srr}'/'
    cd ${path_to_contigs}${srr}'/'
    
    #remove the "/"
    sra_id=${srr%/}
    echo $sra_id

    #--------------

    if [ "$use_unmapped_contigs" = "FALSE" ] ; then
        contig=${path_to_contigs}${sra_id}/"trinity_out_dir.Trinity.fasta"
    else 
        contig=${path_to_contigs}${sra_id}/Bowtie2_mapping_hg38/trinity_unmapped.Trinity.fasta
    fi

    out=${sra_id}_Blast_query_VRC.txt

    #--------------

    if [ "$restart_Blast" = "TRUE" ] ; then

        seqkit grep --pattern-file ${path_to_contigs}${sra_id}/diff_list.txt \
            ${path_to_contigs}${sra_id}/"trinity_out_dir.Trinity.fasta" > ${path_to_contigs}${sra_id}/restarted.Trinity.fasta

        contig=${path_to_contigs}${sra_id}/restarted.Trinity.fasta
        out=${sra_id}_Blast_query_VRC_p2.txt
    fi



    echo $contig

    #--------------
    # ## BLAST RUN ##
    # VRC Parameters + reduction of max hsps and target_seq # ==> Ne fais pas ce que c'est sensé faire
    # avant : srun blastn ...
    blastn -db /shared/bank/nt/nt_2024-04-12/blast/nt \
       -query $contig \
       -outfmt "7 qseqid sseqid evalue pident qcovhsp sacc length stitle staxids salltitles" \
       -out $out \
       -evalue 0.000000000001 \
       -max_target_seqs 20 \
       -max_hsps 20 \
       -perc_identity 60 \
       -qcov_hsp_perc 60 \
       -num_threads 12
    
    echo "Blast: Done"

    #Splitting between comment line (= log file) and the rest: Prevent error when subject sequence has a '#' in title
    sed -n '/^#/p' ${sra_id}_Blast_query_VRC.txt > ${sra_id}_Blast_query_VRC_comment.txt  
    sed -n '/^#/!p' ${sra_id}_Blast_query_VRC.txt > ${sra_id}_Blast_query_VRC_no_comment.txt


    if [ $restart_Blast = TRUE ] ; then
        cat $out >> ${sra_id}_Blast_query_VRC.txt
    fi  
    
    ##########################
    # Changement: 19/04/24: Taxonomy is no longer done with python package but with taxonkit (= conda activate kmcp)
    ##########################

    # Format query
    python /shared/projects/microbiome_translocation/fastq_scripts/Blast_utils/Get_classification_from_query.py \
       --dir ${path_to_contigs}${sra_id}/ \
       --FileName ${sra_id}_Blast_query_VRC_no_comment.txt \
       --QueryType blast \
       --OutName Blast_query.formated 
    
    cd $path_to_contigs


    # Remove colnames and classify based on taxid (6th column) 
    # k: superkingdom, K:kingdom, etc.

    tail -n +2 ${path_to_contigs}${sra_id}/${sra_id}_Blast_query.formated.tsv \
        | taxonkit reformat -I 6 -F -f "{k}|{K}|{p}|{c}|{o}|{f}|{g}|{s}|{t}"  > ${path_to_contigs}${sra_id}/${sra_id}_Blast_query_classif.formated.tsv 


    echo "Python: Done"

    mkdir -p ${path_to_contigs}${sra_id}/Quantification/

    ## Assess a unique taxon to each contig ##
    Rscript --no-save --no-restore \
        /shared/projects/microbiome_translocation/fastq_scripts/Blast_utils/Filter_reads.r ${sra_id} ${path_to_contigs} \
        /shared/projects/microbiome_translocation/outputFile_${sra_id}.Rout 2>&1



    rm ${path_to_contigs}${sra_id}/${sra_id}_Blast_query_VRC_no_comment.txt

    echo $sra_id >> ${path_to_contigs}/job_completed.txt


    echo "-----"
done



echo "Date end:" `date`
