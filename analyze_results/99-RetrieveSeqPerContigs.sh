#!/bin/bash
#############################
#SBATCH -J SAreads
#SBATCH --array=0-19%10
#############################

module purge
module load seqkit
module load samtools
module load r/4.2.3


MAIN_PATH="${1:-/shared/projects/microbiome_translocation/data/Douek_Cleveland/}"
PATH_RES="${2:-/shared/projects/microbiome_translocation/results/Douek_Cleveland/Contigs/}"
ID="${3:-/shared/projects/microbiome_translocation/data/Douek_Cleveland/sra_list_RNA.txt}"


# cd ${PATH_RES}Bowtie2_mapping_hg38/

# seqkit seq -n trinity_unmapped.Trinity.fasta > Contig_names.txt 


declare -a SRA_IDs  
for folder in $(cat $ID) ;do
    SRA_IDs=("${SRA_IDs[@]}" "$folder")
done

SRA_ID=${SRA_IDs[$SLURM_ARRAY_TASK_ID]}

echo $SRA_ID

# Create file listing for each contigs per kingdom/phylum
Rscript --no-save --no-restore \
    /shared/projects/microbiome_translocation/fastq_scripts/analyze_results/99-SplitContigsPerClassif.R ${SRA_ID} ${PATH_RES} Split \
    /shared/projects/microbiome_translocation/outputFile_${SRA_ID}.Rout 2>&1



PATH_classif=${PATH_RES}${SRA_ID}/contigs_classification/
PATH_quant=${PATH_RES}${SRA_ID}/Quantification/phylum/

cd $$PATH_classif

# Index/sort SAM file produced by Salmon
samtools view -@ 20 -u ${PATH_quant}/phylum_quant_${SRA_ID}.sam > ${PATH_quant}/phylum_quant_${SRA_ID}.bam
samtools sort -@ 20 ${PATH_quant}/phylum_quant_${SRA_ID}.bam > ${PATH_quant}/phylum_quant_${SRA_ID}.sorted
samtools index -@ 20 ${PATH_quant}/phylum_quant_${SRA_ID}.sorted    
BAM_file=${PATH_quant}/phylum_quant_${SRA_ID}.sorted 

touch ${PATH_classif}ReadsToContigs.txt

for sk in ${PATH_classif}*.txt ; do 
        
    # Retrieve sequences from contigs
    filename=$(basename $sk)
    filename=${filename//.txt/}
    echo $filename


    # Separate each fasta file to isolate contigs that belongs to each superkingdom
    seqkit grep --pattern-file ${sk} \
        ${PATH_RES}${SRA_ID}/trinity_out_dir.Trinity.fasta > \
        ${PATH_classif}${filename}.fasta


       # samtools view -@ 20 $BAM_file $(cat $sk) | awk '{print $1 "\t" $3}' > ${PATH_classif}${filename}_readsToContigs.txt

    phylum=${filename//Contigs/}

    #for each phylum do...
    for phy in ${PATH_classif}${phylum}/*.txt ; do 

        # Retrieve sequences from contigs
        filename2=$(basename $phy)
        filename2=${filename2//.txt/}
        echo $filename2
        
        
        seqkit grep --pattern-file ${phy} \
            ${PATH_RES}${SRA_ID}/trinity_out_dir.Trinity.fasta > \
            ${PATH_classif}/${phylum}/${filename2}.fasta

        samtools view -@ 20 $BAM_file $(cat $phy) | awk '{print $1 "\t" $3}' > ${PATH_classif}/${phylum}/${filename2}_readsToContigs.txt
        cat ${PATH_classif}/${phylum}/${filename2}_readsToContigs.txt >> ${PATH_classif}ReadsToContigs.txt



    done

    Rscript --no-save --no-restore \
        /shared/projects/microbiome_translocation/fastq_scripts/analyze_results/99-SplitContigsPerClassif.R ${SRA_ID} ${PATH_RES} Merge \
        /shared/projects/microbiome_translocation/outputFile_${SRA_ID}.Rout 2>&1
    
done

    # samtools view -u test.sam > test.bam
    # samtools sort test.bam > test.sorted
    # samtools index test.sorted     
    # samtools view test.sorted $(cat idtest.txt) | awk '{print $1 "\t" $3}'






