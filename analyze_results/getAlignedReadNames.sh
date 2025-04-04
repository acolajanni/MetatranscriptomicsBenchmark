#!/bin/bash
#############################
#SBATCH -J get_readNames
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


PATH_DATA="${1:-/shared/projects/microbiome_translocation/data/Douek_Cleveland/}"
prefix="${2:-SRR}"
compressed="${3:-FALSE}"
ID="${4:-/shared/projects/microbiome_translocation/data/Douek_Cleveland/sra_list_RNA.txt}"
data_type="${5:-Simulated}"

script_dir=/shared/projects/microbiome_translocation/fastq_scripts/analyze_results/

# Storing folder names in an array 
declare -a SRA_IDs  
for folder in $(cat $ID) ;do
    SRA_IDs=("${SRA_IDs[@]}" "$folder")
done


SRA_ID=${SRA_IDs[$SLURM_ARRAY_TASK_ID]}


# Step 1- Retrieve read names
if [[ $data_type == "Simulated" ]]; then

    folders=(Bowtie2_mapping_hg38 STAR_mapping_hg38 Bowtie2_mapping_chm13 STAR_mapping_chm13 Bowtie2_mapping_hg19 STAR_mapping_hg19)
    previous_folders=(raw_reads_200 Bowtie2_mapping_hg38 STAR_mapping_hg38 Bowtie2_mapping_chm13 STAR_mapping_chm13 Bowtie2_mapping_hg19)

    bash ${script_dir}getReadNames.sh ${PATH_DATA}/raw_reads_200/ \
        Simulated TRUE FALSE $SRA_ID


else
   folders=(Trimmed_reads Bowtie2_mapping_hg38 STAR_mapping_hg38 Bowtie2_mapping_chm13 STAR_mapping_chm13 Bowtie2_mapping_hg19 STAR_mapping_hg19)
   previous_folders=(raw_reads Trimmed_reads Bowtie2_mapping_hg38 STAR_mapping_hg38 Bowtie2_mapping_chm13 STAR_mapping_chm13 Bowtie2_mapping_hg19)
fi






PATH_aligned_ReadsNames=${PATH_DATA}/aligned_ReadsNames/
mkdir -p $PATH_aligned_ReadsNames

# Loop through the array using indexes
for index in "${!folders[@]}"; do


    folder=${folders[index]}
    previous_folder=${previous_folders[index]}

    path_reads=${PATH_DATA}${folder}/
    raw_read=FALSE
    compressed=FALSE


    if [[ $folder == "raw_reads" || $folder == "Trimmed_reads" || $folder == "raw_reads_200" ]]; then
        raw_read=TRUE
        compressed=TRUE
    fi 

    if [[ $data_type == "Simulated" ]] ; then
        raw_reads=FALSE
    fi

    # # @1: global read path
    # # @2: prefix of read name in the fastq
    # # @3: Boolean, fastq.gz or fastq
    # # @4: Boolean, raw reads or not
    # # @5: char, SRA id // identifiant of fastq file
    bash ${script_dir}getReadNames.sh \
       $path_reads $prefix $compressed $raw_read $SRA_ID

   
    if [[ $folder != "raw_reads" || $folder != "raw_reads_200" ]] ; then
        mkdir -p ${PATH_aligned_ReadsNames}${folder}


        # Step 2: Find IDs in second file that are NOT in first file
        #grep -v -F -x -f \
        #    ${PATH_DATA}${previous_folder}/ReadNames/${SRA_ID}_ReadNames.txt \
        #    ${PATH_DATA}${folder}/ReadNames/${SRA_ID}_ReadNames.txt > ${PATH_aligned_ReadsNames}/${folder}/${SRA_ID}_missing.txt
        # Find lines in the second file that are NOT in the first file

        echo First file: ${PATH_DATA}${folder}/ReadNames/${SRA_ID}_ReadNames.txt
        echo $(wc -l ${PATH_DATA}${folder}/ReadNames/${SRA_ID}_ReadNames.txt)

        echo Second file: ${PATH_DATA}${previous_folder}/ReadNames/${SRA_ID}_ReadNames.txt
        echo $(wc -l ${PATH_DATA}${previous_folder}/ReadNames/${SRA_ID}_ReadNames.txt)

        
        diff --new-line-format="" --unchanged-line-format="" \
            ${PATH_DATA}${previous_folder}/ReadNames/${SRA_ID}_ReadNames.txt \
            ${PATH_DATA}${folder}/ReadNames/${SRA_ID}_ReadNames.txt > ${PATH_aligned_ReadsNames}/${folder}/${SRA_ID}_missing.txt

    echo '-------------------------'

    fi
done


### Adding the last one
folder=STAR_mapping_hg19
previous_folder=Bowtie2_mapping_hg19
mkdir -p ${PATH_aligned_ReadsNames}${folder}

# Step 2: Find IDs in second file that are NOT in first file
# grep -v -F -x -f \
#     ${PATH_DATA}${previous_folder}/ReadNames/${SRA_ID}_ReadNames.txt \
#     ${PATH_DATA}${folder}/ReadNames/${SRA_ID}_ReadNames.txt > ${PATH_aligned_ReadsNames}/${folder}/${SRA_ID}_missing.txt


diff --new-line-format="" --unchanged-line-format="" \
    ${PATH_DATA}${previous_folder}/ReadNames/${SRA_ID}_ReadNames.txt \
    ${PATH_DATA}${folder}/ReadNames/${SRA_ID}_ReadNames.txt > ${PATH_aligned_ReadsNames}/${folder}/${SRA_ID}_missing.txt



### Synthetizing information: concatenate
mkdir ${PATH_aligned_ReadsNames}/HumanReads/

# Output file
output_file=${PATH_aligned_ReadsNames}/HumanReads/${SRA_ID}_human.txt

# Empty the output file before appending
> "$output_file"

# Loop through files and add source information
for file in ${PATH_aligned_ReadsNames}/*/${SRA_ID}_missing.txt; do
    source_name=$(basename "$(dirname "$file")")  # Get the folder name as the source
    awk -v source="$source_name" '{print $0 "\t" source}' "$file" >> "$output_file"
done