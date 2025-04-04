#! /bin/bash

################################################################################
# > september 2023
# > Script : Bowtie align.sh
# > Function : Using Bowtie2 to align reads
# @ COLAJANNI Antonin
################################################################################

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


PATH_RES="${1:-/shared/projects/microbiome_translocation/data/Douek_cell2021/}"
SEQ_TYPE="${2:-PE}"
ID="${3:-/shared/projects/microbiome_translocation/data/Douek_cell2021/sra_subset_in_manuscript.txt}"
GENOME="${4:-hg38}"
PATH_INPUT="${5:-/shared/projects/microbiome_translocation/data/Douek_cell2021/Trimmed_reads/}"
all_reads_in_one_folder="${6:-TRUE}"
output_if_empty="${7:-FALSE}"
compressed="${8:-TRUE}"

module purge
module load samtools/1.15.1
module load bowtie2/2.5.1
module load bedtools/2.30.0 

PATH_DATA=/shared/projects/microbiome_translocation/data/
#PATH_RES=$(echo $PATH_DATA"Douek_cell2021/")

#############################################################################################################################
# Align reads to genome : BowTie2

FOLDERS=${PATH_RES}${ID}*/

if [ $GENOME = "hg38" ] ; then 
    #Bowtie_index=/shared/bank/homo_sapiens/GRCh38.p14/latest_ensembl/indexes/bowtie2-2.5.1
    Bowtie_index=/shared/projects/microbiome_translocation/database/homo_sapiens/GRCh38.p14/bowtie2-2.5.1/GRCh38.p14
    FileNamePrefix=Bowtie2_hg38

elif [ $GENOME = "hg19" ] ; then 
    Bowtie_index=/shared/bank/homo_sapiens/hg19/bowtie2/hg19
    Bowtie_index=/shared/projects/microbiome_translocation/database/homo_sapiens/GRCh37.p13/bowtie2-2.5.1/GRCh37.p13

    FileNamePrefix=Bowtie2_hg19

elif [ $GENOME = "chm13" ] ; then 
    Bowtie_index=/shared/projects/microbiome_translocation/database/homo_sapiens/T2T-CHM13v2.0/bowtie2-2.5.1/T2T-CHM13v2.0
    FileNamePrefix=Bowtie2_CHM13
fi


if [ $compressed = "TRUE" ] ; then 
    extension=.fastq.gz
else 
    extension=.fastq 
fi



mkdir -p ${PATH_RES}Bowtie2_mapping_${GENOME}/

ID_list=$(awk '{ print $1 }' ${ID})
for SRA_ID in $ID_list ; do
#for path_to_folders in $FOLDERS ; do                   # loop through the path that ends with the SRA id (SRRxxxxxx)
#    path_to_folders="${path_to_folders%/}"             # strip trailing slash (if any)
#    SRA_ID="${path_to_folders##*/}"                    # We extrat the last folder name : the SRA ID
#    SRA_ID="${SRA_ID%\*}" 
    if [ $SRA_ID = "aaa" ] ; then 
        continue
    fi

    echo $SRA_ID
    
    if [ $SEQ_TYPE = "PE" ] ; then 

        mkdir -p ${PATH_RES}Bowtie2_mapping_${GENOME}/${SRA_ID}/
        cd ${PATH_RES}Bowtie2_mapping_${GENOME}/${SRA_ID}/



        echo - Paired-end -


        if [ $all_reads_in_one_folder = "FALSE" ] ; then 
            PATH_reads=${PATH_INPUT}${SRA_ID}/
        else 
            PATH_reads=${PATH_INPUT}
        fi

        # No mixed: always search for both read for a pair 
        # no discordant, don't match reads if both reads of the pair are not aligned
        bowtie2 -p 20 --no-mixed --no-discordant \
            -x $Bowtie_index \
            -1 ${PATH_reads}${SRA_ID}*_1*${extension} \
            -2 ${PATH_reads}${SRA_ID}*_2*${extension} | samtools view -bS - >  ${SRA_ID}_${FileNamePrefix}_mapping.bam
                


        samtools view -b -@ 20 -f 4 ${SRA_ID}_${FileNamePrefix}_mapping.bam \
            | bamToFastq -i - \
                -fq ${SRA_ID}_${FileNamePrefix}_unmapped_1.fastq \
                -fq2 ${SRA_ID}_${FileNamePrefix}_unmapped_2.fastq

        if [ $output_if_empty = "TRUE" ] ; then 

            # Check if the file is empty
            if [ ! -s "${SRA_ID}_${FileNamePrefix}_unmapped_1.fastq" ]; then
                echo "File is empty, copying previous fastq file in here"
                cp ${PATH_reads}${SRA_ID}*_1*.fastq.gz ${PATH_RES}Bowtie2_mapping_${GENOME}/${SRA_ID}/
                cp ${PATH_reads}${SRA_ID}*_2*.fastq.gz ${PATH_RES}Bowtie2_mapping_${GENOME}/${SRA_ID}/

                gunzip ${PATH_RES}Bowtie2_mapping_${GENOME}/${SRA_ID}/*_1*.fastq.gz
                gunzip ${PATH_RES}Bowtie2_mapping_${GENOME}/${SRA_ID}/*_2*.fastq.gz

            fi
        fi

    elif [ $SEQ_TYPE = "trinity" ] ; then 

        # Specially for trinity contigs - remove contigs that map on human genome
        mkdir -p ${PATH_RES}${SRA_ID}/Bowtie2_mapping_${GENOME}/
        cd ${PATH_RES}${SRA_ID}/Bowtie2_mapping_${GENOME}

        PATH_reads=${PATH_INPUT}${SRA_ID}/trinity_out_dir.Trinity.fasta

        bowtie2 -p 20 --no-mixed --no-discordant \
            -x $Bowtie_index \
            -f ${PATH_reads} | samtools view -bS - >  contigs_${FileNamePrefix}_mapping.bam
                
        
        samtools fasta -f 4 contigs_${FileNamePrefix}_mapping.bam > trinity_unmapped.Trinity.fasta


    else 

        echo - Single-end -
        cd $PATH_RES"Bowtie2_mapping/"

        GSM=$(echo $PATH_RES"Trimmed_reads/GSM*")
        for merged_reads in $GSM ; do

            path_to_read="${merged_reads%/}"             # strip trailing slash (if any)
            GSM_ID="${path_to_read##*/}"                 # remove everything before the last '/'
            GSM_ID="${GSM_ID%%.*}"                        # Remove everything after the point(s) (doing it twice)

            echo $GSM_ID
            mkdir -p $(echo $PATH_RES"Bowtie2_mapping/$GSM_ID") 
            cd $PATH_RES"Bowtie2_mapping/"$GSM_ID

            bowtie2 -p 20 -x $Bowtie_index \
                -U ${merged_reads} | samtools view -bS - >  ${GSM_ID}_Bowtie2_mapping.bam

            gzip $merged_reads

            # Retrive BAM parts for unaligned reads 
            samtools view -u  -f 4 ${GSM_ID}_Bowtie2_mapping.bam  > unmapped.bam  # Unmapped reads 


            bamToFastq -i unmapped.bam \
                -fq ${GSM_ID}_${FileNamePrefix}_unmapped.fastq 

        done
        # Break the loop of SRA_id 
        break
    fi

    echo " --- "


done

echo "----------------------------------"
echo " BOWTIE2 - END "


