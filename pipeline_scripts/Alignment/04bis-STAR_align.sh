#! /bin/sh

################################################################################
# > september 2023
# > Script : STAR build align.sh
# > Function : Using STAR aligner to build genome index and align reads to it
# @ COLAJANNI Antonin
################################################################################


PATH_RES="${1:-/shared/projects/microbiome_translocation/data/}"
SEQ_TYPE="${2:-PE}"
ID="${3:-/shared/projects/microbiome_translocation/data/sra_list_RNA.txt}"
GENOME_input="${4:-hg19}"
PATH_INPUT="${5:-/shared/projects/microbiome_translocation/data/Bowtie2_mapping/}"
GENOMIC_MAT="${6-RNA}"
Index="${7:hg38}"

PATH_DATA=/shared/projects/microbiome_translocation/data/



if [ $Index = "hg38" ] ; then 
    STAR_index=/shared/projects/microbiome_translocation/database/homo_sapiens/GRCh38.p14/star-2.7.11a
elif [ $Index = "chm13" ] ; then 
    STAR_index=/shared/projects/microbiome_translocation/database/homo_sapiens/T2T-CHM13v2.0/star-2.7.11a
elif [ $Index = "hg19" ] ; then 
    STAR_index=/shared/projects/microbiome_translocation/database/homo_sapiens/GRCh37.p13/star-2.7.11a
fi

echo $STAR_index



if [ $GENOME_input = "hg38" ] ; then 
    FileNamePrefix_in=${SRA_ID}_Bowtie2_hg38_unmapped

elif [ $GENOME_input = "hg19" ] ; then 
    FileNamePrefix_in=${SRA_ID}_Bowtie2_hg19_unmapped

elif [ $GENOME_input = "chm13" ] ; then 
    FileNamePrefix_in=${SRA_ID}_Bowtie2_CHM13_unmapped
    
fi

echo $FileNamePrefix_in

# Align reads to genome : STAR
FileNamePrefix_out=_STAR

# Loop through reads
FOLDERS=${PATH_INPUT}${ID}*/

mkdir -p ${PATH_RES}STAR_mapping_${GENOME_input}/
cd ${PATH_RES}STAR_mapping_${GENOME_input}/



ID_list=$(awk '{ print $1 }' ${ID})
echo $ID_list
for SRA_ID in $ID_list ; do
    if [ $SRA_ID = "aaa" ] ; then 
        continue
    fi
    echo $SRA_ID

    mkdir -p ${PATH_RES}STAR_mapping_${GENOME_input}/${SRA_ID}/
    cd  ${PATH_RES}STAR_mapping_${GENOME_input}/${SRA_ID}/


    if [ $SEQ_TYPE = "PE" ] ; then 
        echo - Paired-end -

        echo ${PATH_INPUT}${SRA_ID}/${SRA_ID}${FileNamePrefix_in}"_1.fastq.gz"


        # Check and decompress the second file
        if [[ -f "${PATH_INPUT}${SRA_ID}/${SRA_ID}${FileNamePrefix_in}_2.fastq.gz" ]]; then
            gunzip "${PATH_INPUT}${SRA_ID}/${SRA_ID}${FileNamePrefix_in}_2.fastq.gz"
        fi

        # Check and decompress the second file
        if [[ -f "${PATH_INPUT}${SRA_ID}/${SRA_ID}${FileNamePrefix_in}_1.fastq.gz" ]]; then
            gunzip "${PATH_INPUT}${SRA_ID}/${SRA_ID}${FileNamePrefix_in}_1.fastq.gz"
        fi

        read_path1=${PATH_INPUT}${SRA_ID}/${SRA_ID}${FileNamePrefix_in}"_1.fastq"
        read_path2=${PATH_INPUT}${SRA_ID}/${SRA_ID}${FileNamePrefix_in}"_2.fastq"  


        if [ $GENOMIC_MAT = "RNA" ] ; then 

            STAR --runThreadN 20 \
                --readFilesIn $read_path1 $read_path2  \
                --genomeDir $STAR_index \
                --outSAMtype BAM SortedByCoordinate \
                --outFileNamePrefix ${SRA_ID}${FileNamePrefix_out}_mapping_ \
                --outReadsUnmapped Fastx \
                --outSAMstrandField intronMotif 1>star_stdout.txt 2>star_stderr.txt # strand derived from the intron motif. Reads with inconsistent and/or non-canonical introns are filtered out.

        else
                # Remove the possibility of splicing since the size of an intron should be less than 1bp and 2bp minimal (= impossible)
                STAR --runThreadN 20 \
                    --readFilesIn $read_path1 $read_path2  \
                    --genomeDir $STAR_index \
                    --alignIntronMax 1 --alignIntronMin 2 \
                    --outSAMtype BAM SortedByCoordinate \
                    --outFileNamePrefix ${SRA_ID}${FileNamePrefix_out}_mapping_ \
                    --outReadsUnmapped Fastx \
                    --outSAMstrandField intronMotif 1>star_stdout.txt 2>star_stderr.txt # strand derived from the intron motif. Reads with inconsistent and/or non-canonical introns are filtered out.
        fi


        # Once mapped, we can recompress fatsq 
        #gzip --force ${PATH_INPUT}${SRA_ID}/${SRA_ID}${FileNamePrefix_in}_1.fastq 
        #gzip --force ${PATH_INPUT}${SRA_ID}/${SRA_ID}${FileNamePrefix_in}_2.fastq

        FileNamePrefix=_STAR_

        # 0:N: 00 ==> No match for read 0 with read '00'
        awk '{gsub(/ 0:N:  [0-1][0-1]/,"/1"); print}' ${SRA_ID}${FileNamePrefix_out}_mapping_Unmapped.out.mate1 > ${SRA_ID}${FileNamePrefix}unmapped_1.fastq
        # 1:N: ==> no match for read 1 with read 00
        awk '{gsub(/ 1:N:  [0-1][0-1]/,"/2"); print}' ${SRA_ID}${FileNamePrefix_out}_mapping_Unmapped.out.mate2 > ${SRA_ID}${FileNamePrefix}unmapped_2.fastq    

        rm ${SRA_ID}${FileNamePrefix_out}_mapping_Unmapped.out.mate1
        rm ${SRA_ID}${FileNamePrefix_out}_mapping_Unmapped.out.mate2
        
        gzip --force -c ${SRA_ID}${FileNamePrefix}unmapped_1.fastq > ${SRA_ID}${FileNamePrefix}unmapped_1.fastq.gz
        gzip --force -c ${SRA_ID}${FileNamePrefix}unmapped_2.fastq > ${SRA_ID}${FileNamePrefix}unmapped_2.fastq.gz


    elif [ $SEQ_TYPE = "SE" ] ; then
        echo - Single-end -

        # Expecting one file e.g.  SRAID${FileNamePrefix_in}.fastq.gz
        input_gz=${PATH_INPUT}${SRA_ID}/${SRA_ID}${FileNamePrefix_in}.fastq.gz
        input_fastq=${PATH_INPUT}${SRA_ID}/${SRA_ID}${FileNamePrefix_in}.fastq

        echo "$input_gz"

        # Decompress if needed
        if [[ -f "$input_gz" ]]; then
            gunzip "$input_gz"
        fi

        read_path=${input_fastq}

        # ============= STAR MAPPING =============
        STAR --runThreadN 12 \
            --readFilesIn $read_path \
            --genomeDir $STAR_index \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${SRA_ID}${FileNamePrefix_out}_mapping_ \
            --outReadsUnmapped Fastx \
            --outSAMstrandField intronMotif \
            1>star_stdout.txt 2>star_stderr.txt



        # ========================================
        # STAR produces:  ${SRA_ID}${FileNamePrefix_out}_mapping_Unmapped.out.mate1
        # for single-end, only mate1 is used

        FileNamePrefix=_STAR_

        ### Fix name tags in unmapped reads
        # awk '{gsub(/ 0:N:  [0-1][0-1]/,""); print}' \
        #     ${SRA_ID}${FileNamePrefix_out}_mapping_Unmapped.out.mate1 \
        #     > ${SRA_ID}${FileNamePrefix}unmapped.fastq

        sed 's/ 0:N: //' ${SRA_ID}${FileNamePrefix_out}_mapping_Unmapped.out.mate1 > ${SRA_ID}${FileNamePrefix}unmapped.fastq

        ### Remove original extra STAR output
        rm ${SRA_ID}${FileNamePrefix_out}_mapping_Unmapped.out.mate1

        # Compress resulting unmapped read file
        # gzip --force -c ${SRA_ID}${FileNamePrefix}unmapped.fastq \
        #     > ${SRA_ID}${FileNamePrefix}unmapped.fastq.gz

    fi







    # else 
    #     echo - Single-end -
    #     # map single-end reads to genome
    #     STAR --runThreadN 20 \
    #         --readFilesIn ${PATH_INPUT}/${SRA_ID}/${SRA_ID}${FileNamePrefix_in}".fastq" \
    #         --genomeDir $STAR_index \
    #         --outSAMtype BAM SortedByCoordinate \
    #         --outFileNamePrefix ${SRA_ID}${FileNamePrefix_out}_mapping_ \
    #         --outReadsUnmapped Fastx \
    #         --outSAMstrandField intronMotif 1>star_stdout.txt 2>star_stderr.txt # strand derived from the intron motif. Reads with inconsistent and/or non-canonical introns are filtered out.

    #     FileNamePrefix=_STAR_

    #     mv ${SRA_ID}${FileNamePrefix_out}_mapping_Unmapped.out.mate1 ${SRA_ID}${FileNamePrefix}unmapped.fastq 
    #     gzip --force -c ${SRA_ID}${FileNamePrefix}unmapped.fastq > ${SRA_ID}${FileNamePrefix}unmapped.fastq.gz

    # fi
    # echo " --- "


done

echo "----------------------------------"
echo " STAR - END "


#############################################################################################################################
# Build Genome Index / Not needed


#mkdir -p $PATH_DATA"HG38_STAR/"
#mkdir -p $PATH_GENECODE"index_files/"
#PATH_GENECODE=$(echo $PATH_DATA"HG38_STAR/")

#cd $PATH_GENECODE

#wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/GRCh38.primary_assembly.genome.fa.gz
#wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v44.primary_assembly.annotation.gtf.gz
#gzip -d *.gz

#GENOME=GRCh38.primary_assembly.genome.fa
#GTF=gencode.v44.primary_assembly.annotation.gtf


## Files available in databank 

# for genomeSAindexNbases :  log2(GenomeLength)/2 - 1) human:3.4Gbases ==> 14.85
#STAR --runThreadN 20 --runMode genomeGenerate  \
#    --genomeSAindexNbases 15 \
#    --genomeDir $PATH_GENECODE \
#    --genomeFastaFiles ${GENOME} \
#    --sjdbOverhang 150 \
#    --sjdbGTFfile ${GTF} \
#    --outFileNamePrefix ${PATH_GENECODE}"index_files/Human_genome"


#############################################################################################################################
#############################################################################################################################
#### map previously unmapped paired-end reads to genome (old)
#STAR --runThreadN 20 \
#    --readFilesIn ${PATH_RES}Bowtie2_mapping/${SRA_ID}/${SRA_ID}${FileNamePrefix_in}"_1.fastq" ${PATH_RES}Bowtie2_mapping/${SRA_ID}/${SRA_ID}${FileNamePrefix_in}"_2.fastq" \
#    --genomeDir $STAR_index \
#    --outSAMtype BAM SortedByCoordinate \
#    --outFileNamePrefix ${SRA_ID}${FileNamePrefix_out}_mapping_ \
#    --outSAMunmapped Within

#############################################################################################################################
# Reads Quality check (FASTQC + MultiQC)

#PATH_MAIN=/shared/projects/microbiome_translocation/
#PATH_RES=$(echo $PATH_MAIN"results/")
#PATH_DATA=$(echo $PATH_MAIN"data/"Mapped_reads)

#cd $PATH_RES
#fastqc ${PATH_DATA}/*/*_Pre_ART_PlasmaAligned.sortedByCoord.out.bam --format bam --outdir ${PATH_RES} 
#multiqc ${PATH_RES}/*/*_Pre_ART_PlasmaAligned.sortedByCoord* -n "Pre_ART_Plasma"





############################################################################################################################# (old)

# Retrive BAM parts for unaligned reads (old)
#samtools view -u  -f 4 -F 264 ${SRA_ID}${FileNamePrefix}Aligned.sortedByCoord.out.bam  > tmps1.bam  # Unmapped reads whose mate is mapped
#samtools view -u -f 8 -F 260 ${SRA_ID}${FileNamePrefix}Aligned.sortedByCoord.out.bam  > tmps2.bam   # mapped reads whose mate is unmapped
#samtools view -u -f 12 -F 256 ${SRA_ID}${FileNamePrefix}Aligned.sortedByCoord.out.bam > tmps3.bam   # Both reads of the pair are unmapped
#samtools merge -u - tmps[123].bam | samtools sort -n -o unmapped.bam



# BAM ==> fastq
#bamToFastq -i unmapped.bam \
#    -fq ${SRA_ID}${FileNamePrefix_out}_unmapped_1.fastq \
#    -fq2 ${SRA_ID}${FileNamePrefix_out}_unmapped_2.fastq