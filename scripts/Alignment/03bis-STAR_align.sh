#! /bin/sh

################################################################################
# > september 2023
# > Script : STAR build align.sh
# > Function : Using STAR aligner to build genome index and align reads to it
# @ COLAJANNI Antonin
################################################################################


PATH_RES="${1:-/shared/projects/microbiome_translocation/data/Douek_cell2021/}"
SEQ_TYPE="${2:-PE}"
ID="${3:-/shared/projects/microbiome_translocation/data/Douek_cell2021/sra_subset_in_manuscript.txt}"
GENOME_input="${4:-hg19}"
PATH_INPUT="${5:-/shared/projects/microbiome_translocation/data/Douek_cell2021/Bowtie2_mapping/}"
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

    else 
        echo - Single-end -
        # map single-end reads to genome
        STAR --runThreadN 20 \
            --readFilesIn ${PATH_INPUT}/${SRA_ID}/${SRA_ID}${FileNamePrefix_in}".fastq" \
            --genomeDir $STAR_index \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${SRA_ID}${FileNamePrefix_out}_mapping_ \
            --outReadsUnmapped Fastx \
            --outSAMstrandField intronMotif 1>star_stdout.txt 2>star_stderr.txt # strand derived from the intron motif. Reads with inconsistent and/or non-canonical introns are filtered out.

        FileNamePrefix=_STAR_

        mv ${SRA_ID}${FileNamePrefix_out}_mapping_Unmapped.out.mate1 ${SRA_ID}${FileNamePrefix}unmapped.fastq 

    fi
    echo " --- "


done

echo "----------------------------------"
echo " STAR - END "

