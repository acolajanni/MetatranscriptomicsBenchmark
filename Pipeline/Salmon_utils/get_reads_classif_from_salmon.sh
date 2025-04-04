#!/bin/bash


PATH_RES="${1:-/shared/projects/microbiome_translocation/results/Simulation_old/Contigs/}"
SRA_ID="${2:-A5}"

PATH_RESULTS=${PATH_RES}${SRA_ID}/

#path_results=/home/acolajanni/Documents/work/fastq_scripts/test_salmon_toreadsclassif/


TAXONOMY=( strain species genus family order class phylum superkingdom )

index=0
ncontig_max=0
index_to_keep=0

different_lowest_lvl=( )

# Counting contigs in each file
# Get files where there are new contigs not in the precedent 
for lvl in ${TAXONOMY[@]} ; do  
    filename=contig_ID_${lvl}.txt
    ncontig=$(cat ${PATH_RESULTS}Quantification/${filename} | wc -l )

    if [[ $ncontig -gt $ncontig_max ]] ; then

        ncontig_max=$ncontig        
        ## Add level where there is more contigs classified
        # always 'strain' lvl 
        different_lowest_lvl+=( $lvl ) 

    elif [[ $ncontig -eq $ncontig_max ]] ; then
        echo 'same'
    fi

    index=$(($index+1))
    

done


echo ${different_lowest_lvl[@]}
printf "%s\n" "${different_lowest_lvl[@]}" > ${PATH_RESULTS}Quantification/tax_lvl_to_map.txt



#### Fetch reads associated with each contigs    
# @1: path to results (~results/dataset/Contigs/ID_XXX/)
# @2: IDentifiant
# @3: taxonomic lvl
# @4: Files with contig ID
# bash /shared/projects/microbiome_translocation/fastq_scripts/Salmon_utils/retrieveReadsFromContigs.sh \
#     $PATH_RES $SRA_ID strain $SRA_ID ${PATH_RESULTS}Quantification/contig_ID_strain.txt


# depending on files detected, fetch 'new' contigs 
if [[ ${#different_lowest_lvl[@]} -ge 2 ]] ; then

    length=$((${#different_lowest_lvl[@]} - 1))
    # loop through indices
    for i in ${!different_lowest_lvl[@]} ; do 

        if [[ $i -eq $length ]] ; then 
            continue
        fi

        next_i=$(( $i + 1 ))


        next_lvl=${different_lowest_lvl[$next_i]}
        lvl=${different_lowest_lvl[$i]}

        # Keep new contig id found at the new tax lvl
        sort ${PATH_RESULTS}Quantification/contig_ID_$lvl.txt ${PATH_RESULTS}Quantification/contig_ID_$next_lvl.txt \
            | uniq -u > ${PATH_RESULTS}Quantification/contig_ID_diff_$next_lvl.txt

        echo $next_lvl
        echo "next script ----------------------------"

        #### Fetch reads associated with each contigs    
        # @1: path to results (~results/dataset/Contigs/ID_XXX/)
        # @2: IDentifiant
        # @3: taxonomic lvl
        # @4: Files with contig ID
        bash /shared/projects/microbiome_translocation/fastq_scripts/Salmon_utils/retrieveReadsFromContigs.sh \
            $PATH_RES $SRA_ID $next_lvl ${PATH_RESULTS}Quantification/contig_ID_diff_$next_lvl.txt


        echo ____________________________ sort uniq _______________________
        echo $lvl
        echo $next_lvl
        echo _______________________
        # Sort uniq
    done

fi

# Retrieve all the pairs reads ==> contigs and dedplicate them
cat ${PATH_RESULTS}ContigsToReads/*_readsToContigs.txt \
    | awk '!a[$0]++' > ${PATH_RESULTS}ContigsToReads/reads_to_contig.txt


### Handle Unclassied and Human reads
#### Fetch reads associated with each contigs    
# @1: path to results (~results/dataset/Contigs/ID_XXX/)
# @2: IDentifiant
# @3: taxonomic lvl
bash /shared/projects/microbiome_translocation/fastq_scripts/Salmon_utils/retrieveReadsFromContigs.sh \
    $PATH_RES $SRA_ID unclassified $SRA_ID 
    
bash /shared/projects/microbiome_translocation/fastq_scripts/Salmon_utils/retrieveReadsFromContigs.sh \
    $PATH_RES $SRA_ID human $SRA_ID