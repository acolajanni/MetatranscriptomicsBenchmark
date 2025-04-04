#!/bin/bash

database_path=/shared/projects/microbiome_translocation/database/
# cd /shared/projects/microbiome_translocation/database/
# cd Simulation_database/

taxid_array=$(awk -F'\t' '{print $2}' ${database_path}Sequence_database.tsv)

file_array=$(awk -F'\t' '{print $3}' ${database_path}Sequence_database.tsv)

# change into a bash array
file_array=(`echo ${file_array}`)

counter=1
for i in $taxid_array ; do 
    taxid=$i
    file=${file_array[$counter]}
    echo $file

    filepath=${database_path}Simulation_database/${file}
    filepath_out=${database_path}Simulation_format/${file}
    # Remove .gz because seqkit will unzip the file 
    filepath_out=${filepath_out%.gz}

    # Add the taxid to the sequence identifiant
    seqkit replace -p '(.+)' -r taxid:${taxid}_'${1}' $filepath > ${filepath_out}



    counter=$((counter+1))
done

