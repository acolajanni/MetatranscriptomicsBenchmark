# Select a subset of file

library(stringr)

path = "~/"

path_data = file.path(path,"data/Douek_cell2021/")
path_data_DNA = file.path(path,"data/Douek_cell2021_DNA/")
path_data_2 = file.path(path,"data/Douek_Montreal/")
path_data_3 = file.path(path,"data/Douek_Cleveland/")


metadata = read.csv2(file.path(path_data,"metadata.csv"),header = TRUE, sep=",")

table(metadata$subject_classification)                      
table(metadata$source_name)                      


metadata_subset = metadata[metadata$source_name == "Plasma" &
                             metadata$subject_classification == "Pre_ART_HIV_Subject" & 
                             metadata$LibrarySource == "TRANSCRIPTOMIC" &
                             metadata$cd8_count == "Data_in_Manuscript" ,]


metadata_subset2 = metadata[metadata$source_name == "Plasma" &
                             metadata$subject_classification == "Pre_ART_HIV_Subject" & 
                             metadata$LibrarySource == "GENOMIC" & 
                             metadata$cd8_count == "Data_in_Manuscript",]


metadata_subsetfull = metadata[metadata$source_name == "Plasma" &
                              metadata$subject_classification == "Pre_ART_HIV_Subject" & 
                              metadata$LibrarySource %in% c("GENOMIC","TRANSCRIPTOMIC") & 
                              metadata$cd8_count == "Data_in_Manuscript",]


## Montreal anbd cleveland cohort
# Correspondance SRA ID - paper ID
id_correspondance=read.csv( file.path(path_data,"metadata_metagenome.csv"))
id_correspondance=id_correspondance[,c("Title","Accession")]


metadata=merge(id_correspondance,metadata_subsetfull, by="Accession")

write.table(metadata[,c("indiv", "Run")],file = paste0(path,"data/Douek_cell2021/key_correspondance.txt"),quote = FALSE, row.names=FALSE,sep = "\t")


metadata=merge(id_correspondance,metadata, by.y="Sample.Name", by.x = "Accession")
metadata$indiv=unlist(lapply(metadata$Title, function(x){return(str_remove(x, "DNA_|RNA_"))}))

metadata_subset3 = metadata[metadata$source_name == "Plasma" &
                            metadata$geographical_location == "Canada: Montreal",]

metadata_subset4 = metadata[metadata$source_name == "Plasma" &
                              metadata$geographical_location == "USA: Cleveland, Ohio",]




table(metadata_subset3$indiv)
table(metadata_subset4$indiv)


# metadata_subset3 = metadata[metadata$source_name == "Plasma" &
#                               metadata$geographical_location == "USA: Cleveland, Ohio",]
# 
# table(metadata_subset3$indiv)

# 26 SRA file
writeLines(metadata_subset$Run, file.path(path_data, "sra_subset.txt"))


# 11 in Manuscript 
metadata_subset = metadata[metadata$source_name == "Plasma" &
                             metadata$subject_classification == "Pre_ART_HIV_Subject" & 
                             metadata$LibrarySource == "TRANSCRIPTOMIC" &
                             metadata$cd8_count == "Data_in_Manuscript" ,]


write.csv(metadata_subset,  file.path(path_data, "sra_metadata_subset.csv") )
writeLines(metadata_subset$Run, file.path(path_data, "sra_subset_in_manuscript.txt"))
writeLines(metadata_subset2$Run, file.path(path_data_DNA, "sra_subset_in_manuscript.txt"))
writeLines(metadata_subset3$Run, file.path(path_data_2, "sra_subset_in_manuscript.txt"))
writeLines(metadata_subset4$Run, file.path(path_data_3, "sra_list.txt"))

write.table(metadata_subset3[,c("Title","Run")],  file.path(path_data_2, "sra_metadata.tsv"), quote=F, row.names=FALSE )
write.table(metadata_subset4[,c("Title","Run")],  file.path(path_data_3, "sra_metadata.tsv"), quote=F, row.names=FALSE,sep = "\t" )



# Retrieve all times SRA id
tmp = metadata[metadata$Organism == "human blood metagenome" & metadata$cd8_count == "Data_in_Manuscript" ,]
write.csv(tmp,  file.path(path_data, "metadata_alltimes_inpaper.csv") )

# Correspondance SRA ID - paper ID
id_correspondance=read.csv( file.path(path_data,"metadata_metagenome.csv"))
id_correspondance=id_correspondance[,c("Title","Accession")]


merged_table=merge(tmp,id_correspondance, by.x="Sample.Name", by.y = "Accession")
write.csv(merged_table,  file.path(path_data, "metadata_subset_ID.csv") )
