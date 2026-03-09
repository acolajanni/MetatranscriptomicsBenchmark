library(stringr)
library(dplyr)

db_path="/home/acolajanni/Documents/work/database/"
db_path="~//"


refseq=read.delim(paste0(db_path,"database_clean/assembly_summary_refseq.txt"), 
                  fill = TRUE, as.is = TRUE, skip=1, quote="")

refseq_lineage=read.delim(paste0(db_path,"database_clean/refseq_summary_lineage.tsv"),
                          col.names = c("taxid","taxon_full_name0","taxon_full_name"))

refseq_lineage$taxon_full_name0 = NULL

refseq_lineage[c("D","P","C","O","F","G","S",'T')] = str_split_fixed(refseq_lineage$taxon_full_name, fixed('|'), 8)
refseq_lineage$taxon_full_name = NULL
refseq_full=merge(refseq, refseq_lineage, by="taxid")
refseq_full[,c("bioproject","refseq_category" ,"biosample","wgs_master","isolate","version_status","release_type",
               "asm_name","asm_submitter","gbrs_paired_asm","paired_asm_comp","relation_to_type_material",
               "asm_not_live_date","assembly_type","genome_size","genome_size_ungapped","gc_percent",
               "replicon_count","scaffold_count","contig_count","annotation_provider",
               "annotation_name","pubmed_id")] = NULL

#Remove annotation that are too recent compared to databases (nt database in 11 april 2024) 
### annotation date but maybe looking at "seq_rel_date"
refseq_full = refseq_full %>% filter(annotation_date <= "2024/04/11")



################################################################################

realist_simu_path="~//data/Simulation/realist/"

genus_to_select=read.table(paste0(realist_simu_path, "genus_to_select_lineage.txt"), sep="\t",
                           col.names = c("G", "taxid_g","K","P","C","O","F"))

genus_to_select = genus_to_select[genus_to_select$G != "Homo", ]



refseq_selected_sequences=refseq_full[refseq_full$G %in% genus_to_select$G , ]
table(refseq_selected_sequences$G)

################################################################################


refseq_selected_sequences$ftp_path = sapply(refseq_selected_sequences$ftp_path,
                                            function(x){ paste0(str_replace(x,"https://","ftp://"),"/")   })

refseq_selected_sequences=refseq_selected_sequences[refseq_selected_sequences$S != "",]
refseq_selected_sequences=unique(refseq_selected_sequences)

### For species with multiple genome representation, take only complete genomes
g_selected=as.data.frame(table(refseq_selected_sequences$G))

# duplicated_sp=sp_selected[sp_selected$Freq > 1,]$Var1
# duplicated=sp_selected[sp_selected$Freq > 1,]
# single_sp=sp_selected[sp_selected$Freq <= 1,]$Var1
# single=sp_selected[sp_selected$Freq == 1,]


refseq_selected_sequences$Full_id=sapply(refseq_selected_sequences$ftp_path, function(x){
  tmp=str_split(x,fixed("/"))[[1]]
  tmp=tmp[tmp != ""]
  return(tmp[length(tmp)])
})


# refseq_single=refseq_selected_sequences[refseq_selected_sequences$S %in% single_sp,]
n_replicate=50
random_drawing=list()

for (i in c(0:(n_replicate-1))) {
  
  set.seed(i)
  
  rep=parallel::mclapply(g_selected$Var1, function(g){
    tmp_df=refseq_selected_sequences[refseq_selected_sequences$G == g  ,]
    
    # Favor genomes not excluded from refseq, and complete genomes
    testing_genome=nrow(tmp_df[tmp_df$excluded_from_refseq == "na",])
    if(testing_genome > 0){
      tmp_df = tmp_df [ tmp_df$excluded_from_refseq == "na", ] }
    
    testing_genome=nrow(tmp_df[tmp_df$assembly_level == "Complete Genome",])
    if(testing_genome > 0){
      tmp_df = tmp_df [ tmp_df$assembly_level == "Complete Genome", ] }
    
    # Take one randomly
    tmp_df=tmp_df[sample(1:nrow(tmp_df),1),]
    return(tmp_df)
  },mc.cores = 12)
  rep=ldply(rep)
  
  random_drawing[[ toString(i) ]] = rep
}

save(random_drawing, file = "~//database_clean/selected_genomes/realist/selected_genomes_list.RData")
load(file = "~//database_clean/selected_genomes/realist/selected_genomes_list.RData")


#### Save each dataframe 


dir_to_save = "~//database_clean/selected_genomes/realist/"
n_replicate=50
for (i in c(0:(n_replicate-1))){
  print(i)
  current_path=paste0(dir_to_save, toString(i),"/" )
  dir.create(current_path, showWarnings = FALSE)
  
  write.table(random_drawing[[toString(i)]],paste0(current_path,"selected_genomes_realist.tsv"), 
              sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
  
}


full_seq=ldply(random_drawing)
full_seq$.id = NULL
full_seq=unique(full_seq)
write.table(full_seq , paste0(dir_to_save,"selected_genomes_realist.tsv"),
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
