library(stringr)
library(dplyr)

db_path="/shared/projects/microbiome_translocation/"


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

#Remove annotation that are too recent compared to database (nt database in 11 april 2024) 
### annotation date but maybe looking at "seq_rel_date"
refseq_full = refseq_full %>% filter(annotation_date <= "2024/04/11")


### Bacteria ### 
HumGut_bact=read.delim(paste0(db_path,"database_clean/selected_genomes/Bacteria/HumGut_taxid_lineage.tsv"),
                       col.names = c("taxid", "remove", "taxon_full_name"), fill = TRUE)

HumGut_bact$remove=NULL
HumGut_bact[c("D","P","C","O","F","G","S",'T')] = str_split_fixed(HumGut_bact$taxon_full_name, fixed('|'), 8)

refseq_bact = refseq_full[refseq_full$S %in% HumGut_bact$S ,]


length(unique(refseq_bact$S))



### Fungi
FUNOMIC_fungi=read.delim(paste0(db_path,"database_clean/selected_genomes/Fungi/FUNOMIC_taxid_lineage.tsv"),
                         col.names = c("taxid", "remove", "taxon_full_name"), fill = TRUE)


FUNOMIC_fungi$remove=NULL
FUNOMIC_fungi[c("D","P","C","O","F","G","S",'T')] = str_split_fixed(FUNOMIC_fungi$taxon_full_name, fixed('|'), 8)

refseq_fungi = refseq_full[refseq_full$S %in% FUNOMIC_fungi$S ,]

### Other Eukaryotes

eupathdb=read.delim(paste0(db_path,"database_clean/selected_genomes/other_eukaryotes/EuPathDB54_taxid_lineage.tsv"),
                    col.names = c("taxid", "remove", "taxon_full_name"), fill = TRUE)


eupathdb$remove=NULL
eupathdb[c("D","P","C","O","F","G","S",'T')] = str_split_fixed(eupathdb$taxon_full_name, fixed('|'), 8)

refseq_eukaryotes = refseq_full[refseq_full$S %in% eupathdb$S ,]

### Viruses
known_vir=readLines(paste0(db_path,"database_clean/selected_genomes/Viruses/viral_taxa_Lian2021_Bai2022.txt"))
known_vir=known_vir[known_vir != "Virgaviridae"]


remove_nonhuman_vir=function(refseq_df, regex, lvl){
  
  refseq_df2=refseq_df[refseq_df[[lvl]] == regex,]
  refseq_OG=refseq_df[!refseq_df[[lvl]] == regex,]
  
  refseq_df2=refseq_df2[str_detect(refseq_df2$S, "Human") | str_detect(refseq_df2$S, "human"),]
  
  return(rbind(refseq_OG, refseq_df2))
}

refseq_vir=refseq_full[refseq_full$F %in% known_vir | refseq_full$C %in% known_vir ,]


refseq_vir=remove_nonhuman_vir(refseq_vir, "Retroviridae", "F" )
refseq_vir=remove_nonhuman_vir(refseq_vir, "Adenoviridae", "F" )
refseq_vir=remove_nonhuman_vir(refseq_vir, "Lymphocryptovirus", "G" )
refseq_vir=remove_nonhuman_vir(refseq_vir, "Proboscivirus", "G" )
refseq_vir=remove_nonhuman_vir(refseq_vir, "Simplexvirus", "G" )
refseq_vir=remove_nonhuman_vir(refseq_vir, "Circoviridae", "F" )


# Found only in plants
refseq_vir=refseq_vir[refseq_vir$G != "Hordeivirus", ]

to_remove=c("Bovine","Equine", "Canine", "Pig", "Ovine", "Duck", "Fowl", "Simian", "Turkey", "Murine", 
            "Deer","Frog","Snake","Bat","Falcon","Goose","Great tit", "Guinea pig","Lizard","Penguin",
            "Pigeon","Platyrrhini", "Polar Bear", "Pond slider","Psittacine", "Possum","Raptor","Rousettus",
            "Sea lion", "Simian", "Skua", "Skunk", "Squirrel", "Sturgeon", "Tree shrew", "Turkey",
            "feline", "Caprine", "Baboon", "sheep", "Gibbon", "Avian", "Rabbit", "Monkey", "Cucumber",
            "Tomato", "Tobacco", "Pea", "Beet", "parrot","swine" ,"Mink" ,"Equid" ,"Otarine" ,"seal", "Chicken")


pattern <- paste(to_remove, collapse = "|")
refseq_vir=refseq_vir %>%  filter(!grepl(pattern, S, ignore.case = TRUE))


refseq_vir_no_caudo=refseq_vir[refseq_vir$P != "Uroviricota",]

# Useless: ecoli phage = caudoviricetes (in data "refseq_vir")
#refseq_ecoliphage=refseq_full[str_detect(refseq_full$S, "phage") & str_detect(refseq_full$S, "Coli|coli|Escherichia") ,]

refseq_selected_sequences=do.call(rbind, list(refseq_vir, refseq_bact, refseq_eukaryotes, refseq_fungi))

refseq_selected_sequences$ftp_path = sapply(refseq_selected_sequences$ftp_path,
                                            function(x){ paste0(str_replace(x,"https://","ftp://"),"/")   })

refseq_selected_sequences=refseq_selected_sequences[refseq_selected_sequences$S != "",]
refseq_selected_sequences=unique(refseq_selected_sequences)
### For species with multiple genome representation, take only complete genomes
sp_selected=as.data.frame(table(refseq_selected_sequences$S))

duplicated_sp=sp_selected[sp_selected$Freq > 1,]$Var1
duplicated=sp_selected[sp_selected$Freq > 1,]
single_sp=sp_selected[sp_selected$Freq <= 1,]$Var1
single=sp_selected[sp_selected$Freq == 1,]


refseq_selected_sequences$Full_id=sapply(refseq_selected_sequences$ftp_path, function(x){
  tmp=str_split(x,fixed("/"))[[1]]
  tmp=tmp[tmp != ""]
  return(tmp[length(tmp)])
})


refseq_single=refseq_selected_sequences[refseq_selected_sequences$S %in% single_sp,]
n_replicate=10
random_drawing=list()

### 10 different drawing of 1 genome per species
for (i in c(0:(n_replicate-1))) {
  set.seed(i)
  print(i)
  
  df_rep=refseq_single
  
  rep=lapply(duplicated_sp, function(sp){
    tmp_df=refseq_selected_sequences[refseq_selected_sequences$S == sp ,]
    tmp_df=tmp_df[sample(1:nrow(tmp_df),1),]
    return(tmp_df)
  })
  rep[["single"]]=refseq_single
  df_rep=do.call(rbind, rep)
  
  
  random_drawing[[ toString(i) ]] = df_rep 
  
}

random_drawing=lapply(random_drawing, function(x){
  x$P = str_replace_all(x$P, " ", "-")
  return(x)
})




save(random_drawing, file = "/shared/projects/microbiome_translocation/database_clean/selected_genomes/selected_genomes_list.RData")


distant_path="/run/user/1001/gvfs/sftp:host=core.cluster.france-bioinformatique.fr,user=acolajanni/shared/projects/microbiome_translocation/"

lapply(names(random_drawing), function(rep){
  x=random_drawing[[rep]]
  print(rep)
  current_path=paste0(distant_path, "database_clean/")
  dir.create(paste0(current_path,rep),showWarnings = FALSE)
  write.table(x,paste0(current_path,rep,"/selected_genomes.tsv"), 
              sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  } )


stop()

selected_sequences_unique=unique(do.call(rbind,random_drawing)) 

write.table(selected_sequences_unique, file = paste0(distant_path, "database_clean/refseq_selected_genomes.tsv"),
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)


nrow_to_take=round(nrow(selected_sequences_unique)/20)+1
nr <- nrow(selected_sequences_unique)
tmp=split(selected_sequences_unique, rep(1:ceiling(nr/nrow_to_take), each=nrow_to_take, length.out=nr))

for (i in 1:length(tmp)){
  
  df=tmp[[i]]
  write.table(df, file = paste0(distant_path, "database_clean/refseq_selected_genomes_p",i,".tsv"),
              sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  
  
}



# refseq_filtered=refseq_selected_sequences[refseq_selected_sequences$S %in% single_sp,]
# pref_assembly=c("Complete Genome","Chromosome","Scaffold","Contig")
# for (sp in duplicated_sp){
#   
#   tmp_df=refseq_selected_sequences[refseq_selected_sequences$S==sp,]
#   unique(tmp_df$assembly_level)
#   
#   # Keep preferred assembly level genome > Chromosome > scaffold > Contig
#   index_pref=pref_assembly[min(which(pref_assembly %in% unique(tmp_df$assembly_level)))]
#   
#   tmp_df=tmp_df[tmp_df$assembly_level == index_pref,]
#   
#   refseq_filtered=rbind(refseq_filtered, tmp_df)
# } 


length(unique(refseq_filtered$S))
t=as.data.frame(table(refseq_filtered$S))

refseq_selected_sequences=refseq_filtered

nrow_to_take=round(nrow(refseq_selected_sequences)/20)+1
nr <- nrow(refseq_selected_sequences)
tmp=split(refseq_selected_sequences, rep(1:ceiling(nr/nrow_to_take), each=nrow_to_take, length.out=nr))

distant_path="/run/user/1001/gvfs/sftp:host=core.cluster.france-bioinformatique.fr,user=acolajanni/shared/projects/microbiome_translocation/"




refseq_selected_sequences$P=sapply(refseq_selected_sequences$P, function(x){ str_replace_all(x, " ", "-")} )

write.table(refseq_selected_sequences, file = paste0(distant_path, "database_clean/refseq_selected_genomes.tsv"),
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

for (i in 1:length(tmp)){
  
  df=tmp[[i]]
  write.table(df, file = paste0(distant_path, "database_clean/refseq_selected_genomes_p",i,".tsv"),
              sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  
  
}



table(refseq_selected_sequences$D)



test=random_drawing
test=unique(do.call(rbind, test))
test=random_drawing[[1]]
lapply(random_drawing, function(x) {length(unique(x$S))})


table=as.data.frame(table(test$D))
table=table[table$Freq > 0,]




load( "/shared/projects/microbiome_translocation/database_clean/selected_genomes/selected_genomes_list.RData")
test=lapply(random_drawing[1], function(x){
  t=as.data.frame(table(x$F, x$G) ) 
  t=t[t$Freq > 1 ,]
  return(t)
  })
test=test[["0"]]
table(test$Var1)
t2=aggregate(Var2 ~ Var1, test, length)
t2=t2[t2$Var2 > 1 ,]

table(t2$Var1)



sel=random_drawing[["0"]]
sel=sel[sel$F %in% t2$Var1 , ]


refseq_single=refseq_selected_sequences[refseq_selected_sequences$S %in% single_sp,]
n_replicate=10
random_drawing=list()

n_replicate=10
random_drawing=list()
results <- parallel::mclapply(0:(n_replicate-1), function(i) {
  set.seed(i)
  print(i)
  
  df_rep <- refseq_single
  
  rep <- lapply(duplicated_sp, function(sp) {
    tmp_df <- refseq_selected_sequences[refseq_selected_sequences$S == sp, ]
    tmp_df <- tmp_df[sample(1:nrow(tmp_df), 1), ]
    return(tmp_df)
  })
  
  rep[["single"]] <- refseq_single
  df_rep <- do.call(rbind, rep)
  
  return(list(name = toString(i), data = df_rep))
}, mc.cores = 10)

# Store results in random_drawing list
random_drawing <- setNames(lapply(results, `[[`, "data"), sapply(results, `[[`, "name"))

### Taking only families with at least 2 genus
counting_family=parallel::mclapply(random_drawing, function(x){
  # Get unique genus in data
  t=unique(x[,c("F","G")])
  
  # Counting number of different genus per family
  t2=aggregate(G ~ F, t, length)
  t2=t2[t2$G > 1 ,]
  return(t2$F)

},mc.cores = 10)

Family_drawing <- parallel::mclapply(0:(n_replicate-1), function(i) {

  families=counting_family[[as.character(i)]]
  selected_df=random_drawing[[as.character(i)]]
  
  selected_df[selected_df$F %in% families , ]
  
  return(selected_df[selected_df$F %in% families , ])
  
  
}, mc.cores = 10)
names(Family_drawing) = c(0:9)


Family_drawing=lapply(Family_drawing, function(df){ 
  df$F = sapply(df$F, function(x) str_replace_all(x, " ", "-") )
  return(df)
  } )


lapply(names(Family_drawing), function(rep){
  x=Family_drawing[[rep]]
  print(rep)
  current_path=paste0(db_path, "database_clean/")
  dir.create(paste0(current_path,rep),showWarnings = FALSE)
  write.table(x,paste0(current_path,rep,"/selected_genomes_family.tsv"), 
              sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
  
} )



