#!/shared/ifbstor1/software/miniconda/envs/r-4.2.3/bin/Rscript

library(stringr)
library(plyr)
library(ggsankey)
library(ggplot2)
library(dplyr)


Taxonomy_small_color_palette = function(df, reference_colum){

  
  df=df[order(df[[reference_colum]]), ]
  df=df[!str_detect(df$taxon_full_name,"Other"), ]
  
  
  color_df_1 = data.frame("variable" = unique(df[str_detect(df[[reference_colum]], "Bacteria"), reference_colum  ] ) )
  color_df_1$col = colorRampPalette(c("sienna","lightyellow","steelblue3" ))(nrow(color_df_1))
  
  color_df_2 = data.frame("variable" = unique(df[
    str_detect(df[[reference_colum]], "Eukaryota") & !str_detect(df[[reference_colum]], "Fungi|myco"), reference_colum ] ) )
  
  
  color_df_2$col = colorRampPalette(c("plum1","firebrick1"))(nrow(color_df_2))
  
  color_df_3 = data.frame("variable" = unique(df[str_detect(df[[reference_colum]], "Viruses"), reference_colum ] ) )
  color_df_3$col = colorRampPalette(c("yellow2","darkorange"))(nrow(color_df_3))
  
  color_df_4 = data.frame("variable" = unique(df[str_detect(df[[reference_colum]], "Fungi|myco"), reference_colum ] ) )
  color_df_4$col = colorRampPalette(c("seagreen1","seagreen4"))(nrow(color_df_4))
  
  other_df <- data.frame(variable=c('Other Bacteria', 'Other Viruses', "Other Eukaryota","Other Fungi", "Archaea"),
                         col=c("royalblue4", "#D5FF29","firebrick4","darkgreen", "azure"))
  
  color_df = rbind(color_df_1,color_df_2,color_df_3,color_df_4,other_df)
  
  color_dict = color_df$col
  names(color_dict) = color_df$variable
  return(color_dict) }

format_human_taxon=function(classif_df){
  human_df=classif_df[classif_df$class == "Mammalia",]
  human_df$taxon="Homo_sapiens"
  return(human_df)
}

format_unclassified_taxa=function(classif_df){
  unclassified=classif_df[str_detect(classif_df$superkingdom, "unclassified"),]
  
  unclassified$taxon=ifelse(str_detect(unclassified$superkingdom, "unclassified"), 
                            "unclassified", unclassified$superkingdom)
  
  unclassified$taxon=ifelse(str_detect(unclassified$superkingdom, "other entries"), 
                            "artificial sequences", unclassified$taxon)
  unclassified$taxon=ifelse(str_detect(unclassified$superkingdom, "cellular organisms"), 
                            "cellular organisms", unclassified$taxon)
  return(unclassified)
}


format_taxa=function(classif_df, wanted_lvl){
  ranks=c("superkingdom","phylum", "class", "order", "family", "genus", "species","strain")
  
  human_df=format_human_taxon(classif_df)
  unclassified_df=format_unclassified_taxa(classif_df)
  
  classif_df=classif_df[! classif_df$read %in% c(human_df$read, unclassified_df$read),]
  
  col_to_concat=ranks[1:which(ranks==wanted_lvl)]
  classif_df$taxon=sapply(1:nrow(classif_df), function(x){paste(classif_df[x,col_to_concat], collapse = "|")})
  
  
  classif_df=classif_df[,c("read","taxon")]
  human_df=human_df[,c("read","taxon")]
  unclassified_df=unclassified_df[,c("read","taxon")]
  
  classif_merged=do.call(rbind, list(classif_df, human_df, unclassified_df) )
  classif_merged$taxon=ifelse(classif_merged$taxon == "unassembled|unassembled",
                              'Unassembled', classif_merged$taxon)
  return(classif_merged )
}

aggregate_taxa=function(table_df, n=5, replace_taxon="Bacteria"){
  top=head(table_df, n = n)
  low=table_df[!table_df$Var1 %in% top$Var1 , ]
  
  top$final_name=top$Var1
  low$final_name=sapply(low$superkingdom, function(x){paste0("Other ",replace_taxon) })
  
  table_df=rbind(top,low)
  table_df[,c("Freq","superkingdom","phylum")]=NULL
  colnames(table_df)[1]="taxon_full_name"
  
  return(table_df)
}

args = commandArgs(trailingOnly=TRUE)
 
SRA_id = args[1]
dataset = args[2]
# 
# print(dataset)


path = "/shared/projects/microbiome_translocation/"
# path = "/run/user/1001/gvfs/sftp:host=core.cluster.france-bioinformatique.fr,user=acolajanni/shared/projects/microbiome_translocation/"
#path = "/home/acolajanni/Documents/work/"

# SRA_id="SRR14418861"
# dataset="Douek_Cleveland"
# SRA_id="A21"
# dataset="Simulation"
print(dataset)

path_result=paste0(path,"results/",dataset)
path_data=paste0(path,"data/",dataset)

sra_list=readLines(paste0(path_data,"/sra_list_RNA.txt"))
sra_metadata=read.table(paste0(path_data,"/sra_metadata.tsv"),header=TRUE)

patient_id=sra_metadata[sra_metadata$Run==SRA_id,]$Title

load(paste0(path_result,"/Taxon_correspondance.Rdata"))
load(paste0(path_result,"/color_dict.Rdata"))


# result_dir_ku=paste0(path_result,"/krakenuniq/microbialDB/",SRA_id,"/")
result_dir_ku=paste0(path_result,"/kraken/nt/k2uniq/",SRA_id,"/")
result_dir_salmon=paste0(path_result,"/Contigs/",SRA_id,"/")


classif_ku=read.table(paste0(result_dir_ku,"ReadsClassif.txt"), fill = TRUE, sep="\t")
colnames(classif_ku)=c("read","taxid","superkingdom","phylum", "class", "order", "family", "genus", "species")
# If necessary, replace some characters
# classif_ku$read = sapply(classif_ku$read, function(x){str_replace_all(x, fixed("|"), fixed("-") )})
# classif_ku$read = sapply(classif_ku$read, function(x){str_split(x, fixed("/") )[[1]][1] })

classif_ku[classif_ku==""]="unclassified"

classif_ku_raw=classif_ku
classif_ku$taxid = NULL
counting=nrow(classif_ku)
reads=length(unique(classif_ku$read))


classif_sa=read.delim(paste0(result_dir_salmon,"ContigsToReads/ReadsClassif.txt"),fill=TRUE, sep='\t', 
                      header=TRUE)

colnames(classif_sa) = c("read","superkingdom","phylum", "class", "order", "family", "genus", "species","strain")
classif_sa[,c("V1","strain")]=NULL


#### Error in separatorin the dataframe: sep= " "
### problem: ex with Homo sapiens etc
### Temporarary solution: continue with genus only
classif_ku$species = NULL
classif_sa$species = NULL


length(unique(classif_sa$read))



### All the reads absent from Salmon classif are those that are not used by contigs : 
# Now directly in the file "reads classif"

# missing_reads=classif_ku$read[!classif_ku$read %in% classif_sa$read]
# 
# missing_df=classif_sa[NULL,]
# missing_df$read=NULL
# missing_df[1,]=rep("not_in_contigs",6)
# 
# missing_df[1:length(missing_reads),]=missing_df[1,]
# missing_df$read=missing_reads
# 
# classif_sa=rbind(classif_sa, missing_df)



classif_ku_phy=format_taxa(classif_ku, "phylum")
classif_sa_phy=format_taxa(classif_sa, "phylum")

# classif_sa_phy$classification="Blast-Salmon"
# classif_ku_phy$classification="Kraken2Uniq"

colnames(classif_sa_phy)[colnames(classif_sa_phy)=="taxon"]="taxon_Blast"
colnames(classif_ku_phy)[colnames(classif_ku_phy)=="taxon"]="taxon_K2Uniq"

# Bacteria | unclassified Bacteria phylum == > Bacteria | unclassified 
# classif_ku_phy$taxon_K2Uniq = ifelse(str_detect(classif_ku_phy$taxon_K2Uniq, "Bacteria|unclassified"), 
#                                          "Bacteria|unclassified", classif_ku_phy$taxon_K2Uniq )

classif_merged=merge(classif_ku_phy, classif_sa_phy, by=c("read"))


### Couting all the taxa in both classif ==> merge rarer taxa
taxa=unique(c(classif_merged$taxon_K2Uniq, classif_merged$taxon_Blast))
counting_table_ku=as.data.frame(table(classif_merged$taxon_K2Uniq))
counting_table_sa=as.data.frame(table(classif_merged$taxon_Blast))

counting_table=rbind(counting_table_ku,counting_table_sa)
counting_table=aggregate(Freq ~ Var1, data = counting_table, FUN=sum)

### Get classif
counting_table[,c("superkingdom","phylum")]=str_split_fixed(counting_table$Var1, fixed('|'), n=2)




counting_table=counting_table[order(-counting_table$Freq), ]
# Bacteria

bact=aggregate_taxa(counting_table[counting_table$superkingdom == "Bacteria", ], 6, "Bacteria")

# Viruses
vir=aggregate_taxa(counting_table[counting_table$superkingdom == "Viruses", ],5, "Viruses")

# Fungi
fungi=aggregate_taxa(counting_table[str_detect(counting_table$phylum, "myco|Fungi"), ],3 , "Fungi")

# Eukaryotes
euk=aggregate_taxa(counting_table[counting_table$superkingdom == "Eukaryota" & 
                                      !counting_table$Var1 %in% fungi$taxon_full_name , ] , 9 , "Eukaryota")


# Archeae
archaea=counting_table[counting_table$superkingdom == "Archaea",]
colnames(archaea)[1]="taxon_full_name"
archaea[,c("Freq","superkingdom","phylum")]=NULL
archaea$final_name="Archaea"

# Others
key_from_data=do.call(rbind, list(bact,vir,euk,fungi, archaea))

other=counting_table[!counting_table$Var1 %in% key_from_data$taxon_full_name,]
colnames(other)[1]=c("taxon_full_name")

other$final_name=ifelse(str_detect(as.vector(other$taxon_full_name), "unclassified"), "unclassified", as.vector(other$taxon_full_name))
other[,c("Freq","superkingdom","phylum")]=NULL
key_from_data=rbind(key_from_data, other)



### Get color dict
color_dict=Taxonomy_small_color_palette(key_from_data, reference_colum = "final_name")
missing_names=as.vector(key_from_data$final_name[!key_from_data$final_name %in% names(color_dict)])
unique(missing_names)

color_dict[["Not in Contigs"]]="red4"
color_dict[["Homo_sapiens"]]="purple"
color_dict[["unclassified"]]="#777777"
color_dict[["cellular organisms"]]="#777777"
color_dict[["artificial sequences"]]="#777777"
color_dict[["Unassembled"]]="#444444"

#color_dict[["unclassified root superkingdom"]]="green"

dir.create(paste0(path,"results/",dataset,"/rdata/"), showWarnings = FALSE)

save(classif_merged, file=paste0(path,"results/",dataset,"/rdata/",SRA_id,"_classif_merged.rdata") ) 

stop()


df <- classif_merged %>%
  make_long(taxon_Blast, taxon_K2Uniq)

p1=ggplot(df, aes(x = x, 
                  next_x = next_x, 
                  node = node, 
                  next_node = next_node,
                  fill = factor(node) )) + 
  geom_sankey(flow.alpha=0.85,node.color=1,space=150) +
  scale_fill_manual(values = color_dict) + guides(fill = guide_legend(ncol = 1))+
  # scale_fill_discrete(guide="none") +
  xlab(paste("Patient",patient_id))+#  scale_x_discrete(labels = c("Blue pipeline","Red pipeline")) +
  theme_sankey(base_size = 16, ) 

p1


ggsave(p1,filename=paste0(path_result,"/sankey_plots/",SRA_id,"_sankey_truth.png"),
       device = "png",height = 36, width = 60, units = "cm")


# Reducing reads to be displayed: # See how interesting sequences detected by k2 are classified
reduced_classif=classif_merged[!classif_merged$taxon_K2Uniq %in% c(#"unclassified","artificial sequences","cellular organisms", 
                                                                   "Homo_sapiens",
                                                                   "Eukaryota|Chordata","Eukaryota|unclassified Eukaryota phylum",
                                                                   "Eukaryota|Arthropoda","Eukaryota|Streptophyta","Eukaryota|Mollusca") ,]

# Merge the two dataframes based on the matching column
reduced_classif <- merge(reduced_classif, key_from_data, by.x = "taxon_K2Uniq", by.y = "taxon_full_name", all.x = TRUE)
reduced_classif$taxon_K2Uniq=NULL
colnames(reduced_classif)[colnames(reduced_classif) == "final_name"]="taxon_K2Uniq"

reduced_classif <- merge(reduced_classif, key_from_data, by.x = "taxon_Blast", by.y = "taxon_full_name", all.x = TRUE)
reduced_classif$taxon_Blast=NULL
colnames(reduced_classif)[colnames(reduced_classif) == "final_name"]="taxon_Blast"



df2 <- reduced_classif %>%
  make_long(taxon_Blast, taxon_K2Uniq)

p2=ggplot(df2, aes(x = x, 
                  next_x = next_x, 
                  node = node, 
                  next_node = next_node,
                  fill = factor(node) )) + 
  geom_sankey(flow.alpha=0.85,node.color=1,space=150) +
  scale_fill_manual(values = color_dict) + guides(fill = guide_legend(ncol = 1))+
  # scale_fill_discrete(guide="none") +
  xlab(paste("Patient",patient_id))+#  scale_x_discrete(labels = c("Blue pipeline","Red pipeline")) +
  theme_sankey(base_size = 16, ) 

p2

ggsave(p2,filename=paste0(path_result,"/sankey_plots/",SRA_id,"_sankey_no_eukaryotes.png"),
       device = "png",height = 36, width = 60, units = "cm")

# Reducing reads to be displayed: # See how interesting sequences detected by k2 are classified
reduced_classif=classif_merged[!classif_merged$taxon_Blast %in% c("unclassified","artificial sequences","cellular organisms", "Homo_sapiens",
                                                                   "Eukaryota|Chordata","Eukaryota|unclassified Eukaryota phylum",
                                                                   "Eukaryota|Arthropoda","Eukaryota|Streptophyta","Eukaryota|Mollusca",
                                                                   "Not in Contigs") ,]

# Merge the two dataframes based on the matching column
reduced_classif <- merge(reduced_classif, key_from_data, by.x = "taxon_K2Uniq", by.y = "taxon_full_name", all.x = TRUE)
reduced_classif$taxon_K2Uniq=NULL
colnames(reduced_classif)[colnames(reduced_classif) == "final_name"]="taxon_K2Uniq"

reduced_classif <- merge(reduced_classif, key_from_data, by.x = "taxon_Blast", by.y = "taxon_full_name", all.x = TRUE)
reduced_classif$taxon_Blast=NULL
colnames(reduced_classif)[colnames(reduced_classif) == "final_name"]="taxon_Blast"



df3 <- reduced_classif %>%
  make_long(taxon_Blast, taxon_K2Uniq)

p3=ggplot(df3, aes(x = x, 
                   next_x = next_x, 
                   node = node, 
                   next_node = next_node,
                   fill = factor(node) )) + 
  geom_sankey(flow.alpha=0.85,node.color=1,space=150) +
  scale_fill_manual(values = color_dict) + guides(fill = guide_legend(ncol = 1))+
  # scale_fill_discrete(guide="none") +
  xlab(paste("Patient",patient_id))+#  scale_x_discrete(labels = c("Blue pipeline","Red pipeline")) +
  theme_sankey(base_size = 16, ) 

p3

ggsave(p3,filename=paste0(path_result,"/sankey_plots/",SRA_id,"test2.png"),
       device = "png",height = 36, width = 60, units = "cm")

################################################################################

### Trying rebuilding quant_lvl item: quantification for each level of classification
ranks=c("superkingdom","phylum", "class", "order", "family", "genus")


### 1: Make a Function
### 2: Maybe flag human reads ?
classif_df=classif_ku
quant_lvl=list()
for (lvl in ranks){
  print(lvl)
  lvl_df=classif_df
  
  col_to_concat=ranks[1:which(ranks==lvl)]
  
  if(lvl == "superkingdom"){
    lvl_df$taxon_full_name = lvl_df[[lvl]]
  }else{
    lvl_df$taxon_full_name=sapply(1:nrow(lvl_df), function(x){paste(lvl_df[x,col_to_concat], collapse = "|")})
  }
  lvl_df$taxon=lvl_df[[lvl]]

  
  lvl_df=lvl_df[,c("taxon","taxon_full_name")]
  
  # Counting occurences 
  lvl_df <- lvl_df %>%
    group_by(across(everything())) %>%  # Group by all columns
    summarise(NumReads = n()) %>%          # Count the occurrences of each row
    ungroup()
  

  lvl_df$tax_lvl = lvl
  lvl_df$readFreq = lvl_df$NumReads / sum(lvl_df$NumReads) 
  
  quant_lvl[[lvl]] = lvl_df
}









