#!/shared/ifbstor1/software/miniconda/envs/r-4.2.3/bin/Rscript

library(stringr)
library(plyr)
library(ggsankey)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)

###############################################################################
# F U N C T I O N S
################################################################################
ranks=c("superkingdom","phylum", "class", "order", "family", "genus", "species")

aggregate_taxa_count_df=function(table_df, n=5, column_to_filter, Frequency_column,replace_taxon="Bacteria"){
  table_df=table_df[table_df$superkingdom == replace_taxon , ]
  table_df=table_df[order( table_df[[Frequency_column]],decreasing = TRUE ) , ]
  
  if (replace_taxon == "Eukaryota"){
    
    print("eukaryota")
    fungi_df=table_df[str_detect(table_df[[column_to_filter]], "Fungi|fungi|myco" ),]
    fungi_df$superkingdom = "Eukaryota - Fungi"
    
    fun=aggregate_taxa_count_df(fungi_df, 2, column_to_filter, Frequency_column, "Eukaryota - Fungi"  )
    table_df=table_df[! table_df[[column_to_filter]] %in% fungi_df[[column_to_filter]] , ]    
    
  }
  
  top=head(table_df, n = n)
  low=table_df[!table_df[[column_to_filter]] %in% top[[column_to_filter]] , ]
  
  top$final_name=top[[column_to_filter]]
  low$final_name=sapply(low$superkingdom, function(x){paste0("Other ",replace_taxon) })
  
  table_df=rbind(top,low)
  table_df[,c("superkingdom","phylum")]=NULL
  colnames(table_df)[1]="taxon_full_name"
  
  if (replace_taxon == "Eukaryota"){
    table_df=rbind(fun,table_df)
  }
  
  
  return(table_df)
}


Taxonomy_small_color_palette = function(df, reference_colum){
  
  df=df[order(df[[reference_colum]]), ]
  df=df[!str_detect(df[[reference_colum]],"Other"), ]
  
  color_df_1 = data.frame("variable" = unique(df[str_detect(df[[reference_colum]], "Bacteria"), reference_colum  ] ) )
  color_df_1$col = colorRampPalette(c("aliceblue","royalblue4"  ))(nrow(color_df_1))
  
  color_df_2 = data.frame("variable" = unique(df[str_detect(df[[reference_colum]], "Eukaryota") , reference_colum ] ) )
  #color_df_2$col = colorRampPalette(c("thistle1","#960018"))(nrow(color_df_2))
  color_df_2$col = colorRampPalette(rev(c("#084808","#E5FFB2")))(nrow(color_df_2))
  
  
  color_df_3 = data.frame("variable" = unique(df[str_detect(df[[reference_colum]], "Viruses"), reference_colum ] ) )
  color_df_3$col = colorRampPalette(c("#FFFFAA","orange","#A34700"))(nrow(color_df_3))
  
  
  other_df <- data.frame(variable=c('Other Bacteria', 'Other Viruses', "Other Eukaryota","Archaea"),
                         col=c("darkblue", "#5A2E00","#032A05", "#FEEFFA"))
  
  color_df = rbind(color_df_1,color_df_2,color_df_3,other_df)
  
  color_dict = color_df$col
  names(color_dict) = color_df$variable
  return(color_dict) }

import_ReadsClassif=function(sra , path_results, data_name, column_label, return_all=FALSE){ 
  # 1: Import data // minor differences between kraken and salmon
  if(data_name=="kraken"){ 
    rclassif_path=paste0(path_results, sra, "/ReadsClassif.txt") 
    x=read.table(rclassif_path, header=FALSE, fill = TRUE, sep = "\t",quote = "")
    colnames(x)=c("read", "taxid", "superkingdom","phylum", "class", "order", "family", "genus", "species")
    
    x[x$taxid %in% c(0, 1,198431,155900,131567,28384), ranks] = "unclassified"
    x$taxid=NULL
    
  } else if(data_name=="hybrid"){
    rclassif_path=paste0(path_results, sra, "/ReadsClassif.txt") 
    x=read.table(rclassif_path, header=FALSE, fill = TRUE, sep = "\t",quote = "")
    colnames(x)=c("read", "superkingdom","phylum", "class", "order", "family", "genus", "species")
    ### Some errors where one row is compose of twice the SRA id
    x=x[x$superkingdom != "",]
    x[x$read == x$superkingdom , ranks] = "Unclassified"
    
  } else{
    rclassif_path=paste0(path_results, sra, "/ContigsToReads/ReadsClassif.txt") 
    x=read.table(rclassif_path, header=TRUE, fill = TRUE, sep = "\t", quote = "")
    colnames(x)=c("read", "superkingdom","phylum", "class", "order", "family", "genus", "species", "strain")
    x$strain = NULL
  }
  
  # Change classif names (except for viruses: holes in the classification)
  x_vir = x[x$superkingdom == "Viruses",]
  x = x[! x$superkingdom == "Viruses",]
  
  # ex: Replace "unclassified Victivallaceae genus" by "unclassified" and "0" by "unclassified"
  # 0 means that the LCA was used by my own algorithm for blast classification
  x[ranks] <- lapply(x[ranks], function(col) {
    gsub("^unclassified.*", "unclassified", col) })
  
  x[ranks] <- lapply(x[ranks], function(col) {
    gsub("0", "unclassified", col) })
  
  x=rbind(x, x_vir)

  x[x==""]="unclassified"
  
  if (return_all==TRUE){ return(x) }
  
  x$label=ifelse(x$class == "Mammalia" , "Homo sapiens", paste0(x$superkingdom, " - ", x$phylum))
  x$label=ifelse(x$superkingdom == "unclassified" , "Unclassified", x$label)
  x$label=ifelse(x$superkingdom == "unassembled" , "Unassembled", x$label)
  
  
  x[,c("taxid", unlist(ranks))]=NULL
  colnames(x)[2]=column_label
  return(x)
}

################################################################################
# V A R I A B L E S
################################################################################
#path = "/run/user/1001/gvfs/sftp:host=core.cluster.france-bioinformatique.fr,user=acolajanni/shared/projects/microbiome_translocation/"
path = "/shared/projects/microbiome_translocation/"

SRA_id = args[1]
dataset = args[2]


# SRA_id="SRR14418894"
SRA_id="all"
# #SRA_id="A21"
# dataset="Simulation"
dataset="Douek_Cleveland"
# SRA_id="responder"
# SRA_id="SRR14418854"
# SRA_id="all"
print(dataset)
print(SRA_id)

path_result=paste0(path,"results/",dataset,"/")
path_data=paste0(path,"data/",dataset,"/")

sra_list=readLines(paste0(path_data,"/sra_list_RNA.txt"))
sra_metadata=read.table(paste0(path_data,"/sra_metadata.tsv"),header=TRUE)
sra_metadata_class=read.table(paste0(path_data,"/Metadata.csv"),header=TRUE, sep=",")


patient_id=sra_metadata[sra_metadata$Run==SRA_id,]$Title

load(paste0(path_result,"/Taxon_correspondance.Rdata"))
load(paste0(path_result,"/color_dict.Rdata"))
##########

result_dir_ku=paste0(path_result,"/kraken/microbialDB/kuniq/")
result_dir_k2=paste0(path_result,"/kraken/nt/k2uniq/")
result_dir_salmon=paste0(path_result,"/Contigs/")
result_dir_blastKu=paste0(path_result,"hybrid/hybrid_Blast-kuniq/")
result_dir_kuk2=paste0(path_result,"hybrid/hybrid_kuniq-kraken2/")


### Import classifs
if (SRA_id %in% c("all","responder","non_responder") ){
  
  if (SRA_id == "responder"){
    sra_list = sra_metadata_class[sra_metadata_class$subject_classification == "Responder (>3yr on ART)",]$Run
  }else if (SRA_id == "non_responder"){
    sra_list = sra_metadata_class[sra_metadata_class$subject_classification != "Responder (>3yr on ART)",]$Run
  }
  
  classifs_list = parallel::mclapply(sra_list, function(sra){
    return(list(
      #"KU"      = import_ReadsClassif(sra, result_dir_ku, "kraken", "KrakenUniq") , 
      "blast"   = import_ReadsClassif(sra, result_dir_salmon, "salmon", "Blast") ,
      "blastKU" = import_ReadsClassif(sra, result_dir_blastKu, "hybrid", "hybrid_Blast_kuniq") ,
      "KuK2"    = import_ReadsClassif(sra, result_dir_kuk2, "hybrid", "KuK2") 
      #, "kraken2" = import_ReadsClassif(sra, result_dir_k2, "kraken", "Kraken2")
      ))
  }, mc.cores=20)
  
  classifs = list()
  for (method in names(classifs_list[[1]]) ) {
    classifs[[method]] <- do.call(rbind, lapply(classifs_list, function(x) x[[method]]))
  }
  rm(classifs_list)

  
} else {
  classifs= list(
    #"KU"      = import_ReadsClassif(SRA_id, result_dir_ku, "kraken", "KrakenUniq") , 
    "blast"   = import_ReadsClassif(SRA_id, result_dir_salmon, "salmon", "Blast") ,
    "blastKU" = import_ReadsClassif(SRA_id, result_dir_blastKu, "hybrid", "hybrid_Blast_kuniq") ,
    "KuK2"    = import_ReadsClassif(SRA_id, result_dir_kuk2, "hybrid", "KuK2") 
    #,"kraken2" = import_ReadsClassif(SRA_id, result_dir_k2, "kraken", "Kraken2") 
    )
  
}

merged_df <- Reduce(function(x, y) merge(x, y, by = "read", all = TRUE), classifs)
rm(classifs)


### Merge rare taxons
taxons=as.data.frame(table(unlist(merged_df[ , -1])))
colnames(taxons)[1] = "taxon_full_name"

taxons$superkingdom=sapply(as.character(taxons$taxon_full_name), function(x){
  x=strsplit(x, " - ", fixed=TRUE )[[1]][1]
  return(x)
})
taxons$phylum=sapply(as.character(taxons$taxon_full_name), function(x){
  x=strsplit(x, " - ",fixed=TRUE)[[1]][2]
  return(x)
})


### Creating color palette
rename_taxa=do.call(rbind, list(
  "bacteria"=aggregate_taxa_count_df(taxons,5,"taxon_full_name","Freq", "Bacteria"),
  "viruses"=aggregate_taxa_count_df(taxons,4,"taxon_full_name","Freq", "Viruses"),
  "Archaea"=aggregate_taxa_count_df(taxons,1,"taxon_full_name","Freq", "Archaea"),
  "Euk"=aggregate_taxa_count_df(taxons,4,"taxon_full_name","Freq", "Eukaryota") ) )

as.character(taxons$taxon_full_name)[! as.character(taxons$taxon_full_name) %in% rename_taxa$taxon_full_name]


taxons[,c("superkingdom","phylum")]=NULL
taxons=taxons[! taxons$taxon_full_name %in% rename_taxa$taxon_full_name , ]
taxons$final_name=taxons$taxon_full_name
rename_taxa=rbind(taxons, rename_taxa)
row.names(rename_taxa) = NULL

rename_taxa[ , ]$final_name = 
rename_taxa$final_name = ifelse(str_detect(rename_taxa$taxon_full_name, "myco"), "Eukaryota - Fungi" , as.character(rename_taxa$final_name) )

if ( "Unclassified - Unclassified" %in% rename_taxa$taxon_full_name){
rename_taxa[rename_taxa$taxon_full_name == "Unclassified - Unclassified",]$final_name = "Unclassified"
}

rename_taxa[str_detect(rename_taxa$taxon_full_name, "Archaea"),]$final_name = "Archaea"


# Replace values in selected columns using dictionary
cols_to_replace <- c("KrakenUniq", "Blast", "hybrid_Blast_kuniq","KuK2","Kraken2")  # List of columns to modify
cols_to_replace <- c("Blast", "hybrid_Blast_kuniq","KuK2")  # List of columns to modify

# Loop through selected columns and replace values
for (col in cols_to_replace) {
  print(col)
  merged_df[[col]]=rename_taxa$final_name[match(merged_df[[col]], as.character(rename_taxa$taxon_full_name))]
}

color_dict=Taxonomy_small_color_palette(rename_taxa,"final_name")
color_dict[["Unassembled"]]="#444444"
color_dict[["Homo sapiens"]]="darkorchid"
color_dict[["Unclassified"]]="#999999"

tax=names(color_dict)

tax_order=c(
  "Homo sapiens","Unassembled","Unclassified",
  tax[str_detect(tax,"Bacteria") & ! str_detect(tax, "Other")], "Other Bacteria","Archaea",
  tax[str_detect(tax,"Eukaryota") & ! str_detect(tax, "Other") & ! str_detect(tax,"myco")], "Other Eukaryota",
  tax[str_detect(tax,"Viruses") & ! str_detect(tax, "Other")], "Other Viruses")



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# S A N K E Y   P L O T S
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Building dataframe for Sankey
df <- merged_df %>%
  #make_long(Blast, hybrid_Blast_kuniq, KuK2, KrakenUniq)
  make_long(Blast, hybrid_Blast_kuniq, KuK2)

p1=ggplot(df, aes(x = x, 
                  next_x = next_x, 
                  node = node, 
                  next_node = next_node,
                  fill = factor(node) )) + 
  geom_sankey(flow.alpha=0.85,node.color=1,space=150) +
  scale_fill_manual(values = color_dict) +
  #xlab(paste("Patient",patient_id))+  
  #scale_x_discrete(labels = c("Blue pipeline","Red pipeline")) +
  theme_sankey(base_size = 11 ) +
  theme(legend.position = "bottom")
p1


ggsave(p1,filename=paste0("/shared/projects/microbiome_translocation/results/Douek_Cleveland/fig_Benchmark/Sankey_all_",SRA_id,".jpeg"),
       device = "jpeg",width = 30, height = 26,dpi = 350 , units = "cm", create.dir = TRUE)


# Remove rows where both hybrid/combined predicts human
d2 <- merged_df[! (#merged_df$KrakenUniq == "Homo sapiens" &  
                     merged_df$hybrid_Blast_kuniq == "Homo sapiens" & 
                     merged_df$KuK2 == "Homo sapiens"),]%>%
  #make_long(Blast, hybrid_Blast_kuniq, KuK2, KrakenUniq)
  make_long(Blast, hybrid_Blast_kuniq, KuK2)

p2=ggplot(d2, aes(x = x, 
                  next_x = next_x, 
                  node = node, 
                  next_node = next_node,
                  fill = factor(node) )) + 
  geom_sankey(flow.alpha=0.85,node.color=1,space=150) +
  scale_fill_manual(values = color_dict) +
  #xlab(paste("Patient",patient_id))+  
  #scale_x_discrete(labels = c("Blue pipeline","Red pipeline")) +
  theme_sankey(base_size = 11 ) +
  theme(legend.position = "bottom")
p2
ggsave(p2,filename=paste0("/shared/projects/microbiome_translocation/results/Douek_Cleveland/fig_Benchmark/Sankey_nohuman_",SRA_id,".jpeg"),
       device = "jpeg",width = 40, height = 26,dpi = 350 , units = "cm", create.dir = TRUE)


# Remove rows where both hybrid/combined predicts human and unclassified
d3 <- merged_df[! (#merged_df$KrakenUniq %in% c("Homo sapiens","Unclassified") &  
                     merged_df$hybrid_Blast_kuniq %in% c("Homo sapiens","Unclassified")& 
                     merged_df$KuK2  %in% c("Homo sapiens","Unclassified")  ),]%>%
  #make_long(Blast, hybrid_Blast_kuniq, KuK2, KrakenUniq)
  make_long(Blast, hybrid_Blast_kuniq, KuK2)

p3=ggplot(d3, aes(x = x, 
                  next_x = next_x, 
                  node = node, 
                  next_node = next_node,
                  fill = factor(node) )) + 
  geom_sankey(flow.alpha=0.85,node.color=1,space=150) +
  scale_fill_manual(values = color_dict) +
  #xlab(paste("Patient",patient_id))+  
  #scale_x_discrete(labels = c("Blue pipeline","Red pipeline")) +
  theme_sankey(base_size = 11 ) +
  theme(legend.position = "bottom")
p3
ggsave(p3,filename=paste0("/shared/projects/microbiome_translocation/results/Douek_Cleveland/fig_Benchmark/Sankey_discordant_",SRA_id,".jpeg"),
       device = "jpeg",width = 30, height = 26,dpi = 350 , units = "cm", create.dir = TRUE)



d4 <- merged_df[! (#merged_df$KrakenUniq %in% c("Homo sapiens","Unclassified") &  
                     merged_df$hybrid_Blast_kuniq %in% c("Homo sapiens","Unclassified")& 
                     merged_df$KuK2 %in% c("Homo sapiens","Unclassified")  ) & 
                  !(merged_df$KuK2 == "Viruses - Uroviricota" & merged_df$hybrid_Blast_kuniq == "Homo sapiens")
                ,]%>%
  #make_long(Blast, hybrid_Blast_kuniq, KuK2, KrakenUniq)
  make_long(Blast, hybrid_Blast_kuniq, KuK2)

p4=ggplot(d4, aes(x = x, 
                  next_x = next_x, 
                  node = node, 
                  next_node = next_node,
                  fill = factor(node) )) + 
  geom_sankey(flow.alpha=0.85,node.color=1,space=150) +
  scale_fill_manual(values = color_dict) +
  #xlab(paste("Patient",patient_id))+  
  #scale_x_discrete(labels = c("Blue pipeline","Red pipeline")) +
  theme_sankey(base_size = 11 ) +
  theme(legend.position = "bottom")
p4
ggsave(p4,filename=paste0("/shared/projects/microbiome_translocation/results/Douek_Cleveland/fig_Benchmark/Sankey_remaining_",SRA_id,".jpeg"),
       device = "jpeg",width = 30, height = 26,dpi = 350 , units = "cm", create.dir = TRUE)


read_set=unique( c(
  merged_df[! merged_df$hybrid_Blast_kuniq %in% c("Unclassified", "Unassembled","Homo sapiens") ,]$read,
  merged_df[! merged_df$KuK2 %in% c("Unclassified", "Unassembled","Homo sapiens") ,]$read,
  merged_df[! merged_df$hybrid_Blast_kuniq %in% c("Unclassified", "Unassembled","Homo sapiens") ,]$read ))


d5 = merged_df[merged_df$read %in% read_set ,] %>% make_long(Blast, hybrid_Blast_kuniq, KuK2)

rm(list = c("df","d2","d3","d4"))
save.image(file=paste0('/shared/projects/microbiome_translocation/results/Douek_Cleveland/fig_Benchmark/',SRA_id,"_Sankey_env.Rdata"))


p5=ggplot(d5, aes(x = x, 
                  next_x = next_x, 
                  node = node, 
                  next_node = next_node,
                  fill = factor(node) )) + 
  geom_sankey(flow.alpha=0.85,node.color=1,space=150) +
  scale_fill_manual(values = color_dict) +
  #xlab(paste("Patient",patient_id))+  
  #scale_x_discrete(labels = c("Blue pipeline","Red pipeline")) +
  theme_sankey(base_size = 11 ) +
  theme(legend.position = "bottom")
p5
ggsave(p5,filename=paste0("/shared/projects/microbiome_translocation/results/Douek_Cleveland/fig_Benchmark/Sankey_union_classified_",SRA_id,".jpeg"),
       device = "jpeg",width = 40, height = 26,dpi = 350 , units = "cm", create.dir = TRUE)




rm(list = c("df","d2","d3","d4","d5"))

if (SRA_id != "all"){ 
  stop() 
}


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# B A R P L O T S   -   S E T   O F   Super Kingdoms
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#methods=colnames(merged_df[2:ncol(merged_df)])
labels2lines <- c(
  "Blast" = "Contig-based \n Blastn - nt database",
  "hybrid_Blast_kuniq" = "Hybrid \n Blast + KrakenUniq",
  "KuK2" = "Kmer-based \n KrakenUniq + Kraken2")

labels2lines <- c(
  "Blast" = "- Contig-based -\n(nt database)",
  "hybrid_Blast_kuniq" = "- Hybrid -\nContig-based (nt database) +\nKrakenUniq (microbialDB)",
  "KuK2" = "- K-mer-based Combined -\nKrakenUniq (microbialDB) +\nKraken2 (nt database)"
)

label_df=as.data.frame(labels2lines)
label_df$method=rownames(label_df)

methods=c("Blast","hybrid_Blast_kuniq", "KuK2")

merged_df$vir  = apply(merged_df[, methods], 1, function(row) any(str_detect(row, "Viruses")))
merged_df$euk  = apply(merged_df[, methods], 1, function(row) any(str_detect(row, "Eukaryota")))
merged_df$bact = apply(merged_df[, methods], 1, function(row) any(str_detect(row, "Bacteria")))
merged_df$.id <- str_extract(merged_df$read, "^[^.]+")

total=as.data.frame(table(merged_df$.id))
colnames(total)=c(".id","total")

classif_df <- merged_df %>%
  pivot_longer(cols = all_of(methods), names_to = "method", values_to = "taxon_full_name") %>%
  group_by(.id, vir, euk, bact, method, taxon_full_name) %>%
  summarise(count = n(), .groups = "drop")

n=unique(classif_df$taxon_full_name)
setdiff(n,tax_order)



classif_df=merge(classif_df,label_df,by="method")
classif_df=merge(total,classif_df,by=".id")
classif_df$prop=classif_df$count / classif_df$total

tax_order=c(
  "Unassembled","Unclassified","Homo sapiens",
  tax[str_detect(tax,"Bacteria") & ! str_detect(tax, "Other")], "Other Bacteria","Archaea",
  tax[str_detect(tax,"Eukaryota") & ! str_detect(tax, "Other") & ! str_detect(tax,"myco")], "Other Eukaryota",
  tax[str_detect(tax,"Viruses") & ! str_detect(tax, "Other")], "Other Viruses")

classif_df$taxon_full_name = factor(classif_df$taxon_full_name, levels=tax_order)

legend_plot=ggplot(classif_df, aes(x=count,y=.id,fill=taxon_full_name))+
  #geom_bar(position="stack", stat="identity",color="black",width=0.95) +
  geom_bar(position="stack", stat="identity", width = 0.94) + 
  theme_linedraw()+ xlab("Read quantity") +
  scale_fill_manual(values=color_dict, name="")+
  guides(fill = guide_legend(ncol = 1)) +
  theme(
    legend.position = "right",
    legend.key = element_rect(color="black"),
    legend.text=element_text(size=12),
    legend.title=element_text(size=12),
    legend.key.size =unit(.6,"cm") )

legend=cowplot::get_plot_component(legend_plot,'guide-box-right',return_all = TRUE)


total_per_abrev=aggregate()



classif_df$taxon_full_name <- factor(classif_df$taxon_full_name, levels=tax_order )

### Computing total for each superkingdom (change order of row)
total_vir <- classif_df[classif_df$method=="Blast" & classif_df$vir,] %>%
  group_by(.id) %>%summarise(total_vir = sum(count), .groups = "drop")

total_bact <- classif_df[classif_df$method=="Blast" & classif_df$bact,] %>%
  group_by(.id) %>% summarise(total_bact = sum(count), .groups = "drop")

total_euk <- classif_df[classif_df$method=="Blast" & classif_df$euk,] %>%
  group_by(.id) %>% summarise(total_euk = sum(count), .groups = "drop")

total_per_sk = merge(total_euk, 
                     merge(total_vir, total_bact, by=".id")
                     , by=".id")

classif_df=merge(classif_df, total_per_sk, by=".id")

plots <- lapply(c("Bacteria","Eukaryota","Viruses"), function(sk) {

  if      (sk == "Bacteria"){abrev="bact"}
  else if (sk == "Viruses"){abrev="vir"}
  else if (sk == "Eukaryota"){abrev="euk"}
  
  col_total=paste0("total_",abrev)
  
  classif_df$current_total = classif_df[[col_total]]
  classif_df$superkingdom = sk
  
  barplot=ggplot(classif_df[classif_df[[abrev]],], 
         aes(x=count,y=reorder(.id, current_total),fill=taxon_full_name))+
    #geom_bar(position="stack", stat="identity",color="black",width=0.95) +
    geom_bar(position="stack", stat="identity", width = 0.93) + 
    theme_linedraw()+ xlab("Read quantity") +
    scale_fill_manual(values=color_dict, name="")+
    theme(plot.title = element_text(hjust = 0.5, size=16),
          legend.position = "right",
          legend.key = element_rect(color="black"),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          legend.key.size =unit(.6,"cm"),
          axis.text.x=element_text(angle=0,size=8),
          axis.text.y=element_text(angle=0,size=4),
          strip.background = element_rect(fill="#333333", color="black") , 
          strip.text = element_text(size = 12, face = "bold",color="white"),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill = "#CCCCCC"),    
          panel.grid.major = element_line(colour = "#444444",),
          panel.grid.minor = element_line(colour = "#444444",)) +
    facet_grid(vars(superkingdom) , vars(labels2lines), scales="free_x" )

  if (sk!="Bacteria"){ 
  barplot = barplot + theme(
    strip.text.x = element_blank() , 
    strip.background.x = element_blank(),
    plot.margin = unit( c(0,0,0,0) , units = "lines" ) ) }
  
  if (sk!="Eukaryota"){ 
    barplot = barplot + theme(
      axis.title.x = element_blank() ) }
  
  return(barplot)
  
})

# Combine plots and harmonize legend
final_plot <- wrap_plots(plots, ncol = 1) + 
  plot_layout(guides = "collect") & theme(legend.position = "none")


grid_plot=plot_grid(final_plot,legend , rel_widths = c(10,2.5) )

print(grid_plot)
ggsave(grid_plot,filename=paste0("/shared/projects/microbiome_translocation/results/Douek_Cleveland/fig_Benchmark/barplot_shared_",SRA_id,".jpeg"),
       device = "jpeg",width = 40, height = 26,dpi = 350 , units = "cm", create.dir = TRUE)


save.image(file=paste0('/shared/projects/microbiome_translocation/results/Douek_Cleveland/fig_Benchmark/',SRA_id,"_env.Rdata"))
#load(paste0('/shared/projects/microbiome_translocation/results/Douek_Cleveland/fig_Benchmark/',SRA_id,"_env.Rdata"))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# B A R P L O T S   G L O B A L
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
all_methods=colnames(merged_df)[2:6]

labels2lines <- c(
  "Blast" = "- Contig-based -\n(nt database)",
  "hybrid_Blast_kuniq" = "- Hybrid -\nContig-based (nt database) +\nKrakenUniq (microbialDB)",
  "KuK2" = "- K-mer-based Combined -\nKrakenUniq (microbialDB) +\nKraken2 (nt database)")
label_df=as.data.frame(labels2lines)
label_df$method=rownames(label_df)

total=as.data.frame(table(merged_df$.id))
colnames(total)=c(".id","total")

quant_df <- merged_df %>%
  pivot_longer(cols = all_of(all_methods),
               names_to = "method", 
               values_to = "taxon_full_name") %>%
  group_by(.id, method, taxon_full_name) %>%
  summarise(count = n(), .groups = "drop")

quant_df=merge(total,quant_df,by=".id")
quant_df$prop=quant_df$count / quant_df$total
quant_df=merge(quant_df,label_df, by="method")



tax_order=c(
  "Unassembled","Unclassified","Homo sapiens",
  tax[str_detect(tax,"Bacteria") & ! str_detect(tax, "Other")], "Other Bacteria","Archaea",
  tax[str_detect(tax,"Eukaryota") & ! str_detect(tax, "Other") & ! str_detect(tax,"myco")], "Other Eukaryota",
  tax[str_detect(tax,"Viruses") & ! str_detect(tax, "Other")], "Other Viruses")

quant_df$taxon_full_name = factor(quant_df$taxon_full_name, levels=tax_order)
  
### Proportion plots
quantity_plot=ggplot(quant_df, aes(x=count,y=reorder(.id,total),fill=taxon_full_name))+
  geom_bar(position="stack", stat="identity", width = 0.95) + 
  theme_linedraw()+ xlab("Read quantity") +
  scale_fill_manual(values=color_dict, name="")+
  guides(fill = guide_legend(ncol = 1)) +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        legend.position = "right",
        legend.key = element_rect(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.key.size =unit(.6,"cm"),
        axis.text.x=element_text(angle=0,size=8),
        strip.background = element_rect(fill="#333333", color="black") , 
        strip.text = element_text(size = 12, face = "bold",color="white"),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "#CCCCCC"),    
        panel.grid.major = element_line(colour = "#444444",),
        panel.grid.minor = element_line(colour = "#444444",)) +
  facet_wrap(vars(labels2lines))
quantity_plot

prop_plot=ggplot(quant_df, aes(x=prop,y=reorder(.id,total),fill=taxon_full_name))+
  geom_bar(position="stack", stat="identity", width = 0.95) + 
  theme_linedraw()+ xlab("Read proportion") +
  scale_fill_manual(values=color_dict, name="")+
  guides(fill = guide_legend(ncol = 1)) +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        legend.position = "right",
        legend.key = element_rect(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.key.size =unit(.6,"cm"),
        axis.text.x=element_text(angle=0,size=8),
        strip.background = element_rect(fill="#333333", color="black") , 
        strip.text = element_text(size = 12, face = "bold",color="white"),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "#333333"),    
        panel.grid.major = element_line(colour = "#444444",),
        panel.grid.minor = element_line(colour = "#444444",)) +
  facet_wrap(vars(labels2lines))
prop_plot

quant_df_smaller=quant_df[! quant_df$taxon_full_name %in% c("Unassembled","Unclassified","Homo sapiens"),]

prop_plot_non_human=ggplot(quant_df_smaller, aes(x=count,y=reorder(.id,total),fill=taxon_full_name))+
  geom_bar(position="fill", stat="identity", width = 0.94) + 
  theme_linedraw()+ xlab("Read proportion") +
  scale_fill_manual(values=color_dict, name="")+
  guides(fill = guide_legend(ncol = 1)) +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        legend.position = "right",
        legend.key = element_rect(color="black"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.key.size =unit(.6,"cm"),
        axis.text.x=element_text(angle=0,size=8),
        strip.background = element_rect(fill="#333333", color="black") , 
        strip.text = element_text(size = 12, face = "bold",color="white"),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "#333333"),    
        panel.grid.major = element_line(colour = "#444444",),
        panel.grid.minor = element_line(colour = "#444444",)) +
  facet_wrap(vars(labels2lines))

ggsave(quantity_plot,filename=paste0("/shared/projects/microbiome_translocation/results/Douek_Cleveland/fig_Benchmark/barplot_quantity_shared_",SRA_id,".jpeg"),
       device = "jpeg",width = 40, height = 26,dpi = 350 , units = "cm", create.dir = TRUE)
ggsave(prop_plot,filename=paste0("/shared/projects/microbiome_translocation/results/Douek_Cleveland/fig_Benchmark/barplot_prop_shared_",SRA_id,".jpeg"),
       device = "jpeg",width = 40, height = 26,dpi = 350 , units = "cm", create.dir = TRUE)
ggsave(prop_plot_non_human,filename=paste0("/shared/projects/microbiome_translocation/results/Douek_Cleveland/fig_Benchmark/barplot_prop_noHUman_shared_",SRA_id,".jpeg"),
       device = "jpeg",width = 40, height = 26,dpi = 350 , units = "cm", create.dir = TRUE)



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# F I N D   A S S O C I A T I O N S   B E T W E E N   S P (uroviricota-human // Ascomycota-Bacillota)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


merged_df$.id <- str_extract(merged_df$read, "^[^.]+")
total=as.data.frame(table(merged_df$.id))
colnames(total)=c(".id","total")
merged_df=merge(merged_df, total, by=".id")

sk_df=aggregate(read ~ .id + Blast + total, merged_df, length)
sk_df$prop=sk_df$read / sk_df$total
sk_df_u=sk_df[sk_df$Blast %in% c("Homo sapiens", "Unassembled","Unclassified"),]

mean(sk_df_u[sk_df_u$Blast == "Unclassified" , ]$prop )
mean(sk_df_u[sk_df_u$Blast == "Unassembled" , ]$prop )


### Inspect certain association:
FungiBact_df=merged_df[str_detect(merged_df$hybrid_Blast_kuniq,"myco") & str_detect(merged_df$KuK2,"Bacteria") , c("read","hybrid_Blast_kuniq","KuK2")]
UroHuman_df=merged_df[merged_df$hybrid_Blast_kuniq == "Homo sapiens" & merged_df$KuK2 == "Viruses - Uroviricota" , c("read","hybrid_Blast_kuniq","KuK2") ]


freqFungiBact=aggregate(read ~ hybrid_Blast_kuniq + KuK2, FungiBact_df, length)
FungiBact_df=FungiBact_df[FungiBact_df$hybrid_Blast_kuniq == "Eukaryota - Ascomycota" & FungiBact_df$KuK2 == "Bacteria - Bacillota",]
writeLines(FungiBact_df$read, "/shared/projects/microbiome_translocation/results/Douek_Cleveland/FungiBact_reads.txt" ) 
writeLines(UroHuman_df$read, "/shared/projects/microbiome_translocation/results/Douek_Cleveland/HumanUroviricota_reads.txt" ) 

reads_list=unique(c(FungiBact_df$read, UroHuman_df$read))


classifs_full = parallel::mclapply(sra_list, function(sra){
  bku = import_ReadsClassif(sra, result_dir_blastKu, "hybrid", "hybrid_Blast_kuniq",TRUE) 
  kuk2    = import_ReadsClassif(sra, result_dir_kuk2, "hybrid", "KuK2",TRUE) 
  blast  = import_ReadsClassif(sra, result_dir_salmon, "salmon", "Blast",TRUE) 
  
  colnames(kuk2)[2:8] = paste0(ranks,"_KuK2")
  colnames(bku)[2:8] = paste0(ranks,"_hybrid")
  colnames(blast)[2:8] = paste0(ranks,"_blast")
  
  merged_data <- bku %>%
    full_join(kuk2, by = "read") %>%
    full_join(blast, by = "read")
  
  aggregated_data <- aggregate(read ~ ., data = merged_data, FUN = length)
  return(aggregated_data)
  
},  mc.cores=20)

names(classifs_full) = sra_list

################################################################################
blast = parallel::mclapply(sra_list, function(sra){return(classifs_full[[sra]]$blast)},mc.cores=20)
names(blast) = sra_list

blast = parallel::mclapply(sra_list, function(sra){
  x=blast[[sra]]
  x=aggregate(read ~ superkingdom_blast+phylum_blast+class_blast+order_blast+family_blast+genus_blast+species_blast, 
              x, length)
  return(x)},mc.cores=20)
names(blast) = sra_list

blast=ldply(blast)
total=aggregate(read ~ .id, blast, sum)
colnames(total)[2]='total'
blast=merge(blast,total, by=".id")
blast$prop=blast$read / blast$total

blast_Human=blast[ !(!blast$superkingdom_blast %in% c("unclassified","unassembled") & blast$class_blast != "Mammalia"),]


blast_Human=aggregate(read ~ .id + total, blast_Human, sum)
blast_Human$prop=blast_Human$read / blast_Human$total
blast_Human$prop_non_human= 1 - blast_Human$prop

blast_Human$read_non_human=blast_Human$total - blast_Human$read
mean(blast_Human$read_non_human)
mean(blast_Human$prop_non_human)

blast_noHuman=blast[ (!blast$superkingdom_blast %in% c("unclassified","unassembled") & blast$class_blast != "Mammalia") ,]


################################################################################

classifs_full2 = parallel::mclapply(sra_list, function(sra){
  tmp_list=classifs_full[[sra]]
  return(
    merge(tmp_list[[1]],
          merge(tmp_list[[2]],tmp_list[[3]],by="read",all=TRUE ),
          by="read",all=TRUE ) 
    )
  },  mc.cores=20)
names(classifs_full2) = sra_list





classifs_full2=ldply(classifs_full2)

fungi_bact=classifs_full2[classifs_full2$read %in% FungiBact_df$read , ]
uro_Human=classifs_full2[classifs_full2$read %in% UroHuman_df$read , ]


uro=aggregate(read ~ .id + species_KuK2, uro_Human, length)
  
uro_global=aggregate(read ~ species_KuK2, uro, sum)
uro_global$total = nrow(UroHuman_df)
uro_global$prop=uro_global$read/uro_global$total

fungiBact=aggregate(read ~ species_KuK2 + species_hybrid, fungi_bact, length)

fungi=aggregate(read ~  order_hybrid, fungi_bact, length)


### Retrouver le total 

path = "/shared/projects/microbiome_translocation/"




# SRA_id="SRR14418872"
SRA_id="all"
# #SRA_id="A21"
# dataset="Simulation"
dataset="Douek_Cleveland"
# SRA_id="SRR14418854"
# SRA_id="all"
print(dataset)
print(SRA_id)


path_result=paste0(path,"results/",dataset,"/")
path_data=paste0(path,"data/",dataset,"/")

sra_list=readLines(paste0(path_data,"/sra_list_RNA.txt"))
sra_metadata=read.table(paste0(path_data,"/sra_metadata.tsv"),header=TRUE)


patient_id=sra_metadata[sra_metadata$Run==SRA_id,]$Title

load(paste0(path_result,"/Taxon_correspondance.Rdata"))
load(paste0(path_result,"/color_dict.Rdata"))
##########

result_dir_ku=paste0(path_result,"/kraken/microbialDB/kuniq/")
result_dir_k2=paste0(path_result,"/kraken/nt/k2uniq/")
result_dir_salmon=paste0(path_result,"/Contigs/")
result_dir_blastKu=paste0(path_result,"hybrid/hybrid_Blast-kuniq/")
result_dir_kuk2=paste0(path_result,"hybrid/hybrid_kuniq-kraken2/")

### Import classifs

#classifs_full2 = parallel::mclapply(sra_list[c(2,8,10,11,14,16,20,26,28,30,32,34,38,44,46,47,50,52)], function(sra){
classifs_full= parallel::mclapply(sra_list, function(sra){
  bku = import_ReadsClassif(sra, result_dir_blastKu, "hybrid", "hybrid_Blast_kuniq",TRUE) 
  bku$bku=paste0(bku$superkingdom,"|",bku$family,"|",bku$species)
  bku[,ranks]=NULL
  
  kuk2    = import_ReadsClassif(sra, result_dir_kuk2, "hybrid", "KuK2",TRUE) 
  kuk2$kuk2=paste0(kuk2$superkingdom,"|",kuk2$family,"|",kuk2$species)
  kuk2[,ranks]=NULL
  
  ku    = import_ReadsClassif(sra, result_dir_ku, "kraken", "kuniq",TRUE) 
  ku$ku=paste0(ku$superkingdom,"|",ku$family,"|",ku$species)
  ku[,ranks]=NULL
  
  
  k2    = import_ReadsClassif(sra, result_dir_k2, "kraken", "kraken2",TRUE) 
  k2$k2=paste0(k2$superkingdom,"|",k2$family,"|",k2$species)
  k2[,ranks]=NULL
  
  blast  = import_ReadsClassif(sra, result_dir_salmon, "salmon", "Blast",TRUE) 
  blast$blast=paste0(blast$superkingdom,"|",blast$family,"|",blast$species)
  blast[,ranks]=NULL
  
  merged_data <- bku %>%
    full_join(k2, by = "read") %>%
    full_join(ku, by = "read") %>%
    full_join(kuk2, by = "read") %>%
    full_join(blast, by = "read")
  
  #aggregated_data <- aggregate(read ~ ., data = merged_data, FUN = length)
  return(merged_data)
  
},  mc.cores=11)

names(classifs_full) = sra_list
merged_df=ldply(classifs_full)


#merged_df <- Reduce(function(x, y) merge(x, y, by = "read", all = TRUE), classifs_full)
rm(classifs_full)

# merged_df$.id <- str_extract(merged_df$read, "^[^.]+")

#save(merged_df, file="/shared/projects/microbiome_translocation/results/Douek_Cleveland/fig_Benchmark/merged_classifs_sankey.RData")
load("/shared/projects/microbiome_translocation/results/Douek_Cleveland/fig_Benchmark/merged_classifs_sankey.RData")

merged_agg=aggregate(read ~ .,merged_df, length)
total=aggregate(read ~ .id, merged_agg, sum)
colnames(total)[2]="total"
merged_agg2=merge(merged_agg, total, by=".id")

total_trimmed=read.delim(file="/shared/projects/microbiome_translocation/results/Douek_Cleveland/trimmed_reads_stats.tsv")[,c(1,4)]
total_trimmed=total_trimmed[!str_detect(total_trimmed$file, "2.trimmed.fastq.gz") , ]
total_trimmed$file=basename(total_trimmed$file)
total_trimmed$.id = sapply(total_trimmed$file, function(x) str_split(x, "_")[[1]][1] )
total_trimmed$file=NULL

####
total_raw=read.delim(file="/shared/projects/microbiome_translocation/results/Douek_Cleveland/raw_reads_stats.tsv")[,c(1,4)]
total_raw=total_raw[!str_detect(total_raw$file, "_2.fastq.gz") , ]
total_raw$file=basename(total_raw$file)
total_raw$.id = sapply(total_raw$file, function(x) str_split(x, "_")[[1]][1] )
total_raw$file=NULL
####


merged_agg3=merge(merged_agg2, total_trimmed, by=".id")

#merged_global=aggregate(read ~ bku + kuk2 + blast + ku + k2, merged_agg3, sum)


get_summary_numbers=function(merged_df, column){
  
  method=merged_df[,c(".id","read","total","num_seqs",column)]
  method=aggregate(read ~ . , method, sum)
  
  method_human= method[method[[column]] == fixed("unclassified|unclassified|unclassified") | 
                         method[[column]] == fixed("unassembled|unassembled|unassembled") | 
                       str_detect(method[[column]],"Hominidae") ,  ]

  method_human = method_human %>%
    mutate(parsed_column = str_extract(method_human[[column]], "[^|]+$"))
  
  method_human$parsed_column = ifelse((!method_human$parsed_column %in% c("unclassified","unassembled")), "Homo sapiens", method_human$parsed_column)
  
  df_wide = method_human %>%
      pivot_wider(names_from = parsed_column, values_from = read, values_fill = 0) %>%
      select(-{{column}}) %>%
      group_by(.id, total, num_seqs) %>%  # Ensure unique rows
      summarise(across(where(is.numeric), sum), .groups = "drop")
  
  if (column != "blast"){ df_wide$unassembled = 0 }
  df_wide$unclassified_all = df_wide$unassembled + df_wide$unclassified
  
  
  df_wide$removed          = df_wide$`Homo sapiens` + df_wide$unassembled + df_wide$unclassified 
  df_wide$removed_prop     = df_wide$removed / df_wide$total
  df_wide$removed_prop_raw = df_wide$removed / df_wide$num_seqs
  
  df_wide$non_human      = df_wide$total - df_wide$removed 
  df_wide$non_human_prop = df_wide$non_human / df_wide$total 
  df_wide$non_human_raw  = df_wide$total / df_wide$num_seqs

  return(df_wide)
}

blast=get_summary_numbers(merged_agg3,"blast")
blast$method="blast"

kuk2=get_summary_numbers(merged_agg3,"kuk2")
kuk2$method="kuk2"

bku=get_summary_numbers(merged_agg3,"bku")
bku$method="bku"

k2=get_summary_numbers(merged_agg3,"k2")
k2$method="k2"

ku=get_summary_numbers(merged_agg3,"ku")
ku$method="ku"

summary_methods=do.call(rbind, list(blast, kuk2, bku, k2, ku))

summary_methods$prop_filtered = summary_methods$total / summary_methods$num_seqs

df_summary <- summary_methods %>%
  group_by(method) %>%
  summarise(across(where(is.numeric), list(mean = mean, sd = sd), .names = "{.col}_{.fn}"))


save(df_summary, file="/shared/projects/microbiome_translocation/results/Douek_Cleveland/fig_Benchmark/numbers_table2.RData")
#load("/shared/projects/microbiome_translocation/results/Douek_Cleveland/fig_Benchmark/numbers_table2.RData")






blast=merged_agg3[,c(".id","blast","read","total","num_seqs")]
blast=aggregate(read ~ . , blast, sum)

blast_human= blast[blast$blast == "unclassified|unclassified|unclassified" | 
                     blast$blast == "unassembled|unassembled|unassembled" | 
                     str_detect(blast$blast,"Hominidae") ,  ]



blast_human=aggregate(read ~ .id + total + num_seqs , blast_human, sum)
blast_human$non_human=blast_human$total - blast_human$read

blast_human$prop     = blast_human$read / blast_human$total
blast_human$prop_raw = blast_human$read / blast_human$num_seqs


blast_human$prop_nonhuman     = blast_human$non_human / blast_human$total
blast_human$prop_nonhuman_raw = blast_human$non_human / blast_human$num_seqs

mean(blast_human$prop_nonhuman_raw)


