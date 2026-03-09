#!/shared/ifbstor1/software/miniconda/envs/r-4.2.3/bin/Rscript

library(stringr)
library(plyr)
library(ggsankey)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(cowplot)
#library(ggstream)

args = commandArgs(trailingOnly=TRUE)

###############################################################################
# F U N C T I O N S
################################################################################
ranks=c("superkingdom","phylum", "class", "order", "family", "genus", "species")

aggregate_taxa_count_df=function(table_df, n=5, column_to_filter, Frequency_column,replace_taxon="Bacteria"){
  table_df=table_df[table_df$superkingdom == replace_taxon , ]
  table_df=table_df[order( table_df[[Frequency_column]],decreasing = TRUE ) , ]
  
  # if (replace_taxon == "Eukaryota"){
  #   
  #   print("eukaryota")
  #   fungi_df=table_df[str_detect(table_df[[column_to_filter]], "Fungi|fungi|myco" ),]
  #   fungi_df$superkingdom = "Eukaryota - Fungi"
  #   
  #   fun=aggregate_taxa_count_df(fungi_df, 2, column_to_filter, Frequency_column, "Eukaryota - Fungi"  )
  #   table_df=table_df[! table_df[[column_to_filter]] %in% fungi_df[[column_to_filter]] , ]    
  #   
  # }
  
  top=head(table_df, n = n)
  low=table_df[!table_df[[column_to_filter]] %in% top[[column_to_filter]] , ]
  
  top$final_name=top[[column_to_filter]]
  low$final_name=sapply(low$superkingdom, function(x){paste0("Other ",replace_taxon) })
  
  table_df=rbind(top,low)
  table_df[,c("superkingdom","phylum")]=NULL
  colnames(table_df)[1]="taxon_full_name"
  
  # if (replace_taxon == "Eukaryota"){
  #   table_df=rbind(fun,table_df)
  # }
  
  
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


Taxonomy_small_color_palette = function(df, reference_colum){
  
  df=df[order(df[[reference_colum]]), ]
  df=df[!str_detect(df[[reference_colum]],"Other"), ]
  
  color_df_1 = data.frame("variable" = unique(df[str_detect(df[[reference_colum]], "Bacteria"), reference_colum  ] ) )
  color_df_1$col = colorRampPalette(c("#e6f3ff","#000099"  ))(nrow(color_df_1))
  
  color_df_2 = data.frame("variable" = unique(df[str_detect(df[[reference_colum]], "Eukaryota") , reference_colum ] ) )
  #color_df_2$col = colorRampPalette(c("thistle1","#960018"))(nrow(color_df_2))
  color_df_2$col = colorRampPalette(rev(c("#084508", "#dfe9c8")))(nrow(color_df_2))
  
  
  color_df_3 = data.frame("variable" = unique(df[str_detect(df[[reference_colum]], "Viruses"), reference_colum ] ) )
  color_df_3$col = colorRampPalette(c("#FFD700", "#CD853F" , "#7A2E00"))(nrow(color_df_3))
  
  
  other_df <- data.frame(variable=c('Other Bacteria', 'Other Viruses', "Other Eukaryota","Archaea"),
                         col=c("#00004d", "#331400","#021803", "#FEEFFA"))
  
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
path = "~/"

SRA_id = args[1]
dataset = args[2]

SRA_id="SRR14418894"
dataset="Douek_Cleveland"

# SRA_id="all"

print(dataset)
print(SRA_id)

path_result=paste0(path,"results/",dataset,"/")
path_data=paste0(path,"data/",dataset,"/")

sra_list=readLines(paste0(path_data,"/sra_list_RNA.txt"))
sra_metadata=read.table(paste0(path_data,"/sra_metadata.tsv"),header=TRUE)
sra_metadata_class=read.table(paste0(path_data,"/Metadata.csv"),header=TRUE, sep=",")


patient_id=sra_metadata[sra_metadata$Run==SRA_id,]$Title

# load(paste0(path_result,"/Taxon_correspondance.Rdata"))
# load(paste0(path_result,"/color_dict.Rdata"))
##########

result_dir_ku=paste0(path_result,"/kraken/microbialDB/kuniq/")
result_dir_k2=paste0(path_result,"/kraken/nt/k2uniq/")
result_dir_salmon=paste0(path_result,"/Contigs/")
result_dir_blastKu=paste0(path_result,"hybrid/hybrid_Blast-kuniq/")
result_dir_kuk2=paste0(path_result,"hybrid/hybrid_kuniq-kraken2/")

result_dir_spades=paste0(path_result,"/Contigs_rnaSpades/")
result_dir_spadesKU=paste0(path_result,"/hybrid/hybrid_spadesBlast-kuniq/")


### Import classifs
if (SRA_id %in% c("all","responder","non_responder") ){

  if (SRA_id == "responder"){
    sra_list = sra_metadata_class[sra_metadata_class$subject_classification == "Responder (>3yr on ART)",]$Run
  }else if (SRA_id == "non_responder"){
    sra_list = sra_metadata_class[sra_metadata_class$subject_classification != "Responder (>3yr on ART)",]$Run
  }

  classifs_list = parallel::mclapply(sra_list, function(sra){
    list=list(
      "KU"         = import_ReadsClassif(sra, result_dir_ku, "kraken", "KrakenUniq") ,
      # "kraken2" = import_ReadsClassif(sra, result_dir_k2, "kraken", "Kraken2")
      
      "blast"      = import_ReadsClassif(sra, result_dir_salmon, "salmon", "Blast") ,
      "spadesblast"= import_ReadsClassif(sra, result_dir_spades, "salmon", "spadesBlast") ,
      
      "spadesKU"   = import_ReadsClassif(sra, result_dir_spadesKU, "salmon", "spadesKU") ,
      "blastKU"    = import_ReadsClassif(sra, result_dir_blastKu, "hybrid", "hybrid_Blast_kuniq") ,
      
      "KuK2"       = import_ReadsClassif(sra, result_dir_kuk2, "hybrid", "KuK2") 
    )

    merged_data <- list$KU %>%
      
      #full_join(list$kraken2, by = "read") %>%
      full_join(list$blast, by = "read") %>%
      full_join(list$spadesblast, by = "read") %>%
      full_join(list$spadesKU, by = "read") %>%
      full_join(list$blastKU, by = "read") %>%
      full_join(list$KuK2, by = "read")

    # merged_data <- list$blast %>%
    #   # full_join(k2, by = "read") %>%
    #   full_join(list$KU, by = "read") %>%
    #   # full_join(list$blast, by = "read")
    #   full_join(list$kraken2, by = "read") %>%


    return(merged_data)

  }, mc.cores=28)

  names(classifs_list)=sra_list
  merged_df=do.call(rbind, classifs_list)
  rm(classifs_list)


} else {
  
  jobs <- list(
    KU           = list(result_dir_ku,       "kraken", "KrakenUniq"),
    #kraken2      = list(result_dir_k2,       "kraken", "Kraken2"),
    blast        = list(result_dir_salmon,   "salmon", "Blast"),
    spadesblast  = list(result_dir_spades,   "salmon", "spadesBlast"),
    spadesKU     = list(result_dir_spadesKU, "hybrid", "spadesKU"),
    blastKU      = list(result_dir_blastKu,  "hybrid", "hybrid_Blast_kuniq"),
    KuK2         = list(result_dir_kuk2,     "hybrid", "KuK2")
  )
  
  classifs <- parallel::mclapply(
    jobs,
    function(x) {
      import_ReadsClassif(
        SRA_id,
        x[[1]],  # result_dir
        x[[2]],  # method
        x[[3]]   # label
      ) }, mc.cores = 7)
  
  names(classifs) <- names(jobs)
  
  merged_df <- Reduce(function(x, y) merge(x, y, by = "read", all = TRUE), classifs)
  #rm(classifs)

}


tmp=merged_df[,c("read","hybrid_Blast_kuniq","KuK2")]
tmp = tmp[tmp$hybrid_Blast_kuniq != tmp$KuK2,]


read_save_dir=paste0(path_result,"specific_reads/")


### Change format of certain viruses:
# some viruses have no formal "phylum" associated, and tools returns "unclassfied" as their phylum, but are classified at lower level
taxons$phylum = ifelse( taxons$phylum == "unclassified Viruses phylum" , 
                        "Incertae sedis" , taxons$phylum)
# taxons$taxon_full_name = ifelse( taxons$taxon_full_name == "Viruses - unclassified Viruses phylum", 
#                                  "Viruses - Incertae sedis" , taxons$taxon_full_name)


# taxons$taxon_full_name = ifelse( taxons$taxon_full_name == "Viruses - Uroviricota", 
#                                  "Unclassified" , taxons$taxon_full_name)
# taxons$phylum = ifelse( taxons$phylum == "Uroviricota" , 
#                         "Unclassified" , taxons$phylum)
### Creating color palette
rename_taxa=do.call(rbind, list(
  "bacteria"=aggregate_taxa_count_df(taxons,5,"taxon_full_name","Freq", "Bacteria"),
  "viruses"=aggregate_taxa_count_df(taxons,4,"taxon_full_name","Freq", "Viruses"),
  "Archaea"=aggregate_taxa_count_df(taxons,1,"taxon_full_name","Freq", "Archaea"),
  "Euk"=aggregate_taxa_count_df(taxons,4,"taxon_full_name","Freq", "Eukaryota") ) )

as.character(taxons$taxon_full_name)[! as.character(taxons$taxon_full_name) %in% rename_taxa$taxon_full_name]

# Missing a few entries, re adding them bellow:
taxons[,c("superkingdom","phylum")]=NULL
taxons=taxons[! taxons$taxon_full_name %in% rename_taxa$taxon_full_name , ]
taxons$final_name=taxons$taxon_full_name
rename_taxa=rbind(taxons, rename_taxa)
row.names(rename_taxa) = NULL

#rename_taxa$final_name = ifelse(str_detect(rename_taxa$taxon_full_name, "myco"), "Eukaryota - Fungi" , as.character(rename_taxa$final_name) )

# Formating
if ( "Unclassified - Unclassified" %in% rename_taxa$taxon_full_name){
  rename_taxa[rename_taxa$taxon_full_name == "Unclassified - Unclassified",]$final_name = "Unclassified"
  #rename_taxa[rename_taxa$taxon_full_name == "unclassified - unclassified",]$final_name = "Unclassified"
  #rename_taxa[rename_taxa$taxon_full_name == "unassembled - unassembled",]$final_name = "Unassembled"
}

# Regroupe all archaea to a "Achaea" group
rename_taxa[str_detect(rename_taxa$taxon_full_name, "Archaea"),]$final_name = "Archaea"

count_na=function(x){return(table(is.na(x)))}
# apply(merged_df, MARGIN = 2, FUN = count_na)


rename_taxa$final_name = ifelse( rename_taxa$taxon_full_name == fixed("Viruses - unclassified Viruses phylum") , 
                        "Viruses - Incertae sedis" , rename_taxa$final_name)

cols_to_replace=colnames(merged_df)[-1]



load(file = paste0("~/results/color_dict_fig6.rdata"))

names(color_dict)=str_replace(names(color_dict), "-", ' - ')
setdiff(names(color_dict), unique(rename_taxa$final_name))



merged_df_og=merged_df
for (col in cols_to_replace) {
  print(col)
  merged_df[[col]]=rename_taxa$final_name[match(merged_df[[col]], as.character(rename_taxa$taxon_full_name))]
}
#apply(merged_df, MARGIN = 2, FUN = count_na)
# 
# color_dict=Taxonomy_small_color_palette(rename_taxa,"final_name")
# color_dict[["Unassembled"]]="#444444"
# color_dict[["Homo sapiens"]]="darkorchid"
# color_dict[["Unclassified"]]="#999999"
# 
tax=names(color_dict)
# 
# if(SRA_id == "SRR14418894"){
#   color1=color_dict[["Viruses - Incertae sedis"]]
#   color2=color_dict[["Viruses - Kitrinoviricota"]]
#   color_dict[["Viruses - Kitrinoviricota"]]=color1
#   color_dict[["Viruses - Incertae sedis"]]=color2
# }

#load("~/results/color_dict_25-05-25.RData")

tax_order=c(
  "Homo sapiens","Unassembled","Unclassified",
  tax[str_detect(tax,"Bacteria") & ! str_detect(tax, "Other")], "Other Bacteria","Archaea",
  tax[str_detect(tax,"Eukaryota") & ! str_detect(tax, "Other")], "Other Eukaryota",
  tax[str_detect(tax,"Viruses") & ! str_detect(tax, "Other")], "Other Viruses")

rename_taxa$final_name[!rename_taxa$final_name %in% tax_order]
tax_order[!tax_order %in% rename_taxa$final_name]



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# S A N K E Y   P L O T S
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Building dataframe for Sankey



#_______________________________________________________________________________
# 2) Remove Human/unclassified consensual 

to_remove=c("Homo sapiens","Unclassified")
cols_to_replace
df2_all <- merged_df[! (merged_df$hybrid_Blast_kuniq %in% to_remove & merged_df$KuK2 %in% to_remove),]%>% #& merged_df$Blast %in% to_remove),]%>%
  make_long(Blast, 
            #spadesBlast, 
            hybrid_Blast_kuniq, , 
            #spadesKU, 
            KuK2)





p2=ggplot(df2_all, aes(x = x, 
                       next_x = next_x, 
                       node = node, 
                       next_node = next_node,
                       fill = factor(node) )) + 
  geom_sankey(flow.alpha=0.85,node.color=1,space=100,size=0.2) +
  scale_fill_manual(values = color_dict) +
  #xlab(paste("Patient",patient_id))+  
  scale_x_discrete(labels = 
                     c("Trinity-Blast\nnt",
                       #"Spades-Blast\nnt",
                       "Hybrid-Trinity-KUniq\nnt + microbialDB",
                       #"Hybrid-Spades-KUniq\nnt + microbialDB",
                       "KUniq-K2\nmicrobialDB + nt"
                     )) +
  theme_sankey(base_size = 10 ) +
  theme(legend.position = "bottom", axis.text.x=element_text(size=10, face="bold"))
p2

ggsave(p2,filename=paste0("~/results/Simulation/figures_review/fig8.jpeg"),
       device = "jpeg",width = 30, height = 20,dpi = 600 , units = "cm", create.dir = TRUE)


ggsave(p2,
       filename = "~/results/Simulation/figures_review/fig8.svg",
       device = "svg",
       width = 30,
       height = 20,
       dpi = 600,
       units = "cm",
       create.dir = TRUE)

ggsave(p2,
       filename = "~/results/Simulation/figures_review/fig8.tiff",
       device = "tiff",compression = "lzw",
       width = 30, height = 20,
       dpi = 600,
       units = "cm",
       create.dir = TRUE)

