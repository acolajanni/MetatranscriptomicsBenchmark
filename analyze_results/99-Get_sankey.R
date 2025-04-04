#!/shared/ifbstor1/software/miniconda/envs/r-4.2.3/bin/Rscript

library(stringr)
library(plyr)
library(ggsankey)
library(ggplot2)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)

SRA_id = args[1]
dataset = args[2]

       
################################################################################
# F U N C T I O N S
################################################################################
ranks=c("superkingdom","phylum", "class", "order", "family", "genus", "species")


load_ReadsClassif=function(sra , path_results, data_name){ 
  ############
  # 1: Import data // minor differences between kraken and salmon
  if(data_name=="kraken"){ 
    rclassif_path=paste0(path_results, sra, "/ReadsClassif.txt") 
    x=read.table(rclassif_path, header=FALSE, fill = TRUE, sep = "\t",quote = "")
    colnames(x)=c("read", "taxid", "superkingdom","phylum", "class", "order", "family", "genus", "species")
    
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
  
  ############
  # 2: format data: columns are: readname + true taxid + ranks // TAXID lost with salmon/blast
  # taxid=1 ==> "root"
  if(data_name=="kraken"){ 
    x[x$taxid %in% c(1,198431,155900,131567), ranks] = "unclassified"
    x$taxid=NULL
  }
  x$taxid <- sub(".*taxid:([0-9]+)_.*", "\\1", x$read)
  
  return(x)
  
}
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
  color_df_1$col = colorRampPalette(c("sienna","lightyellow","steelblue3" ))(nrow(color_df_1))
  
  color_df_2 = data.frame("variable" = unique(df[
    str_detect(df[[reference_colum]], "Eukaryota") & !str_detect(df[[reference_colum]], "Fungi|myco"), reference_colum ] ) )
  
  
  color_df_2$col = colorRampPalette(c("plum1","firebrick1"))(nrow(color_df_2))
  
  color_df_3 = data.frame("variable" = unique(df[str_detect(df[[reference_colum]], "Viruses"), reference_colum ] ) )
  color_df_3$col = colorRampPalette(c("yellow2","darkorange"))(nrow(color_df_3))
  
  color_df_4 = data.frame("variable" = unique(df[str_detect(df[[reference_colum]], "Fungi|myco"), reference_colum ] ) )
  color_df_4$col = colorRampPalette(c("seagreen4","seagreen1"))(nrow(color_df_4))
  
  other_df <- data.frame(variable=c('Other Bacteria', 'Other Viruses', "Other Eukaryota","Other Eukaryota - Fungi", "Archaea"),
                         col=c("royalblue4", "#D5FF29","firebrick4","#084808", "slategray1"))
  
  color_df = rbind(color_df_1,color_df_2,color_df_3,color_df_4,other_df)
  
  color_dict = color_df$col
  names(color_dict) = color_df$variable
  return(color_dict) }


Taxonomy_small_color_palette = function(df, reference_colum){
  
  df=df[order(df[[reference_colum]]), ]
  df=df[!str_detect(df[[reference_colum]],"Other"), ]
  
  color_df_1 = data.frame("variable" = unique(df[str_detect(df[[reference_colum]], "Bacteria"), reference_colum  ] ) )
  color_df_1$col = colorRampPalette(c("aliceblue","royalblue4"  ))(nrow(color_df_1))
  
  color_df_2 = data.frame("variable" = unique(df[str_detect(df[[reference_colum]], "Eukaryota") , reference_colum ] ) )
  #color_df_2$col = colorRampPalette(c("thistle1","#960018"))(nrow(color_df_2))
  color_df_2$col = colorRampPalette(rev(c("#026e02","darkolivegreen1")))(nrow(color_df_2))
  
  
  color_df_3 = data.frame("variable" = unique(df[str_detect(df[[reference_colum]], "Viruses"), reference_colum ] ) )
  color_df_3$col = colorRampPalette(c("yellow2","darkorange2"))(nrow(color_df_3))
  
  
  other_df <- data.frame(variable=c('Other Bacteria', 'Other Viruses', "Other Eukaryota","Archaea"),
                         col=c("darkblue", "darkorange4","#084808", "#FEEFFA"))
  
  color_df = rbind(color_df_1,color_df_2,color_df_3,other_df)
  
  color_dict = color_df$col
  names(color_dict) = color_df$variable
  return(color_dict) }


################################################################################
# V A R I A B L E S
################################################################################
print(dataset)

path = "/run/user/1001/gvfs/sftp:host=core.cluster.france-bioinformatique.fr,user=acolajanni/shared/projects/microbiome_translocation/"
path = "/shared/projects/microbiome_translocation/"

# SRA_id="SRR14418884"
# SRA_id="all"
# #SRA_id="A21"
# dataset="Simulation"
# dataset="Douek_Cleveland"
print(dataset)

path_result=paste0(path,"results/",dataset,"/")
path_data=paste0(path,"data/",dataset,"/")

sra_list=readLines(paste0(path_data,"/sra_list_RNA.txt"))
sra_metadata=read.table(paste0(path_data,"/sra_metadata.tsv"),header=TRUE)

#sra_list=c("SRR14418868","SRR14418871","SRR14418859" )

patient_id=sra_metadata[sra_metadata$Run==SRA_id,]$Title

load(paste0(path_result,"/Taxon_correspondance.Rdata"))
load(paste0(path_result,"/color_dict.Rdata"))


path_salmon=paste0(path,"/results/",dataset,"/Contigs/")
path_kraken=paste0(path,"/results/",dataset,"/kraken/nt/k2uniq/")


if (SRA_id == "all"){
  print("importing data")
  
  data="salmon"
  tables_sa=parallel::mclapply(sra_list, function(sra){load_ReadsClassif(sra, path_salmon, data)},mc.cores=20)
  names(tables_sa)=sra_list
  # load(paste0(path_result,"Contigs/global_ReadsClassif.Rdata"))
  print("import Blast data: done")
  
  data="kraken"
  tables_ku=parallel::mclapply(sra_list, function(sra){load_ReadsClassif(sra, path_kraken, data)},mc.cores=20)
  #tables_ku=list()
  # for (sra in sra_list){
  #   print(sra)
  #   tables_ku[[sra]] = load_ReadsClassif(sra, path_kraken, data)
  #   print("_______________")
  # }
  names(tables_ku)=sra_list
  # load(paste0(path_result,"kraken/nt/global_ReadsClassif.Rdata"))
  print("import K2uniq data: done")
  
  print(names(tables_ku))
  
  classif_sa=unique(do.call(rbind, tables_sa))
  classif_ku=unique(do.call(rbind, tables_ku))

  classif_ku$taxid=NULL
  classif_sa$taxid=NULL

  save(classif_ku, file = paste0(path_result,"kraken/nt/global_ReadsClassif.Rdata"))
  save(classif_sa, file = paste0(path_result,"Contigs/global_ReadsClassif.Rdata"))

} else{
  
  result_dir_ku=paste0(path_result,"/kraken/microbialDB/kuniq/",SRA_id,"/")
  result_dir_salmon=paste0(path_result,"/Contigs/",SRA_id,"/")
  
  classif_ku=read.table(paste0(result_dir_ku,"ReadsClassif.txt"), fill = TRUE, sep="\t")
  colnames(classif_ku)=c("read","taxid","superkingdom","phylum", "class")
  classif_ku$read = sapply(classif_ku$read, function(x){str_replace_all(x, fixed("|"), fixed("-") )})
  classif_ku$read = sapply(classif_ku$read, function(x){str_split(x, fixed("/") )[[1]][1] })
  
  counting=as.data.frame(nrow(classif_ku))
  colnames(counting)="unmappedHuman"
  counting$SRA = SRA_id
  
  classif_sa=read.table(paste0(result_dir_salmon,"ContigsToReads/ReadsClassif.txt"),fill=TRUE, sep="\t")
  colnames(classif_sa)=c("contig","read","superkingdom","phylum", "class")
  classif_sa_raw=classif_sa
  classif_sa$class=NULL
  classif_sa$read = sapply(classif_sa$read, function(x){str_replace_all(x, fixed("|"), fixed("-") )})
  classif_sa_raw$read = sapply(classif_sa_raw$read, function(x){str_replace_all(x, fixed("|"), fixed("-") )})
  
  classif_sa_raw$read = sapply(classif_sa_raw$read, function(x){str_split(x, fixed("/") )[[1]][1] })
  classif_sa$read = sapply(classif_sa$read, function(x){str_split(x, fixed("/") )[[1]][1] })
  
  #classif_sa$phylum = ifelse(str_detect(classif_sa$phylum,"Fungi|myco"),"Fungi",classif_sa$phylum)
}

# unifying classification
classif_sa$taxon_full_name = paste0(classif_sa$superkingdom,"_",classif_sa$phylum)
classif_ku$taxon_full_name = paste0(classif_ku$superkingdom,"_",classif_ku$phylum)

classif_ku$taxon_full_name=ifelse(classif_ku$class == "Mammalia", "Homo_sapiens", classif_ku$taxon_full_name)
classif_sa$taxon_full_name=ifelse(classif_sa$class == "Mammalia", "Homo_sapiens", classif_sa$taxon_full_name)

classif_ku = classif_ku[,c("read","taxon_full_name")]
classif_sa = classif_sa[,c("read","taxon_full_name")]


## Creating an object with the global occurence of all the predicted taxa by the two methods
t=as.data.frame(table(classif_sa$taxon_full_name))
colnames(t)=c("taxon_full_name", "freq")
t$origin="blast"
t2=as.data.frame(table(classif_ku$taxon_full_name))
colnames(t2)=c("taxon_full_name", "freq")
t2$origin="ku"
t=rbind(t,t2)
t_merge=aggregate(freq ~ taxon_full_name, t, sum)
t_merge <- t_merge %>%
  separate(taxon_full_name, into = c("superkingdom", "phylum"), sep = "_", remove = FALSE)

# test=as.data.frame(rep(t_merge$taxon_full_name, t_merge$freq))
# colnames(test)="taxon_full_name"

tmp=list()
tmp[["Bacteria"]]=aggregate_taxa_count_df(t_merge, 1, "taxon_full_name", "freq", "Bacteria" )
tmp[["Viruses"]]=aggregate_taxa_count_df(t_merge, 2, "taxon_full_name", "freq", "Viruses" )
tmp[["Eukaryota"]]=aggregate_taxa_count_df(t_merge, 1, "taxon_full_name", "freq", "Eukaryota" )
tmp[["Archaea"]]=aggregate_taxa_count_df(t_merge, 1, "taxon_full_name", "freq", "Archaea" )


  
tmp=do.call(rbind, tmp)
tmp$final_name=ifelse( str_detect(tmp$taxon_full_name, "Archaea"), "Archaea" , as.character(tmp$final_name)  )
#tmp$final_name=ifelse( str_detect(tmp$taxon_full_name, "myco"), "Other Eukaryota - Fungi" , as.character(tmp$final_name)  )
tmp$final_name=ifelse( (str_detect(tmp$taxon_full_name, "myco") & str_detect(tmp$taxon_full_name, "Eukaryota")),
                       "Eukaryota - Fungi" , as.character(tmp$final_name)  )




missing_ids=t_merge$taxon_full_name[! t_merge$taxon_full_name %in% tmp$taxon_full_name]
missing_ids=t_merge[t_merge$taxon_full_name %in% missing_ids, c("taxon_full_name", "freq")]
missing_ids$final_name = ifelse(missing_ids$taxon_full_name == "Homo_sapiens", "Homo sapiens",t_merge$taxon_full_name)
missing_ids$final_name = ifelse(missing_ids$taxon_full_name == "unclassified_unclassified", "Unclassified",missing_ids$final_name)
missing_ids$final_name = ifelse(missing_ids$taxon_full_name == "unassembled_unassembled", "Unassembled",missing_ids$final_name)

keys=rbind(missing_ids, tmp)

color_dict=Taxonomy_small_color_palette(keys, "final_name")

color_dict[["Homo_sapiens"]]="purple"
color_dict[["Unclassified"]]="#999999"
color_dict[["Unassembled"]]="#444444"


keys$freq=NULL
keys$final_name = ifelse(keys$taxon_full_name == "Eukaryota_unclassified", "Other Eukaryota", keys$final_name)

unique(classif_sa$taxon_full_name[!classif_sa$taxon_full_name %in% keys$taxon_full_name])
unique(classif_ku$taxon_full_name[!classif_ku$taxon_full_name %in% keys$taxon_full_name])

classif_sa=merge(classif_sa,keys)
classif_ku=merge(classif_ku,keys)


# Creating unique taxonomy column
classif_ku$taxo_KU=classif_ku$final_name
classif_sa$taxo_SA=classif_sa$final_name

# Merging dataframe: comparing prediction of both methods
merged_classif=merge(classif_ku[, c("read", "taxo_KU")] ,
                     classif_sa[, c("read", "taxo_SA")] )


merged_classif$score = ifelse(merged_classif$taxo_KU == merged_classif$taxo_SA, 1, 0)


names(color_dict) = str_replace_all(names(color_dict), "_"," - ")
merged_classif$taxo_KU = str_replace_all(merged_classif$taxo_KU, "_"," - ")
merged_classif$taxo_SA = str_replace_all(merged_classif$taxo_SA, "_"," - ")

merged_classif$taxo_SA=str_replace_all(merged_classif$taxo_SA, "Bacteria - Other Bacteria", "Other Bacteria" )
merged_classif$taxo_KU=str_replace_all(merged_classif$taxo_KU, "Bacteria - Other Bacteria", "Other Bacteria" )

merged_classif$taxo_SA=str_replace_all(merged_classif$taxo_SA, "Eukaryota - Other Bacteria", "Other Eukaryota" )
merged_classif$taxo_KU=str_replace_all(merged_classif$taxo_KU, "Eukaryota - Other Bacteria", "Other Eukaryota" )

merged_classif$taxo_SA=str_replace_all(merged_classif$taxo_SA, "Homo - sapiens", "Homo sapiens" )
merged_classif$taxo_KU=str_replace_all(merged_classif$taxo_KU, "Homo - sapiens", "Homo sapiens" )
names(color_dict) = str_replace_all(names(color_dict), "Homo - sapiens","Homo sapiens")


taxa=unique(c(merged_classif$taxo_KU, merged_classif$taxo_SA))
missing_taxa=unique(taxa[!taxa %in% names(color_dict)])


merged_classif$flaged_reads=ifelse((merged_classif$taxo_KU == merged_classif$taxo_SA & merged_classif$taxo_KU == "Unclassified" ),1,0)
merged_classif$flaged_reads=ifelse((merged_classif$taxo_SA == "Unassembled" & merged_classif$taxo_KU == "Unclassified" ),1, merged_classif$flaged_reads)



#selected_reads=merged_classif[(merged_classif$taxo_KU != "" & merged_classif$taxo_KU != "Unclassified"),]$read
selected_reads=merged_classif[merged_classif$flaged_reads != 1,]$read


merged_classif_Blue=merged_classif[merged_classif$read %in% selected_reads,]

##### building alluvial / sankey plot
merged_classif_Blue$read=NULL
merged_classif_Blue$score=NULL
merged_classif_Blue$taxid=NULL
Freq_classif=as.data.frame(table(merged_classif_Blue))


names(color_dict)

merged_classif$score=NULL
# save(merged_classif, file = paste0(path_result,"/merged_ReadsClassif.Rdata"))
merged_classif$read=NULL

small_classif=merged_classif[ str_detect(merged_classif$taxo_KU, "Viruses|Bacteria") ,]
#save(small_classif, file = paste0(path_result,"/Bacteria_Viruses_ReadsClassif.Rdata"))



### Correcting color


color_dict[["Bacteria - Pseudomonadota"]] = colorRampPalette(c("aliceblue","royalblue4"  ))(6)[3]
color_dict[["Eukaryota - Fungi"]] = colorRampPalette(rev(c("#026e02","darkolivegreen1")))(3)[2]
color_dict[["Archaea"]] = "#FEEFFA"
color_dict[["Homo sapiens"]] = "darkorchid"

# 
# df <- merged_classif %>%
#   make_long(taxo_SA, taxo_KU)
# 
# save(df, file = paste0(path_result,"/sankey_df_ReadsClassif.Rdata"))


df <- merged_classif_Blue %>%
  make_long(taxo_SA, taxo_KU)

# save(df, file=paste0(path_result,"/small_sankey_df_ReadsClassif.Rdata"))


p1=ggplot(df, aes(x = x, 
                  next_x = next_x, 
                  node = node, 
                  next_node = next_node,
                  fill = factor(node) )) + 
  geom_sankey(flow.alpha=0.85,node.color=1,space=150) +
  scale_fill_manual(values = color_dict) +
  xlab(paste("Patient",patient_id))+  
  scale_x_discrete(labels = c("Blue pipeline","Red pipeline")) +
  theme_sankey(base_size = 16 ) 
p1

dir.create(paste0(path_result,"/sankey_plots/"), showWarnings = FALSE)

ggsave(p1,filename=paste0(path_result,"/sankey_plots/",SRA_id,"_sankey_full.png"),
       device = "png",height = 24, width = 45, units = "cm")

stop()


################################################################################
# R E S E R V E
################################################################################


# test=df[!(df$x == "taxo_SA" & df$node == "Unclassified"),]
# 
# df2=df[!(df$x == "taxo_SA" & df$node == "Unassembled"),]
# ggplot(df2, aes(x = x, 
#                next_x = next_x, 
#                node = node, 
#                next_node = next_node,
#                fill = factor(node) )) + 
#   geom_sankey(flow.alpha=0.85,node.color=1,space=150) +
#   scale_fill_manual(values = color_dict) +
#   xlab(paste("Patient",patient_id))+  scale_x_discrete(labels = c("Blue pipeline","Red pipeline")) +
#   theme_sankey(base_size = 16, ) 

# df2 <- merged_classif_Blue[merged_classif_Blue$taxo_SA == "Eukaryota - Fungi",] %>%
#   make_long(taxo_SA, taxo_KU)
# p2=ggplot(df2, aes(x = x, 
#                    next_x = next_x, 
#                    node = node, 
#                    next_node = next_node,
#                    fill = factor(node) )) + 
#   geom_sankey(flow.alpha=0.85,node.color=1,space=150) +
#   scale_fill_manual(values = color_dict) +
#   xlab(paste("Patient",patient_id))+  scale_x_discrete(labels = c("Blue pipeline","Red pipeline")) +
#   theme_sankey(base_size = 16, ) 
# p2




# ggsave(p2,filename=paste0(path_result,"/sankey_plots/",SRA_id,"_fungi_redblue.png"),
#        device = "png",height = 24, width = 40, units = "cm")




classif_ku$taxid=sapply(classif_ku$read, function(x){ str_split(x, pattern = ":|-|_")[[1]][2] })
classif_sa$taxid=sapply(classif_sa$read, function(x){ str_split(x, pattern = ":|-|_")[[1]][2] })

sequence_database=read.delim(paste0(path,"database/Sequence_database.tsv"),header=FALSE)
sequence_database=unique(sequence_database[,c(2,4,5)])
colnames(sequence_database)=c("taxid","true_s","true_p")
sequence_database$taxo_truth=paste0(sequence_database$true_s,"_",sequence_database$true_p)


taxid=merge(classif_ku,sequence_database)
taxid=unique(taxid[,c(1,8,9)])

classif_sa$taxo_SA = classif_sa$taxon_full_name
classif_ku$taxo_KU = classif_ku$taxon_full_name

merge_classif_truth=merge(classif_ku[,c("read","taxid","taxo_KU")], classif_sa[,c("read","taxid","taxo_SA")], by="read")
merge_classif_truth$taxid=merge_classif_truth$taxid.y
merge_classif_truth[,c(2,4)]=NULL

merge_classif_truth=merge(merge_classif_truth, sequence_database[,c(1,4)], by="taxid")
merge_classif_truth=merge_classif_truth[,c(3,4,5)]

colnames(merge_classif_truth) = c("taxo_KU" , "taxo_SA" , "taxo_truth" )
###

color_dict[["Homo_sapiens"]]="purple"
color_dict[["Unclassified"]]="#777777"

#color_dict[["Bacillus cereus m1293"]]="red

names(color_dict) = str_replace_all(names(color_dict), "_"," - ")
merge_classif_truth$taxo_KU = str_replace_all(merge_classif_truth$taxo_KU, "_"," - ")
merge_classif_truth$taxo_SA = str_replace_all(merge_classif_truth$taxo_SA, "_"," - ")
merge_classif_truth$taxo_truth = str_replace_all(merge_classif_truth$taxo_truth, "_"," - ")

merge_classif_truth$taxo_SA=str_replace_all(merge_classif_truth$taxo_SA, "Bacteria - Other Bacteria", "Other Bacteria" )
merge_classif_truth$taxo_KU=str_replace_all(merge_classif_truth$taxo_KU, "Bacteria - Other Bacteria", "Other Bacteria" )
merge_classif_truth$taxo_truth=str_replace_all(merge_classif_truth$taxo_truth, "Bacteria - Other Bacteria", "Other Bacteria" )

merge_classif_truth$taxo_SA=str_replace_all(merge_classif_truth$taxo_SA, "Eukaryota - Other Bacteria", "Other Eukaryota" )
merge_classif_truth$taxo_KU=str_replace_all(merge_classif_truth$taxo_KU, "Eukaryota - Other Bacteria", "Other Eukaryota" )
merge_classif_truth$taxo_truth=str_replace_all(merge_classif_truth$taxo_truth, "Eukaryota - Other Bacteria", "Other Eukaryota" )

merge_classif_truth$taxo_truth=str_replace_all(merge_classif_truth$taxo_truth, "Homo - sapiens", "Homo sapiens" )
merge_classif_truth$taxo_SA=str_replace_all(merge_classif_truth$taxo_SA, "Homo - sapiens", "Homo sapiens" )
merge_classif_truth$taxo_KU=str_replace_all(merge_classif_truth$taxo_KU, "Homo - sapiens", "Homo sapiens" )

#merge_classif_truth$taxo_truth=ifelse(merge_classif_truth$taxo_truth %in% c("Eukaryota - Ascomycota", "Eukaryota - Basidiomycota"), 
#                                      "Eukaryota - Fungi", merge_classif_truth$taxo_truth )



names(color_dict) = str_replace_all(names(color_dict), "Homo - sapiens","Homo sapiens")
color_dict["Eukaryota - Ascomycota"] = "seagreen4"
color_dict["Eukaryota - Basidiomycota"] = "seagreen3"
#color_dict["Eukaryota - Fungi"] = NULL

merge_classif_truth$taxo_KU = str_replace(merge_classif_truth$taxo_KU, "unclassified - unclassified", "Unclassified")
merge_classif_truth$taxo_SA = ifelse(merge_classif_truth$taxo_SA==" - ", "Unclassified", merge_classif_truth$taxo_SA)
  
  


df2 <- merge_classif_truth %>%
  make_long(taxo_SA, taxo_truth, taxo_KU)




p1=ggplot(df2, aes(x = x, 
                   next_x = next_x, 
                   node = node, 
                   next_node = next_node,
                   fill = factor(node) )) + 
  geom_sankey(flow.alpha=0.85,node.color=1,space=150) +
  scale_fill_manual(values = color_dict) +
  xlab(paste("Patient",patient_id))+  scale_x_discrete(labels = c("Blue pipeline","Ground Truth","Red pipeline")) +
  theme_sankey(base_size = 16, ) 

p1

ggsave(p1,filename=paste0(path_result,"/sankey_plots/",SRA_id,"_sankey_truth.png"),
       device = "png",height = 24, width = 40, units = "cm")




# fungi_classif=Freq_classif[Freq_classif$taxo_SA=="Eukaryota - Fungi",]
# 
# fungi_classif <- fungi_classif %>%
#   arrange(desc(taxo_KU)) %>%
#   mutate(lab.ypos = cumsum(Freq) - 0.5*Freq)
# fungi_classif
# 
# ggplot(fungi_classif, aes(x="", y=Freq, fill=taxo_KU))+
#   geom_bar(width = 1, stat = "identity", color = "white") +
#   coord_polar("y", start = 0)+
#   geom_text(aes(y = lab.ypos,label = Freq), color = "white")+
#   scale_fill_manual(values = color_dict) +
#   theme_void()
# 
# ##########################
# interesting_reads=merged_classif_Blue[merged_classif_Blue$taxid == "526973",]$read
# 
# 
# r1=paste0(interesting_reads, "/1")
# r2=paste0(interesting_reads, "/2")
# 
# writeLines(r1, paste0(result_dir_ku,"inconsistent_r1.txt"))
# writeLines(r2, paste0(result_dir_ku,"inconsistent_r2.txt"))
# 
# 
# 
# writeLines(interesting_reads, paste0(result_dir_ku,"inconsistent_reads.txt"))
# 
# 
# classif_ku=read.delim(paste0(result_dir_ku,SRA_id,".output.txt"),header=FALSE)
# classif_ku=classif_ku[classif_ku$V2 %in% merged_fungi$read,]
# classif_ku$V5=NULL
# classif_ku$V4=NULL
# classif_ku$V1=NULL
# colnames(classif_ku)=c("read","taxid")
# classif_ku = merge(merged_classif, classif_ku)
# 
# id=names(table(classif_ku$taxid)[table(classif_ku$taxid) == max(table(classif_ku$taxid))])
# classif_inconsistent = classif_ku[classif_ku$taxid == id, ]
# 
# r1=paste0(classif_inconsistent$read, "/1")
# r2=paste0(classif_inconsistent$read, "/2")
# 
# writeLines(r1, paste0(result_dir_ku,"inconsistent_r1.txt"))
# writeLines(r2, paste0(result_dir_ku,"inconsistent_r2.txt"))
# 
# 
# load("/home/acolajanni/Documents/work/results/Douek_Cleveland/quant_lvl_Cleveland.rdata")
# load("/home/acolajanni/Documents/work/results/Douek_Cleveland/quant_lvl_kuniq_RNA.rdata")



