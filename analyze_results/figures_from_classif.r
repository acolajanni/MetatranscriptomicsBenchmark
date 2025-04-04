#!/shared/ifbstor1/software/miniconda/envs/r-4.2.3/bin/Rscript

library(stringr)
library(plyr)
library(ggsankey)
library(ggplot2)
library(dplyr)
library(cowplot)


rename_taxa=function(classif_df, taxa_column="taxon_phylum", count_column="NumReads", sep="-", nbact=5, nvir=3, neuk=3){
  
  classif_df$NumReads=classif_df[[count_column]]
  classif_df$taxon_phylum=classif_df[[taxa_column]]
  
  t = aggregate(classif_df$NumReads, by = list("taxon_phylum"=classif_df$taxon_phylum), FUN = sum)
  colnames(t)[colnames(t) == 'x'] = "NumReads"
  #t=aggregate(NumReads ~ taxon_phylum, classif_df, FUN = sum)
  t=t[order(t$NumReads, decreasing = TRUE),]
  t[,c("superkingdom","phylum")]=str_split_fixed(t$taxon_phylum,pattern = fixed(sep), n=2)
  bact=aggregate_taxa_count_df(t,n=nbact,"taxon_phylum","NumReads","Bacteria")
  vir=aggregate_taxa_count_df(t,n=nvir,"taxon_phylum","NumReads","Viruses")
  euk=aggregate_taxa_count_df(t,n=neuk,"taxon_phylum","NumReads","Eukaryota")
  
  archaea=t[t$superkingdom == "Archaea",]
  if(nrow(archaea)>0 ){ archaea$final_name="Archaea"}
  
  colnames(archaea)[1] = "taxon_full_name"
  archaea[,c("superkingdom","phylum")]=NULL
  
  key_from_data=do.call(rbind, list(bact,euk,vir,archaea))
  key_from_data$NumReads=NULL
  
  missing_id=unique(classif_df$taxon_phylum[!classif_df$taxon_phylum %in% key_from_data$taxon_full_name ])
  
  hs = data.frame("taxon_full_name" = missing_id, "final_name" = missing_id)
  key_from_data=rbind(key_from_data, hs)
  
  return(key_from_data)
}

aggregate_taxa_count_df=function(table_df, n=5, column_to_filter, Frequency_column,replace_taxon="Bacteria"){
  table_df=table_df[table_df$superkingdom == replace_taxon , ]
  table_df=table_df[order( table_df[[Frequency_column]],decreasing = TRUE ) , ]
  
  
  if (replace_taxon == "Eukaryota"){
    
    print("eukaryota")
    fungi_df=table_df[str_detect(table_df[[column_to_filter]], "Fungi|fungi|myco" ),]
    fungi_df$superkingdom = "Eukaryota - Fungi"
    
    fun=aggregate_taxa_count_df(fungi_df, 4, column_to_filter, Frequency_column, "Eukaryota - Fungi"  )
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


format_read_table=function(read_table, column, tax_lvl="phylum", split_char="|"){ 
  ranks=c("superkingdom","phylum", "class", "order", "family", "genus", "species")
  ranks_to_keep=ranks[1:which(ranks==tax_lvl)]
  
  df=as.data.frame(table(read_table[[column]]))
  colname=paste0("taxon_",tax_lvl)
  colnames(df)=c(colname,"NumReads")
  df[,ranks_to_keep] = str_split_fixed(df[[colname]], fixed(split_char), n=length(ranks_to_keep) )
  return(df)
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

format_human_taxon=function(classif_df){
  human_df=classif_df[classif_df$class == "Mammalia",]
  human_df$taxon="Homo_sapiens"
  return(human_df)
}


Taxonomy_small_color_palette2 = function(df, reference_colum){
  
  df=df[order(df[[reference_colum]]), ]
  df=df[!str_detect(df[[reference_colum]],"Other"), ]
  
  color_df_1 = data.frame("variable" = unique(df[str_detect(df[[reference_colum]], "Bacteria"), reference_colum  ] ) )
  color_df_1$col = colorRampPalette(c("aliceblue","royalblue4" ))(nrow(color_df_1))
  
  color_df_2 = data.frame("variable" = unique(df[str_detect(df[[reference_colum]], "Eukaryota") , reference_colum ] ) )
  #color_df_2$col = colorRampPalette(c("thistle1","#960018"))(nrow(color_df_2))
  color_df_2$col = colorRampPalette(rev(c("#026e02","darkolivegreen1")))(nrow(color_df_2))
  
  
  
  color_df_3 = data.frame("variable" = unique(df[str_detect(df[[reference_colum]], "Viruses"), reference_colum ] ) )
  color_df_3$col = colorRampPalette(c("yellow","darkorange2"))(nrow(color_df_3))
  
  #color_df_4$col = colorRampPalette(c("seagreen4","seagreen1"))(nrow(color_df_4))
  
  other_df <- data.frame(variable=c('Other Bacteria', 'Other Viruses', "Other Eukaryota", "Archaea"),
                         col=c("darkblue","darkorange4","#084808","#FEEFFA"))
  
  color_df = rbind(color_df_1,color_df_2,color_df_3,other_df)
  
  color_dict = color_df$col
  names(color_dict) = color_df$variable
  return(color_dict) }


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
  classif_merged$taxon=ifelse(classif_merged$taxon == "not_in_contigs|not_in_contigs",
                              'Not in Contigs', classif_merged$taxon)
  return(classif_merged )
}

aggregate_taxa_repeated=function(table_df, n=5, replace_taxon="Bacteria"){
  top=head(table_df, n = n)
  low=table_df[!table_df$Var1 %in% top$Var1 , ]
  
  top$final_name=top$Var1
  low$final_name=sapply(low$superkingdom, function(x){paste0("Other ",replace_taxon) })
  
  table_df=rbind(top,low)
  table_df[,c("Freq","superkingdom","phylum")]=NULL
  colnames(table_df)[1]="taxon_full_name"
  
  return(table_df)
}


dataset="Douek_Cleveland"
path = "/home/acolajanni/Documents/work/"
path="/run/user/1001/gvfs/sftp:host=core.cluster.france-bioinformatique.fr,user=acolajanni/shared/projects/microbiome_translocation/"
path="/shared/projects/microbiome_translocation/"
print(dataset)

path_result=paste0(path,"results/",dataset)
path_data=paste0(path,"data/",dataset)

sra_list=readLines(paste0(path_data,"/sra_list_RNA.txt"))
sra_metadata=read.table(paste0(path_data,"/sra_metadata.tsv"),header=TRUE)


load(paste0(path_result,"/Taxon_correspondance.Rdata"))
load(paste0(path_result,"/color_dict.Rdata"))

# à corriger
#sra_list=sra_list[! sra_list %in% c("SRR14418861","SRR14418873", "SRR14418888", "SRR14418889", "SRR14418904") ]




classif=parallel::mclapply(sra_list, function(sra){  
  path_rdata=paste0(path_result,"/rdata/")
  print(sra)
  load(paste0(path_rdata, sra, "_classif_merged.rdata"))  
  return(classif_merged)  
  #rm(classif_merged)
},mc.cores=10)
names(classif) = sra_list

# classif=list()
# for (sra in sra_list){
#   path_rdata=paste0(path_result,"/rdata/")
#   print(sra)
#   load(paste0(path_rdata, sra, "_classif_merged.rdata"))  
#   classif[[sra]]=classif_merged  
#   rm(classif_merged)
# }

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Distribution of the number of sequences remaining after filtration

n_reads=t(as.data.frame(lapply(classif, function(x) nrow(x))))
colnames(n_reads)= "read_number"

ggplot(n_reads, aes(x = read_number)) +
  geom_histogram(aes(y = ..density..), bins=12,
                 colour = 1, fill = "steelblue", alpha=0.7) +
  geom_density(linewidth=1.5)+  
  labs(
    title = "Simulated Distribution Based on Summary Statistics", 
    x = "Values", 
    y = "Frequency"
  ) +
  theme_minimal() +  # A clean and minimal theme
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(color = "gray", size = 0.2),
    panel.grid.minor = element_blank()  # Remove minor grid lines
  ) + 
  ggtitle("Distribution of the number of unmapped reads") + 
  scale_x_continuous(breaks = seq(0, 4500000, by = 500000)) +  # Custom x-axis breaks
  scale_y_continuous(expand = c(0, 0))  # Remove space at the bottom of the plot

ggplot(n_reads, aes(x = read_number)) +
  geom_histogram(bins=16,colour = 1, fill = "steelblue", alpha=0.7) +
  geom_vline(aes(xintercept = median(read_number)), color = "#000000", size = 1.25) +
  geom_vline(aes(xintercept = median(read_number) + IQR(read_number)), color = "#000000", size = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = median(read_number) - IQR(read_number)), color = "#000000", size = 1, linetype = "dashed") + 
  annotate("text" , 
           x=median(n_reads[,1])*1.25, 
           y = 16 , 
           label=paste0( "median: ", '\n', as.character(round(median(n_reads[,1])))), size=10 ) +
  labs(
    x = "Number of unmapped reads against human genomes", 
    y = "Frequency"
  ) +
  theme_bw() +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(color = "gray", size = 0.2),
    panel.grid.minor = element_blank()  # Remove minor grid lines
  ) + 
  scale_x_continuous(labels = scales::comma) +
  ggtitle("Distribution of the number of unmapped reads against human genomes", 
          subtitle = "Annotation: median and median +/- IQR") + 
  scale_x_continuous(breaks = seq(0, 4500000, by = 500000)) +  # Custom x-axis breaks
  scale_y_continuous(expand = c(0, 0))  # Remove space at the bottom of the plot









### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

### counting where both methods agrees

agree_score=lapply(classif, function(x) {
  x$score = ifelse(x$taxon_K2Uniq == x$taxon_Blast, yes=1, no = 0) 
  # Counting occurences 
  x$read = NULL
  x <- x %>%
    group_by(across(everything())) %>%  # Group by all columns
    summarise(NumReads = n()) %>%          # Count the occurrences of each row
    ungroup()
  return(x) })

agree_score=lapply(agree_score, function(x) {
  x$transition = paste0(x$taxon_K2Uniq, " --> ",x$taxon_Blast)
  return(x)
})

agree_score=ldply(agree_score)
#agree_score$.id = NULL

agg=aggregate(NumReads ~ transition + score,data=agree_score, FUN = sum)
agg$Freq=agg$NumReads / sum(agg$NumReads)
agg$score = as.factor(agg$score)
agg=agg[order(-agg$Freq),]

ggplot(head(agg,20), aes(y=reorder(transition,Freq), x=Freq, fill=score) ) + 
  geom_bar(stat="identity") + 
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size=12),
    legend.position = c(0.85, 0.1),
    panel.grid.major = element_line(color = "gray", size = 0.2),
    panel.grid.minor = element_blank()  # Remove minor grid lines
  ) + labs(
    title = "Frequency of classification at Phylum level by the two methods", 
    x = "Frequency of Read", 
    y = "Taxon predction: K2Uniq ==> Blast" ) + 
  scale_fill_manual(
    values = c("0" = "red3", "1" = "steelblue"),   # Define colors for 0 and 1
    labels = c("0" = "Different prediction", "1" = "Same Prediction")  ,    # Labels for 0 and 1 in the legend
    name=""
  )

### Regroup by category : A ==> B and B ==> A
agg$A=sapply(agg$transition, function(x) str_split(x, " --> ")[[1]][1]  )
agg$B=sapply(agg$transition, function(x) str_split(x, " --> ")[[1]][2]  )


#agg$transition_name=ifelse(agg$score == "1", yes= agg$A, no = agg$transition)

agg$transition2=paste0(agg$B, " --> ", agg$A)
taxa=unique(c(agg$A, agg$B))




taxa_df=agg[NULL, c("score","NumReads","A","B","transition")]
pair=vector()
i=1
for (taxon in taxa ){
  for (taxon2 in taxa ){
    
    
    # has the pair already been seen
    if (paste0(taxon2,'-',taxon) %in% pair){ next }
    
    if (taxon == taxon2){
      #tmp=agg[agg$A==taxon & agg$B==taxon,c("score","NumReads","A","B")]
      next
    } 
    else{
      
      #print(c(taxon,taxon2))
      solutions=c(paste0(taxon, " --> ", taxon2),
                  paste0(taxon2, " --> ", taxon))
      
      tmp=agg[agg$transition %in% solutions,c("score","NumReads","A","B","transition")]
      
      
      if (nrow(tmp) == 0 ) { 
        #print("unexistent combination")
        next } # no interesting rows
      
      tmp=aggregate(NumReads ~ score, tmp, FUN = sum)
      tmp$A = taxon
      tmp$B = taxon2
      tmp$transition=paste0(tmp$A," --> ",tmp$B)
    }
    
    taxa_df= rbind(taxa_df,tmp)
    
    pair[i]=paste0(taxon,'-',taxon2)
    i=i+1
  }
}

taxa_df=unique(taxa_df)
taxa_df=rbind(taxa_df, agg[agg$score == "1", c("score","NumReads","A","B","transition") ])

taxa_df=taxa_df[order(-taxa_df$NumReads),]
taxa_df$Freq = taxa_df$NumReads / sum(taxa_df$NumReads)


ggplot(head(taxa_df,20), aes(y=reorder(transition,Freq), x=Freq, fill=score) ) + 
  geom_bar(stat="identity") + 
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size=12),
    legend.position = c(0.85, 0.1),
    panel.grid.major = element_line(color = "gray", size = 0.2),
    panel.grid.minor = element_blank()  # Remove minor grid lines
  ) + labs(
    title = "Frequency of classification at Phylum level by the two methods", 
    x = "Frequency of Read", 
    y = "Taxon predction: Most frequent pair" ) + 
  scale_fill_manual(
    values = c("0" = "red3", "1" = "steelblue"),   # Define colors for 0 and 1
    labels = c("0" = "Different prediction", "1" = "Same Prediction")  ,    # Labels for 0 and 1 in the legend
    name=""
  )

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

### Same but without some eukaryotes, unclassified, etc: don't take the taxa_df


reduced=agg[!agg$A %in% 
              c("unclassified","artificial sequences","cellular organisms"
                ,"Eukaryota|Chordata","Eukaryota|unclassified Eukaryota phylum",
                "Eukaryota|Arthropoda","Eukaryota|Streptophyta","Eukaryota|Mollusca") ,]

t=aggregate(NumReads ~ A, reduced, FUN = sum)
t=t[order(t$NumReads, decreasing = TRUE),]
t[,c("superkingdom","phylum")]=str_split_fixed(t$A,pattern = fixed("|"), n=2)



bact=aggregate_taxa_count_df(t,n=5,"A","NumReads","Bacteria")
vir=aggregate_taxa_count_df(t,n=3,"A","NumReads","Viruses")
euk=aggregate_taxa_count_df(t,n=5,"A","NumReads","Eukaryota")
archaea=t[t$superkingdom == "Archaea",]
archaea$final_name="Archaea"
colnames(archaea)[1] = "taxon_full_name"
archaea[,c("superkingdom","phylum")]=NULL

key_from_data=do.call(rbind, list(bact,euk,vir,archaea))
key_from_data$NumReads=NULL
hs = data.frame("taxon_full_name" = "Homo_sapiens", "final_name" = "Homo_sapiens")
key_from_data=rbind(key_from_data, hs)

#count_df_reduced=aggregate(NumReads ~ final_name, count_df_reduced, FUN = sum)



reduced=merge(key_from_data, reduced, by.x='taxon_full_name',by.y="A")

reduced=reduced[!reduced$B %in% 
                  c("unclassified","artificial sequences","cellular organisms", "Not in Contigs"
                    ,"Eukaryota|Chordata","Eukaryota|unclassified Eukaryota phylum",
                    "Eukaryota|Arthropoda","Eukaryota|Streptophyta","Eukaryota|Mollusca") ,]

reduced$transition = paste0(reduced$final_name," --> " ,reduced$B)
reduced=reduced[order(reduced$Freq, decreasing = TRUE),]

library(ggbreak) 
library(patchwork)

ggplot(head(reduced,25), aes(y=reorder(transition,Freq), x=Freq, fill=score) ) + 
  geom_bar(stat="identity") + 
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0, size = 16,face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size=12),
    #legend.position = c(0, 0.1),
    panel.grid.major = element_line(color = "gray", size = 0.2),
    panel.grid.minor = element_blank()  # Remove minor grid lines
  ) + 
  scale_x_break( c(0.005, 0.015), scale=.45 ) +
  labs(
    title = "Frequency of classification at Phylum level by the two methods (with unclassified, and spurious eukaryota removed)", 
    x = "Frequency of Read", 
    y = "Taxon predction: K2Uniq ==> Blast" ) + 
  scale_fill_manual(
    values = c("0" = "red3", "1" = "steelblue"),   # Define colors for 0 and 1
    labels = c("0" = "Different prediction", "1" = "Same Prediction")  ,    # Labels for 0 and 1 in the legend
    name=""
  )



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Proportions as in the poster

### Trying rebuilding quant_lvl item: quantification for each level of classification
ranks=c("superkingdom","phylum", "class", "order", "family", "genus")
ranks=c("superkingdom","phylum")


classif_smaller_list=lapply(classif,function(x){
  k2=format_read_table(x, "taxon_K2Uniq")  
  blast=format_read_table(x, "taxon_Blast")
  l=list("K2uniq"=k2, "Blast"=blast)
  l=ldply(l)
  colnames(l)[1]="Method"
  return(l)
})

##############################################
## replace "unclassified bacteria phylum" 
## by "unclassified"
classif_smaller=ldply(classif_smaller_list)

classif_smaller$phylum=ifelse(str_detect(classif_smaller$phylum, "unclassified") , "unclassified", classif_smaller$phylum)
classif_smaller$taxon_phylum = ifelse(classif_smaller$phylum == "", 
                                      yes=classif_smaller$superkingdom,
                                      no = paste0(classif_smaller$superkingdom, "-",classif_smaller$phylum))
# Collapse read number by the "new" phylum
classif_smaller=aggregate(NumReads ~ .id + Method + taxon_phylum + superkingdom + phylum, classif_smaller, sum)

classif_fullnames=classif_smaller
################################################################################
################################################################################
key_from_data=rename_taxa(classif_smaller, nbact = 5,nvir = 3,neuk = 5)

# t=aggregate(NumReads ~ taxon_phylum, classif_smaller, FUN = sum)
# t=t[order(t$NumReads, decreasing = TRUE),]
# t[,c("superkingdom","phylum")]=str_split_fixed(t$taxon_phylum,pattern = fixed("-"), n=2)
#   
# bact=aggregate_taxa_count_df(t,n=5,"taxon_phylum","NumReads","Bacteria")
# vir=aggregate_taxa_count_df(t,n=3,"taxon_phylum","NumReads","Viruses")
# euk=aggregate_taxa_count_df(t,n=5,"taxon_phylum","NumReads","Eukaryota")
# 
# archaea=t[t$superkingdom == "Archaea",]
# archaea$final_name="Archaea"
# colnames(archaea)[1] = "taxon_full_name"
# archaea[,c("superkingdom","phylum")]=NULL
# 
# key_from_data=do.call(rbind, list(bact,euk,vir,archaea))
# key_from_data$NumReads=NULL
# 
# # Adding missing keys to the dataframe
# missing_id=unique(classif_smaller$taxon_phylum[!classif_smaller$taxon_phylum %in% key_from_data$taxon_full_name ])
# 
# hs = data.frame("taxon_full_name" = missing_id, "final_name" = missing_id)
# key_from_data=rbind(key_from_data, hs)
# 
# unique(key_from_data)
################################################################################


# Creating color_palette:
color_dict=Taxonomy_small_color_palette(key_from_data, "final_name")
key_from_data$final_name[!key_from_data$final_name %in% names(color_dict)]


# color_dict[["Other Bacteria"]]="royalblue4"
# color_dict[["Other Viruses"]]="#D5FF29"
# color_dict[["Other Eukaryota"]]="firebrick4"
# color_dict[["Other Eukaryota - Fungi"]]="#084808"
# color_dict[["Archaea"]]="azure"
color_dict[["Unassembled"]]="#444444"
color_dict[["Homo_sapiens"]]="darkorchid"
color_dict[["unclassified"]]="#999999"
color_dict[["cellular organisms"]]="#777777"
color_dict[["artificial sequences"]]="#CCCCCC"


# change classification column 
classif_smaller=merge(classif_smaller, key_from_data, by.x = "taxon_phylum", by.y="taxon_full_name")
classif_smaller=aggregate(NumReads ~ .id + Method + final_name, classif_smaller, sum)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### barplots


ku=classif_smaller[classif_smaller$Method=="K2uniq",]

classif_smaller <- classif_smaller %>%
  group_by(.id,Method) %>%
  mutate(total = sum(NumReads)) %>%
  ungroup()  # Remove grouping for further operations
classif_smaller$Freq = classif_smaller$NumReads/classif_smaller$total

classif_smaller=unique(classif_smaller)

tax=unique(classif_smaller$final_name)
tax_order=c(
  "artificial sequences" , "unclassified","cellular organisms", "Unassembled", "Homo_sapiens",
  tax[str_detect(tax,"Bacteria") & ! str_detect(tax, "Other")], "Other Bacteria",
  tax[str_detect(tax,"myco") & ! str_detect(tax, "Other")], "Other Eukaryota - Fungi",
  tax[str_detect(tax,"Eukaryota") & ! str_detect(tax, "Other") & ! str_detect(tax,"myco")], "Other Eukaryota",
  "Archaea",
  tax[str_detect(tax,"Viruses") & ! str_detect(tax, "Other")], "Other Viruses")

classif_smaller$final_name <- factor(classif_smaller$final_name, levels=tax_order)


order_patient=rownames(n_reads)[order(n_reads[,1], decreasing = FALSE)]
classif_smaller$.id <- factor(classif_smaller$.id, levels=order_patient)


ggplot(classif_smaller, aes(x=.id, y=Freq, fill = final_name) ) +
  geom_bar(position="stack", stat="identity", width = 0.95)+#color="darkgrey") +
  theme_linedraw() + scale_fill_manual(values = color_dict) + 
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
  guides(fill=guide_legend("Superkingdom - Phylum" , ncol=1)) +
  xlab("Individual") + 
  ylab(paste("Read frequency")) + 
  theme(plot.title = element_text(hjust = 0.5, size=16), strip.text = element_text(size=18), 
        legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key.size =unit(.6,"cm"),
        axis.text.x=element_text(angle=30,size=8),
        panel.background = element_rect(fill = "#111111"),    
        panel.grid.major = element_line(colour = "grey",),
        panel.grid.minor = element_line(colour = "grey",)) +
  coord_flip() + 
  facet_wrap(~ Method, nrow = 1, scales = 'free_x' )

counting=aggregate(NumReads ~ Method+final_name, classif_smaller[classif_smaller$Method == "Blast",], sum)
### Mean number of reads remaining after removing unassembled sequences
sum(counting[counting$final_name!="Unassembled",]$NumReads)/55
sum(counting$NumReads)/55

ggplot(classif_smaller, aes(x=.id, y=NumReads, fill = final_name) ) +
  geom_bar(position="stack", stat="identity", width = 0.95)+#color="darkgrey") +
  theme_linedraw() + scale_fill_manual(values = color_dict) + 
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
  guides(fill=guide_legend("Superkingdom - Phylum" , ncol=1)) +
  xlab("Individual") + 
  ylab(paste("Read frequency")) + 
  theme(plot.title = element_text(hjust = 0.5, size=16), strip.text = element_text(size=18), 
        legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key.size =unit(.6,"cm"),
        axis.text.x=element_text(angle=30,size=8),
        panel.background = element_rect(fill = "#111111"),    
        panel.grid.major = element_line(colour = "grey",),
        panel.grid.minor = element_line(colour = "grey",)) +
  coord_flip() + 
  facet_wrap(~ Method, nrow = 1, scales = 'free_x' )






#### Change in colors: 
unimportant_reads=classif_smaller[classif_smaller$Method=="Blast",]
unimportant_reads=unimportant_reads[unimportant_reads$final_name %in% c("Unassembled","unclassified"),]

total_reads=unique(unimportant_reads[,c(".id","total"),])
t=aggregate(NumReads ~ .id, unimportant_reads, sum)
unimportant_reads=merge(t,total_reads, by=".id")
unimportant_reads$Freq=unimportant_reads$NumReads / unimportant_reads$total


unimportant_reads=classif_smaller[classif_smaller$Method=="K2uniq",]
unimportant_reads=unimportant_reads[unimportant_reads$final_name %in% c("cellular organisms","unclassified"),]

total_reads=unique(unimportant_reads[,c(".id","total"),])
t=aggregate(NumReads ~ .id, unimportant_reads, sum)
unimportant_reads=merge(t,total_reads, by=".id")
unimportant_reads$Freq=unimportant_reads$NumReads / unimportant_reads$total

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### barplots - Filtered



classif_red=classif_fullnames[!classif_fullnames$taxon_phylum %in% 
                                c("unclassified","artificial sequences","cellular organisms", "Unassembled", "Homo_sapiens",
                                  "Eukaryota-unclassified","Eukaryota-Chordata","Eukaryota-Arthropoda","Eukaryota-Streptophyta",
                                  "Eukaryota-Mollusca", "Eukaryota-Bryozoa") ,]

key_from_data=rename_taxa(classif_red, nbact = 4,nvir = 3,neuk = 3)
key_from_data$final_name = ifelse(str_detect(key_from_data$taxon_full_name,"myco"), 'Eukaryota-Fungi',key_from_data$final_name)


color_dict=Taxonomy_small_color_palette2(key_from_data, "final_name")
key_from_data$final_name[!key_from_data$final_name %in% names(color_dict)]

# change classification column 
classif_red=merge(classif_red, key_from_data, by.x = "taxon_phylum", by.y="taxon_full_name")
classif_red=aggregate(NumReads ~ .id + Method + final_name, classif_red, sum)

classif_red <- classif_red %>%
  group_by(.id,Method) %>%
  mutate(total = sum(NumReads)) %>%
  ungroup()  # Remove grouping for further operations
classif_red$Freq = classif_red$NumReads/classif_red$total

classify_reads=aggregate(NumReads ~ .id + Method, classif_red, sum)
median(classify_reads[classify_reads$Method=="Blast",]$NumReads)
median(classify_reads[classify_reads$Method=="K2uniq",]$NumReads)


tax=unique(classif_red$final_name)
tax_order=c(
  "artificial sequences" , "unclassified","cellular organisms", "Unassembled", "Homo_sapiens",
  tax[str_detect(tax,"Bacteria") & ! str_detect(tax, "Other")], "Other Bacteria", "Archaea",
  #tax[str_detect(tax,"myco") & ! str_detect(tax, "Other")], "Other Eukaryota - Fungi",
  tax[str_detect(tax,"Eukaryota") & ! str_detect(tax, "Other") & ! str_detect(tax,"myco")], "Other Eukaryota",
  tax[str_detect(tax,"Viruses") & ! str_detect(tax, "Other")], "Other Viruses")

classif_red$final_name <- factor(classif_red$final_name, levels=tax_order)

order_patient=rownames(n_reads)[order(n_reads[,1], decreasing = FALSE)]
classif_red$.id <- factor(classif_red$.id, levels=order_patient)


classif_red_blast=classif_red[classif_red$Method == "Blast" & classif_red$final_name == "Viruses-Kitrinoviricota", ]
blast_order=classif_red_blast[order(classif_red_blast$Freq) ,]$.id
classif_red$.id = factor(classif_red$.id, levels = blast_order)

p1=ggplot(classif_red, aes(x=.id, y=Freq, fill = final_name) ) +
  geom_bar(position="stack", stat="identity", width = 0.95)+#color="darkgrey") +
  theme_linedraw() + 
  scale_fill_manual(values = color_dict) + 
  guides(fill=guide_legend("Superkingdom - Phylum" , ncol=1)) +
  xlab("Individual") + 
  ylab(paste("Read frequency")) + 
  theme(plot.title = element_text(hjust = 0.5, size=16), strip.text = element_text(size=18), 
        legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key.size =unit(.6,"cm"),
        axis.text.x=element_text(angle=30,size=8),
        panel.background = element_rect(fill = "#111111"),    
        panel.grid.major = element_line(colour = "grey",),
        panel.grid.minor = element_line(colour = "grey",)) +
  coord_flip() + 
  facet_wrap(~ Method, nrow = 1, scales = 'free_x' )
p1


p2=ggplot(classif_red, aes(y=.id, x=NumReads, fill = final_name) ) +
  geom_bar(position="stack", stat="identity", width = 0.95)+#color="darkgrey") +
  theme_linedraw() + 
  scale_fill_manual(values = color_dict) + 
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
  guides(fill=guide_legend("Superkingdom - Phylum" , ncol=1)) +
  ylab("Individual") + 
  xlab(paste("Read Number")) + 
  theme(plot.title = element_text(hjust = 0.5, size=16), strip.text = element_text(size=18), 
        legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key.size =unit(.6,"cm"),
        axis.text.x=element_text(angle=30,size=8),
        panel.background = element_rect(fill = "#111111"),    
        panel.grid.major = element_line(colour = "grey",),
        panel.grid.minor = element_line(colour = "grey",)) +
  scale_x_continuous(labels = scales::comma) +
  facet_wrap(~ Method, nrow = 1)#+
#ggbreak::scale_x_break(c(1e6, 1e6), scales = 0.5)


p1
p2


counting2=aggregate(NumReads ~ Method+.id ,classif_red, sum)
counting3=aggregate(NumReads ~ Method ,counting2, mean)

ggsave(p1,filename=paste0("/shared/projects/microbiome_translocation/results/Douek_Cleveland/figures/distribution_proportion_Cleveland.png"),
       device = "png",width = 48, height = 27, units = "cm", create.dir = TRUE)
ggsave(p2,filename=paste0("/shared/projects/microbiome_translocation/results/Douek_Cleveland/figures/distribution_total_Cleveland.png"),
       device = "png",width = 48, height = 27, units = "cm")



ggplot(classif_red, aes(x=.id, y=NumReads, fill = final_name) ) +
  geom_bar(position="stack", stat="identity", width = 0.95)+#color="darkgrey") +
  theme_linedraw() + scale_fill_manual(values = color_dict) + 
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
  guides(fill=guide_legend("Superkingdom - Phylum" , ncol=1)) +
  xlab("Individual") + 
  ylab(paste("Read frequency")) + 
  theme(plot.title = element_text(hjust = 0.5, size=16), strip.text = element_text(size=18), 
        legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key.size =unit(.6,"cm"),
        axis.text.x=element_text(angle=30,size=8),
        panel.background = element_rect(fill = "#111111"),    
        panel.grid.major = element_line(colour = "grey",),
        panel.grid.minor = element_line(colour = "grey",)) +
  coord_flip(ylim = c(0,8.5e+5)) + 
  facet_wrap(~ Method, nrow = 1, scales = 'free_x' )


#### Composite figures: Retrieve reads predicted as Bacteria by either one of the two methods

classif_per_taxon=function(classif_list, taxon){
  ### Select reads that are classified as Bacteria, viruses by either one of the 2 methods
  classif_list=parallel::mclapply(classif_list,function(x){ 
    x=x[ (str_detect(x$taxon_K2Uniq, taxon) | str_detect(x$taxon_Blast, taxon)) ,  ]
    k2=format_read_table(x, "taxon_K2Uniq")  
    blast=format_read_table(x, "taxon_Blast")
    l=list("Kraken2"=k2, "Blast"=blast)
    l=ldply(l)
    colnames(l)[1]="Method"
    return(l)
  }, mc.cores=20)
  
  classif_df=ldply(classif_list)
  
  
  classif_df$phylum=ifelse(str_detect(classif_df$phylum, "unclassified") , "unclassified", classif_df$phylum)
  classif_df$taxon_phylum = ifelse(classif_df$phylum == "", 
                                   yes=classif_df$superkingdom,
                                   no = paste0(classif_df$superkingdom, "-",classif_df$phylum))
  
  
  classif_df$taxon_phylum = ifelse(classif_df$taxon_phylum %in% c("artificial sequences","cellular organisms"), 
                                   yes="unclassified",
                                   no = classif_df$taxon_phylum)
  
  
  # Collapse read number by the "new" phylum
  classif_df=aggregate(NumReads ~ .id + Method + taxon_phylum + superkingdom + phylum, classif_df, sum)
  
  
  if(taxon == "Bacteria"){
    key_from_data=rename_taxa(classif_df, nbact = 5,nvir = 3,neuk = 1)
  }else if(taxon == "Viruses") {
    key_from_data=rename_taxa(classif_df, nbact = 1,nvir = 3,neuk = 1)
  } else if(taxon == "Eukaryota") {
    key_from_data=rename_taxa(classif_df, nbact = 1,nvir = 3,neuk = 5)
  }else{
    return("error, please give a valid taxon (Eukaryota, Bacteria, Viruses)")
  }
  
  
  key_from_data$final_name = ifelse(str_detect(key_from_data$taxon_full_name,"myco") & str_detect(key_from_data$taxon_full_name,"Eukaryota")
                                    , 'Eukaryota-Fungi',key_from_data$final_name)
  
  key_from_data$final_name = ifelse(key_from_data$final_name == "Eukaryota-unclassified" , 
                                    "Other Eukaryota" , key_from_data$final_name )
  
  key_from_data$final_name = ifelse(key_from_data$final_name == "Bacteria-unclassified" , 
                                    "Other Bacteria" , key_from_data$final_name )
  
  color_dict=Taxonomy_small_color_palette2(key_from_data, "final_name")
  color_dict[["Unassembled"]]="#444444"
  color_dict[["Homo_sapiens"]]="darkorchid"
  color_dict[["unclassified"]]="#999999"
  
  
  
  if(taxon != "Eukaryota"){
    color_dict[["Eukaryota-Fungi"]]=colorRampPalette(rev(c("#026e02","darkolivegreen1")))(6)[3]
    color_dict[["Other Eukaryota"]]=colorRampPalette(rev(c("#026e02","darkolivegreen1")))(5)[5]
  }
  
  
  key_from_data$final_name[!key_from_data$final_name %in% names(color_dict)]
  # change classification column 
  classif_df=merge(classif_df, key_from_data, by.x = "taxon_phylum", by.y="taxon_full_name")
  classif_df=aggregate(NumReads ~ .id + Method + final_name, classif_df, sum)
  
  classif_df <- classif_df %>%
    group_by(.id,Method) %>%
    mutate(total = sum(NumReads)) %>%
    ungroup()  
  classif_df$Freq = classif_df$NumReads/classif_df$total
  
  bacteria_classif=unique(classif_df)
  
  tax=unique(classif_df$final_name)
  tax_order=c(
    "artificial sequences" , "unclassified","cellular organisms", "Unassembled", "Homo_sapiens",
    tax[str_detect(tax,"Bacteria") & ! str_detect(tax, "Other")], "Other Bacteria","Archaea",
    #tax[str_detect(tax,"myco") & ! str_detect(tax, "Other")], "Other Eukaryota - Fungi",
    tax[str_detect(tax,"Eukaryota") & ! str_detect(tax, "Other") & ! str_detect(tax,"myco")], "Other Eukaryota",
    tax[str_detect(tax,"Viruses") & ! str_detect(tax, "Other")], "Other Viruses")
  
  classif_df$final_name <- factor(classif_df$final_name, levels=tax_order)
  
  return(list("classif_df"=classif_df, "Keys"=key_from_data, "colordict"=color_dict))
}



# bacteria_classif_list=parallel::mclapply(classif,function(x){ 
#   x=x[ (str_detect(x$taxon_K2Uniq, "Bacteria") | str_detect(x$taxon_Blast, "Bacteria")) ,  ]
#   return(x)
#   
#   }, mc.cores=20)
# 
# bacteria_classif_list=parallel::mclapply(bacteria_classif_list,function(x){
#   k2=format_read_table(x, "taxon_K2Uniq")  
#   blast=format_read_table(x, "taxon_Blast")
#   l=list("K2uniq"=k2, "Blast"=blast)
#   l=ldply(l)
#   colnames(l)[1]="Method"
#   return(l)
# }, mc.cores=20)
# 
# bacteria_classif=ldply(bacteria_classif_list)
# 
# 
# bacteria_classif$phylum=ifelse(str_detect(bacteria_classif$phylum, "unclassified") , "unclassified", bacteria_classif$phylum)
# bacteria_classif$taxon_phylum = ifelse(bacteria_classif$phylum == "", 
#                                       yes=bacteria_classif$superkingdom,
#                                       no = paste0(bacteria_classif$superkingdom, "-",bacteria_classif$phylum))
# 
# 
# bacteria_classif$taxon_phylum = ifelse(bacteria_classif$taxon_phylum %in% c("artificial sequences","cellular organisms"), 
#                                        yes="unclassified",
#                                        no = bacteria_classif$taxon_phylum)
# 
# 
# # Collapse read number by the "new" phylum
# bacteria_classif=aggregate(NumReads ~ .id + Method + taxon_phylum + superkingdom + phylum, bacteria_classif, sum)
# key_from_data=rename_taxa(bacteria_classif, nbact = 5,nvir = 3,neuk = 1)
# key_from_data$final_name = ifelse(str_detect(key_from_data$taxon_full_name,"myco") & str_detect(key_from_data$taxon_full_name,"Eukaryota")
#                                   , 'Eukaryota-Fungi',key_from_data$final_name)
# 
# key_from_data$final_name = ifelse(key_from_data$final_name == "Eukaryota-unclassified" , 
#                                   "Other Eukaryota" , key_from_data$final_name )
# 
# key_from_data$final_name = ifelse(key_from_data$final_name == "Bacteria-unclassified" , 
#                                   "Other Bacteria" , key_from_data$final_name )
# 
# 
# # Creating color_palette:
# color_dict=Taxonomy_small_color_palette2(key_from_data, "final_name")
# key_from_data$final_name[!key_from_data$final_name %in% names(color_dict)]
# 
# 
# 
# 
# color_dict[["Unassembled"]]="#444444"
# color_dict[["Homo_sapiens"]]="darkorchid"
# color_dict[["unclassified"]]="#999999"
# 
color_dict[["Eukaryota-Fungi"]]=colorRampPalette(rev(c("#026e02","darkolivegreen1")))(6)[3]
color_dict[["Other Eukaryota"]]=colorRampPalette(rev(c("#026e02","darkolivegreen1")))(5)[5]
# 
# 
# 
# 
# # change classification column 
# bacteria_classif=merge(bacteria_classif, key_from_data, by.x = "taxon_phylum", by.y="taxon_full_name")
# bacteria_classif=aggregate(NumReads ~ .id + Method + final_name, bacteria_classif, sum)
# 
# ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# ### barplots
# 
# bacteria_classif_list=classif_per_taxon(classif, "Bacteria")
# virus_classif_list=classif_per_taxon(classif, "Viruses")
# euk_classif_list=classif_per_taxon(classif, "Eukaryota")
# 
# 
# bacteria_classif = bacteria_classif_list$classif_df
# ku=bacteria_classif[bacteria_classif$Method=="K2uniq",]
# 
# bacteria_classif <- bacteria_classif %>%
#   group_by(.id,Method) %>%
#   mutate(total = sum(NumReads)) %>%
#   ungroup()  
# bacteria_classif$Freq = bacteria_classif$NumReads/bacteria_classif$total
# 
# bacteria_classif=unique(bacteria_classif)
# 
# tax=unique(bacteria_classif$final_name)
# tax_order=c(
#   "artificial sequences" , "unclassified","cellular organisms", "Unassembled", "Homo_sapiens",
#   tax[str_detect(tax,"Bacteria") & ! str_detect(tax, "Other")], "Other Bacteria","Archaea",
#   #tax[str_detect(tax,"myco") & ! str_detect(tax, "Other")], "Other Eukaryota - Fungi",
#   tax[str_detect(tax,"Eukaryota") & ! str_detect(tax, "Other") & ! str_detect(tax,"myco")], "Other Eukaryota",
#   tax[str_detect(tax,"Viruses") & ! str_detect(tax, "Other")], "Other Viruses")
# 
# bacteria_classif$final_name <- factor(bacteria_classif$final_name, levels=tax_order)
# 
# color_dict
# 


bacteria_classif_list=classif_per_taxon(classif, "Bacteria")
virus_classif_list=classif_per_taxon(classif, "Viruses")
euk_classif_list=classif_per_taxon(classif, "Eukaryota")



p1=ggplot(bacteria_classif_list$classif_df,#[bacteria_classif$.id %in% c("SRR14418861", "SRR14418871"),], 
          aes(x= reorder (.id, total), y=NumReads, fill = final_name) ) +
  geom_bar(position="stack", stat="identity", width = 0.9)+#color="darkgrey") +
  theme_linedraw() + scale_fill_manual(values = bacteria_classif_list$colordict) + 
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
  guides(fill=guide_legend("Superkingdom - Phylum" , ncol=1)) +
  xlab("") + 
  ylab(paste("Read number")) + 
  theme(plot.title = element_text(hjust = 0.5, size=16), strip.text = element_text(size=18), 
        legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key.size =unit(.6,"cm"),
        axis.text.x=element_text(angle=0,size=8),
        axis.text.y=element_text(angle=0,size=6),
        panel.background = element_rect(fill = "#000000"),    
        panel.grid.major = element_line(colour = "grey",),
        panel.grid.minor = element_line(colour = "grey",)) +
  coord_flip() + 
  facet_wrap(~ Method, nrow = 1, scales = 'free_x' )
#facet_wrap(~ .id, scales = 'free_y' )


p2=ggplot(virus_classif_list$classif_df,#[bacteria_classif$.id %in% c("SRR14418861", "SRR14418871"),], 
          aes(x= reorder (.id, total), y=NumReads, fill = final_name) ) +
  geom_bar(position="stack", stat="identity", width = 0.9)+#color="darkgrey") +
  theme_linedraw() + scale_fill_manual(values = virus_classif_list$colordict) + 
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
  guides(fill=guide_legend("Superkingdom - Phylum" , ncol=1)) +
  xlab("") + 
  ylab(paste("Read number")) + 
  theme(plot.title = element_text(hjust = 0.5, size=16), strip.text = element_text(size=18), 
        legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key.size =unit(.6,"cm"),
        axis.text.x=element_text(angle=0,size=8),
        axis.text.y=element_text(angle=0,size=6),
        panel.background = element_rect(fill = "#000000"),    
        panel.grid.major = element_line(colour = "grey",),
        panel.grid.minor = element_line(colour = "grey",)) +
  #coord_cartesian(xlim = c(0,5e+5))+
  coord_flip() + 
  facet_wrap(~ Method, nrow = 1, scales = 'free_x' )
p2


p3=ggplot(euk_classif_list$classif_df,#[bacteria_classif$.id %in% c("SRR14418861", "SRR14418871"),], 
          aes(x= reorder (.id, total), y=NumReads, fill = final_name) ) +
  geom_bar(position="stack", stat="identity", width = 0.9)+#color="darkgrey") +
  theme_linedraw() + scale_fill_manual(values = euk_classif_list$colordict) + 
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
  guides(fill=guide_legend("Superkingdom - Phylum" , ncol=1)) +
  xlab("") + 
  ylab(paste("Read number")) + 
  theme(plot.title = element_text(hjust = 0.5, size=16), strip.text = element_text(size=18), 
        legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key.size =unit(.6,"cm"),
        axis.text.x=element_text(angle=0,size=8),
        axis.text.y=element_text(angle=0,size=6),
        panel.background = element_rect(fill = "#000000"),    
        panel.grid.major = element_line(colour = "grey",),
        panel.grid.minor = element_line(colour = "grey",)) +
  coord_flip() + 
  facet_wrap(~ Method, nrow = 1, scales = 'free_x' )

p1
p2
p3





ggsave(p1,filename=paste0("/shared/projects/microbiome_translocation/results/Douek_Cleveland/figures/Cleveland_bact.jpeg"),
       device = "jpeg",width = 60, height = 18,dpi = 700 , units = "cm", create.dir = TRUE)
ggsave(p2,filename=paste0("/shared/projects/microbiome_translocation/results/Douek_Cleveland/figures/Cleveland_vir.jpeg"),
       device = "jpeg",width = 60, height = 18,dpi = 700 , units = "cm", create.dir = TRUE)
ggsave(p3,filename=paste0("/shared/projects/microbiome_translocation/results/Douek_Cleveland/figures/Cleveland_euk.jpeg"),
       device = "jpeg",width = 45, height = 13.5,dpi = 700 , units = "cm", create.dir = TRUE)




# ggplot(bacteria_classif,#[bacteria_classif$.id %in% c("SRR14418861", "SRR14418871"),], 
#        aes(x=Method, y=NumReads, fill = final_name) ) +
#   geom_bar(position="stack", stat="identity", width = 0.95)+#color="darkgrey") +
#   theme_linedraw() + scale_fill_manual(values = color_dict) + 
#   theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
#   guides(fill=guide_legend("Superkingdom - Phylum" , ncol=1)) +
#   xlab("Individual") + 
#   ylab(paste("Read frequency")) + 
#   theme(plot.title = element_text(hjust = 0.5, size=16), strip.text = element_text(size=18), 
#         legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key.size =unit(.6,"cm"),
#         axis.text.x=element_text(angle=30,size=8),
#         panel.background = element_rect(fill = "#111111"),    
#         panel.grid.major = element_line(colour = "grey",),
#         panel.grid.minor = element_line(colour = "grey",)) +
#   #coord_flip() + 
#   #facet_wrap(~ .id, nrow = 1, scales = 'free_x' )
#   facet_wrap(~ .id, scales = 'free_y' )
#facet_grid(~ .id, scales = 'free', independent = "y" )









### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### %age agreement between methods


to_remove= c("unclassified","artificial sequences","cellular organisms", "Unassembled", "Homo_sapiens",
             "Eukaryota-unclassified","Eukaryota-Chordata","Eukaryota-Arthropoda","Eukaryota-Streptophyta",
             "Eukaryota-Mollusca", "Eukaryota-Bryozoa") 

replace_phylum_names=function(df,merged_column){
  
  df[,c("sk","p")]=str_split_fixed(df[[merged_column]], fixed("|"), n=2 )
  df$p=ifelse(str_detect(df$p, "unclassified"), "unclassified", df$p)
  df[[merged_column]] = paste0(df$sk, "-", df$p)
  df[[merged_column]] = ifelse(df[[merged_column]]=="unclassified-unclassified", "unclassified", df[[merged_column]] )
  df[[merged_column]] = ifelse(df$p == "", df$sk , df[[merged_column]] )
  df[,c("sk","p")]=NULL
  return(df[[merged_column]] ) }

classif = lapply(classif, function(x) {
  x$read = NULL
  return(x) })

classif_aggreement = lapply(classif, function(x){
  # Format columns
  x$taxon_K2Uniq = replace_phylum_names(x, "taxon_K2Uniq")
  x$taxon_Blast = replace_phylum_names(x, "taxon_Blast")
  
  # soft remove uninteresting reads
  x$taxon_K2Uniq = ifelse(x$taxon_K2Uniq %in% to_remove, "k2_unclassified", x$taxon_K2Uniq )
  x$taxon_Blast = ifelse(x$taxon_Blast %in% to_remove, "blast_unclassified", x$taxon_Blast )
  
  x$score = ifelse(x$taxon_K2Uniq == x$taxon_Blast, 1, 0)
  
  
  return(x)
})


# Counting agree score + dividing it with total classified reads
classif_score=ldply(lapply(classif_aggreement, function(x) return(sum(x$score)) )  )
colnames(classif_score)[2] = "agree_score"

classify_reads=aggregate(NumReads ~ .id + Method, classif_red, sum)
classify_reads_K2 = classify_reads[classify_reads$Method == "K2uniq",]
classify_reads_blast = classify_reads[classify_reads$Method == "Blast",]




classif_Blast = merge(classify_reads_blast, classif_score, by=".id")
classif_K2    = merge(classify_reads_K2, classif_score, by=".id")
classif_score = rbind(classif_Blast, classif_K2)

classif_score$agree_prop=classif_score$agree_score / classif_score$NumReads 

order_p=classif_score[order(classif_score[classif_score$Method=="Blast",]$agree_prop),]$.id
classif_score$.id = factor(classif_score$.id, levels = order_p)

ggplot(classif_score, aes(y=agree_prop, x=.id, fill = Method, color=Method)) + 
  geom_bar(stat="identity", position = "dodge" , width = 0.55) + 
  scale_fill_manual(values=c("#1BA1E2", "#E51400")) +
  scale_color_manual(values=c("#1BA1E2", "#E51400")) + 
  theme_minimal_hgrid() + 
  xlab("") + ylab("Proportion of reads where both methods agrees") + 
  theme(plot.title = element_text(hjust = 0.5, size=16), strip.text = element_text(size=18), 
        legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key.size =unit(.6,"cm"),
        axis.text.x=element_text(angle=30,size=10, hjust = 1, vjust=1)) +
  ggtitle("Proportion of reads where both method agrees at phylum level depending on the amount of non-human-classified reads by each methods")



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
## Online: get contig length and classif


distant_path="/run/user/1001/gvfs/sftp:host=core.cluster.france-bioinformatique.fr,user=acolajanni/shared/projects/microbiome_translocation/"
distant_path="/shared/projects/microbiome_translocation/"

contig_length=lapply(sra_list, function(x){
  print(x)
  t=read.delim(paste0(distant_path, "results/Douek_Cleveland/Contigs/",x,"/Contigs_name_length.txt"), 
               header = FALSE)
  
  t$contig_name <- sub("(TRINITY[^ ]+).*", "\\1", t$V1)
  t$length <- sub(".*len=([0-9]+)\\s.*", "\\1", t$V1)
  t$V1 = NULL
  
  
  contig_classif=read.delim(paste0(distant_path, "results/Douek_Cleveland/Contigs/",x,"/ContigsToReads/contigs_classif.txt"),
                            header = FALSE, col.names = c("contig_name","taxid","superkingdom"))
  contig_classif$contig_name = str_remove(contig_classif$contig_name, " ")
  
  merged=merge(t, contig_classif, by.x="contig_name", by.y="contig_name")
  merged$length = as.numeric(merged$length)
  
  return(merged) })

names(contig_length) = sra_list
save(contig_length, file = "/home/acolajanni/Documents/work/results/Douek_Cleveland/Contigs/contigs_length_classif.RData")

length_df=do.call(rbind, contig_length)
length_df$superkingdom = str_remove(length_df$superkingdom, " ")

n_contigs_per_sk=aggregate(length ~ superkingdom, length_df, length)
median_contigs=aggregate(length ~ superkingdom, length_df, median)

n_contigs_per_sk$label = paste0("n = ",n_contigs_per_sk$length)

length_df$superkingdom <- factor(length_df$superkingdom , levels=c("Bacteria", "Homo_sapiens", "unclassified", "Eukaryota", "Viruses"))

ggplot(length_df, aes(x=superkingdom, y=length, fill=superkingdom)) + 
  geom_hline(yintercept = 200) +
  geom_text(data=n_contigs_per_sk, aes(x=superkingdom, label=label), size=8, y=150) +
  #geom_violin(alpha=.5) + 
  geom_boxplot(width=.55, outlier.shape = NA)  + 
  theme_bw() +
  ylab("Contig length") + xlab(" \n Classification") + 
  theme(plot.title = element_text(hjust = 0.5, size=16), strip.text = element_text(size=18), 
        legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key.size =unit(.6,"cm"),
        axis.text.x=element_text(size=16, family = "bold", color="black")) + 
  scale_y_continuous(breaks=c(0, 200, 257, 280, 363  ,434, 500, 1000))+
  coord_cartesian(ylim = c(100,1000))+
  scale_fill_manual(values=c( "#A7C9D6", "purple", "#777777","#FF7597","#EEEE00")) + 
  ggtitle("Contig length distribution across all patients regroup by classification status")

library(ggpubr)

# my_comparisons <- list( c("Bacteria", "Homo_sapiens"), c("Bacteria", "unclassified"), c("Bacteria", "Eukaryota"),c("Bacteria", "Viruses"), 
#                         c("Homo_sapiens", "unclassified"), c("Homo_sapiens", "Eukaryota"),c("Homo_sapiens", "Viruses") ,
#                         c("unclassified", "Eukaryota"), c("unclassified", "Viruses"),
#                         c("Eukaryota", "Viruses"))
# 
# ggboxplot(length_df, x = "superkingdom", y = "length",
#           color = "superkingdom", palette = "jco")+ 
#   stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
#   stat_compare_means(label.y = 50)     # Add global p-value


kruskal_test <- kruskal.test(length ~ superkingdom, data = length_df)

pairwise_test <- pairwise.wilcox.test(length_df$length, length_df$superkingdom, p.adjust.method = "BH")
pairwise_test

# View the result



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
## Searching specific patient
library(cowplot)

p=classif$SRR14418878

pk=p[,c("read","taxon_K2Uniq")]
pb=p[,c("read","taxon_Blast")]
# pk$Method="K2Uniq"
# pb$Method="Blast"
# colnames(pb)[2]="taxon_phylum"
# colnames(pk)[2]="taxon_phylum"

p=list("K2Uniq"=pk, "Blast"=pb)

rm(pb)  
rm(pk)

p=lapply(p, function(x) {
  x[,c("superkingdom","phylum")]=str_split_fixed(x[,2], fixed("|"), n=2)
  x$phylum=ifelse(str_detect(x$phylum, "unclassified") , "unclassified", x$phylum)
  x[,2] = ifelse(x$phylum == "", 
                 yes=x$superkingdom,
                 no = paste0(x$superkingdom, "-",x$phylum))
  return(x)
})

count=as.data.frame(table(c(p$K2Uniq$taxon_K2Uniq, p$Blast$taxon_Blast)))
colnames(count)=c("taxon_phylum","NumReads")
key_from_data=rename_taxa(count, nbact = 5,nvir = 3,neuk = 5)

p=lapply(names(p), function(name) {
  if(name == "K2Uniq") {colname="taxon_K2Uniq"}
  else {colname="taxon_Blast"}
  print(name)
  x=p[[name]]
  
  x=merge(x,key_from_data, by.x= colname, by.y="taxon_full_name")
  x[[colname]]=NULL
  colnames(x)[colnames(x) == "final_name"] = colname 
  return(x)
})
names(p) = c("K2Uniq","Blast")

classif_merged=merge(p$K2Uniq,p$Blast, by="read")


color_dict=Taxonomy_small_color_palette(key_from_data, "final_name")
key_from_data$final_name[!key_from_data$final_name %in% names(color_dict)]


color_dict[["Not in Contigs"]]="red4"
color_dict[["Homo_sapiens"]]="purple"
color_dict[["unclassified"]]="#777777"
color_dict[["cellular organisms"]]="#444444"
color_dict[["artificial sequences"]]="#AAAAAA"





legend_df=data.frame("taxon"=unique(as.character(c(classif_merged$taxon_Blast, classif_merged$taxon_K2Uniq))))
legend_df$Freq=1

tax=unique(legend_df$taxon)
tax_order=c(
  "artificial sequences" , "unclassified","cellular organisms", "Not in Contigs", "Homo_sapiens",
  tax[str_detect(tax,"Bacteria") & ! str_detect(tax, "Other")], "Other Bacteria",
  tax[str_detect(tax,"myco") & ! str_detect(tax, "Other")], "Other Eukaryota - Fungi",
  tax[str_detect(tax,"Eukaryota") & ! str_detect(tax, "Other") & ! str_detect(tax,"myco")], "Other Eukaryota",
  "Archaea",
  tax[str_detect(tax,"Viruses") & ! str_detect(tax, "Other")], "Other Viruses")

tax_order

legend_df$taxon=factor(legend_df$taxon, levels=tax_order)



df <- classif_merged %>%
  make_long(taxon_Blast, taxon_K2Uniq,)

p1=ggplot(df, aes(x = x, 
                  next_x = next_x, 
                  node = node, 
                  next_node = next_node,
                  fill = factor(node) )) + 
  geom_sankey(flow.alpha=0.85,node.color=1,space=300,show.legend = FALSE ) +
  scale_fill_manual(values = color_dict) + 
  guides(fill = guide_legend(ncol = 1))+
  xlab("SRR14418878")+#  scale_x_discrete(labels = c("Blue pipeline","Red pipeline")) +
  theme_sankey(base_size = 16, ) 

plot=ggplot(legend_df, aes(x=taxon, y=Freq, fill=taxon))+
  geom_bar(stat = "identity") + scale_fill_manual(values=color_dict)+
  theme(plot.title = element_text(hjust = 0.5, size=16), strip.text = element_text(size=18), 
        legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key.size =unit(.6,"cm"),
        axis.text.x=element_text(angle=30,size=8),
        panel.background = element_rect(fill = "#111111"),    
        panel.grid.major = element_line(colour = "grey",),
        panel.grid.minor = element_line(colour = "grey",),
        legend.position="right")+
  guides(fill=guide_legend("Superkingdom - Phylum" , ncol=1)) 
plot

legend=cowplot::get_plot_component(plot,'guide-box-right',return_all = TRUE)
cowplot::plot_grid(p1,legend , rel_widths = c(.9,.25) )

classif_merged$transition=paste0(classif_merged$taxon_Blast, " --> ", classif_merged$taxon_K2Uniq)
classif_merged_pisu=classif_merged[str_detect(classif_merged$transition, "Viruses"),]



df2 <- classif_merged_pisu %>%
  make_long(taxon_Blast, taxon_K2Uniq,)

p2=ggplot(df2, aes(x = x, 
                   next_x = next_x, 
                   node = node, 
                   next_node = next_node,
                   fill = factor(node) )) + 
  geom_sankey(flow.alpha=0.85,node.color=1,space=300,show.legend = FALSE ) +
  scale_fill_manual(values = color_dict) + 
  guides(fill = guide_legend(ncol = 1))+
  xlab("SRR14418878")+#  scale_x_discrete(labels = c("Blue pipeline","Red pipeline")) +
  theme_sankey(base_size = 16 ) 

p2=cowplot::plot_grid(p2,legend , rel_widths = c(.9,.25) )

ggsave(p2,filename=paste0("/home/acolajanni/Documents/work/results/Douek_Cleveland/sankey_plots/SRR14418878_viruses.png"),
       device = "png",height = 24, width = 40, units = "cm")



################################################################################
# Global sankey plot


classif_merged=lapply(classif, function(x){ 
  x$read=NULL
  return(x) }) 
classif_merged=do.call(rbind, classif_merged)


rm(classif)
rm(df)

df3 <- classif_merged %>%
  make_long(taxon_Blast, taxon_K2Uniq)

p3=ggplot(df3, aes(x = x, 
                   next_x = next_x, 
                   node = node, 
                   next_node = next_node,
                   fill = factor(node) )) + 
  geom_sankey(flow.alpha=0.85,node.color=1,space=300,show.legend = FALSE ) +
  scale_fill_manual(values = color_dict) + 
  guides(fill = guide_legend(ncol = 1))+
  xlab("Sanley plot aggregated on all patients")+#  scale_x_discrete(labels = c("Blue pipeline","Red pipeline")) +
  theme_sankey(base_size = 16 ) 
################################################################################


ranks=c("superkingdom","phylum", "class", "order", "family")

# load(paste0(path,"results/Douek_Cleveland/kraken/nt/global_ReadsClassif.Rdata"))  

kraken_dir=paste0(path,"results/Douek_Cleveland/kraken/nt/k2uniq/")
ku_dir=paste0(path,"results/Douek_Cleveland/kraken/microbialDB/kuniq/")
ku_hybrid_dir=paste0(path,"results/Douek_Cleveland/hybrid/hybrid_Blast-kuniq/")

classif_ku=parallel::mclapply(sra_list, function(sra){
  x=read.table(paste0(ku_hybrid_dir,sra,"/ReadsClassif.txt"), header=FALSE, fill = TRUE, sep = "\t",quote = "")
  #colnames(x)=c("read", "taxid", "superkingdom","phylum", "class", "order", "family", "genus", "species")
  colnames(x)=c("read", "superkingdom","phylum", "class", "order", "family", "genus", "species")
  x=x[,! colnames(x) %in% c("genus", "species")]
  x[x$taxid %in% c(1,198431,155900,131567), ranks] = "unclassified"
  x[x==""]="unclassified"
  x=x[x$family != "unclassified" ,]
  return(x)
},mc.cores = 20)
names(classif_ku) = sra_list
# classif_ku=ldply(classif_ku)
# classif_ku$taxid=NULL
# classif_ku_agg=aggregate(read ~ . , classif_ku, length)

classification=lapply(classif_ku, function(x){
  x$read=NULL
  x$taxid=NULL
  x=unique(x)
  return(x)
})

classification=ldply(classification)
classification$.id=NULL
classification=unique(classification)
write.table(classification, file="/shared/projects/microbiome_translocation/results/Douek_Cleveland/kraken/nt/classification.csv",
            quote=FALSE, row.names=FALSE, sep=",")


classif_ku2=parallel::mclapply(classif_ku, function(x){
  x$classification <- apply(x[, ranks], 1, paste0, collapse = "|")
  t=aggregate(superkingdom ~ classification, x, length )
  colnames(t) = c("classification","freq")
  t$freq = t$freq / sum(t$freq)
  return(t)
},mc.cores = 10)

names(classif_ku2) = sra_list
classif_list=ldply(classif_ku2)

test=classif_list %>% tidyr::pivot_wider(names_from = classification, values_from = freq)
test=as.data.frame(t(test))
test$family = row.names(test)
colnames(test) = test[1,]
test=test[2:nrow(test),]
test[is.na(test)] = 0
test=test[,c(".id",sra_list)]
colnames(test)[1]="classification"


metadata=read.csv("/shared/projects/microbiome_translocation/data/Douek_Cleveland/Metadata.csv")
metadata$tmp=sapply(metadata$subject_classification, function(x) str_split(x, fixed(" "))[[1]][1] )
metadata$patient_name=paste0(metadata$tmp , "|", metadata$Run)


# Create a named vector for mapping
column_mapping <- setNames(metadata$patient_name, metadata$Run)

# Rename columns in df2
colnames(test) <- column_mapping[colnames(test)]
colnames(test)[1]="classification"

test[, -1] <- lapply(test[, -1], as.numeric)


write.table(test, file="/shared/projects/microbiome_translocation/results/Douek_Cleveland/hybrid/hybrid_Blast-kuniq/Hybrid_results.csv",
            quote=FALSE, row.names=FALSE, sep=",")


family=test$family
t=as.data.frame(table(family))
t2=as.data.frame(table(classification$family))
redundant=t2[t2$Freq>1,]$Var1


classification$classif <- apply(classification[, ranks], 1, paste0, collapse = "|")





classification=unique(classification)

