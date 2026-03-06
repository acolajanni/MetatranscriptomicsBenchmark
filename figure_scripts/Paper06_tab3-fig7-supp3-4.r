#!/shared/ifbstor1/software/miniconda/envs/r-4.2.3/bin/Rscript
#install.packages("ggh4x")
#devtools::install_github("davidsjoberg/ggsankey")

library(stringr)
library(plyr)
library(ggsankey)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggbreak) 
library(patchwork)
library(tidyr)
library(viridis)

label_df <- as.data.frame(
  list(
    Method = c(
      # Contig solo
      "trinityBlast","spadesBlast" ,
      # Hybrid Trinity
      #"hybrid_Blast-kraken2", 
      "trinityKu", 
      # Hybrid spades
      #"hybrid_SpadesBlast-kraken2",
      "spadesKu",
      # double kmer
      #"hybrid_kraken2-kuniq", 
      "kuk2", 
      # kmer solo
      "kraken2", "kuniq" ),
    
    labels2lines = c(
      # Contig solo
      "- Trinity-Blast -\nnt",
      "- Spades-Blast -\nnt",
      # Hybrid Trinity
      #"Hybrid-Trinity-K2\nnt",
      "- Hybrid-Trinity-KUniq -\nnt + microbialDB",
      # Hybrid spades
      #"Hybrid-Spades-K2\nnt",
      "- Hybrid-Spades-KUniq -\nnt + microbialDB",
      # double kmer
      #"K2-KUniq\nnt+microbialDB",
      "- KUniq-K2 -\nmicrobialDB+nt",
      # kmer solo
      "- K2 -\nnt",
      "- KUniq -\nmicrobialDB" ),
    
    family_method = c(
      "Assembly-based","Assembly-based",
      "Hybrid", "Hybrid",
      #"Hybrid", "Hybrid",
      "Combined assembly-free",#"Combined assembly-free",
      "Assembly-free","Assembly-free" )
    

  ))



label_df$family_method = factor(label_df$family_method, levels = c("Assembly-based", "Assembly-free" ,"Hybrid", "Combined assembly-free") )

label_df$Method = factor(label_df$Method, levels = c("trinityBlast" , "spadesBlast" , 
                                                     "kraken2" , "kuniq", 
                                                     #"hybrid_SpadesBlast-kraken2" , 
                                                     #"Hybrid - Blast + Kraken2" ,
                                                     "spadesKu" , 
                                                     "trinityKu" , 
                                                     #"Hybrid - Kraken2 + KrakenUniq",
                                                     "kuk2") )







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
  
  
  # if (replace_taxon == "Eukaryota"){
  #   
  #   print("eukaryota")
  #   fungi_df=table_df[str_detect(table_df[[column_to_filter]], "Fungi|fungi|myco" ),]
  #   fungi_df$superkingdom = "Eukaryota - Fungi"
  #   
  #   fun=aggregate_taxa_count_df(fungi_df, 4, column_to_filter, Frequency_column, "Eukaryota - Fungi"  )
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
  # 
  
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
path="~/"
print(dataset)
ncores=28

path_result=paste0(path,"results/",dataset)
path_data=paste0(path,"data/",dataset)

sra_list=readLines(paste0(path_data,"/sra_list_RNA.txt"))
sra_metadata=read.table(paste0(path_data,"/sra_metadata.tsv"),header=TRUE)


load(paste0(path_result,"/Taxon_correspondance.Rdata"))
load(paste0(path_result,"/color_dict.Rdata"))


classif=parallel::mclapply(sra_list, function(sra){  
  path_rdata=paste0(path_result,"/rdata/reviews/")
  print(sra)
  load(paste0(path_rdata, sra, "_classif_merged.rdata"))  
  return(classif_merged)  
  #rm(classif_merged)
},mc.cores=ncores)
names(classif) = sra_list



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Distribution of the number of sequences remaining after filtration

n_reads=t(as.data.frame(lapply(classif, function(x) nrow(x))))
colnames(n_reads)= "read_number"

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
theme_fig6=function(plot){
  return(plot + theme(strip.background = element_rect(fill="#333333") , 
                      strip.text = element_text(size = 7, face = "bold"),
                      axis.title.y= element_text(size = 8, face = "bold"),
                      axis.title.x= element_text(size = 8),
                      axis.text.x = element_text(size = 7),
                      axis.text.y = element_text(size = 7),
                      
                      legend.position = "right",
                      legend.key = element_rect(color="black"),
                      legend.text=element_text(size=8),
                      legend.title = element_blank(),
                      legend.key.width = unit(0.55, "cm"),
                      legend.key.height = unit(0.55, "cm"),  
                      
                      panel.grid.major = element_line(colour = "#454545"),
                      panel.grid.minor = element_line(colour = "#454545") ))
}
### counting where both methods agrees

normalize_labels <- function(x) {
  x <- tolower(x)
  x[x %in% c("unassembled", "unclassified")] <- "unclassified"
  x
}

compute_agreement <- function(df) {
  # remove 'read' column if present
  df2 <- df[, setdiff(names(df), "read"), drop = FALSE]
  
  # get all column pairs
  cols <- names(df2)
  pairs <- combn(cols, 2, simplify = FALSE)
  
  # compute agreement score for each pair
  scores <- sapply(pairs, function(p) {
    mean( normalize_labels(df2[[p[1]]]) == normalize_labels(df2[[p[2]]]),
          na.rm = TRUE )})
  
  # name the results
  names(scores) <- sapply(pairs, function(p) paste(p, collapse = " vs "))
  
  scores
}


compute_agreement_2 <- function(df) {
  df2 <- df[, setdiff(names(df), "read"), drop = FALSE]
  
  cols <- names(df2)
  pairs <- combn(cols, 2, simplify = FALSE)
  
  scores <- sapply(pairs, function(p) {
    a <- normalize_labels(df2[[p[1]]])
    b <- normalize_labels(df2[[p[2]]])
    
    # rows to exclude
    drop_idx <- (a == "unclassified" & b == "unclassified") |
      (a == "Homo_sapiens" & b == "Homo_sapiens")
    
    # keep only informative rows
    keep <- !drop_idx & !is.na(a) & !is.na(b)
    
    if (!any(keep)) return(NA_real_)
    
    mean(a[keep] == b[keep])
  })
  
  names(scores) <- sapply(pairs, function(p) paste(p, collapse = " vs "))
  scores
}

# apply to all elements of classif
#all_scores <- lapply(classif, compute_agreement)

scores_df <- parallel::mclapply(names(classif), function(nm) {
  s <- compute_agreement_2(classif[[nm]])
  data.frame(
    sample = nm,
    comparison = names(s),
    agreement = as.numeric(s),
    row.names = NULL
  )
},mc.cores=ncores) |> bind_rows()

#### Square matrix of agreement


score_heatmap = aggregate(agreement ~ comparison, scores_df, mean)
score_heatmap[, c("method1", "method2")] <- do.call(rbind, strsplit(score_heatmap$comparison, " vs ", fixed = TRUE) )
score_heatmap$comparison = NULL
score_heatmap2=score_heatmap

score_heatmap2$method1 = score_heatmap$method2
score_heatmap2$method2 = score_heatmap$method1

score_heatmap=rbind(score_heatmap2, score_heatmap)

score_triangle=score_heatmap

score_triangle <- score_heatmap %>%
  bind_rows(
    tibble(
      method1 = unique(c(score_heatmap$method1, score_heatmap$method2)),
      method2 = unique(c(score_heatmap$method1, score_heatmap$method2)),
      agreement = 1
    ) %>%
      filter(method1 == method2)
  )


order_labels=c("- Trinity-Blast -\nnt",
  "- Spades-Blast -\nnt",
  "- K2 -\nnt",
  "- KUniq -\nmicrobialDB",
  "- Hybrid-Trinity-KUniq -\nnt + microbialDB",
  "- Hybrid-Spades-KUniq -\nnt + microbialDB",
  "- KUniq-K2 -\nmicrobialDB+nt")

order_methods=c("trinityBlast" , "spadesBlast" , 
  "kraken2" , "kuniq", 
  "trinityKu" , 
  "spadesKu" , 
  "kuk2") 

score_triangle=merge(label_df, score_triangle, by.x="Method", by.y="method1")
score_triangle=merge(label_df, score_triangle, by.x="Method", by.y="method2")


score_triangle$Method.y = factor(score_triangle$Method.y, levels = order_methods)
score_triangle$Method = factor(score_triangle$Method, levels = order_methods)


score_triangle$labels2lines.x = factor(score_triangle$labels2lines.x, 
                                       levels = order_labels)
score_triangle$labels2lines.y = factor(score_triangle$labels2lines.y, 
                                       levels = order_labels)


p=ggplot(score_triangle, aes(labels2lines.x, labels2lines.y, fill = agreement)) +
  geom_tile(color = "white") +
  geom_text(aes(
    label = sprintf("%.3f", agreement),
    color = ifelse(agreement > 0.5, "white", "black")  # conditional text color
  ), size = 3, fontface="bold") +
  scale_color_manual(values=c("white","black"))+
  guides(color="none")+
  scale_fill_viridis(option = "F") +
  theme_linedraw() +
  labs(x = NULL,y = NULL, fill = "Agreement") + 
  ggh4x::facet_nested(cols=vars(family_method.x,labels2lines.x), scales="free")

#### Square matrix of agreement - with or without human/unclassified

scores_df2 <- parallel::mclapply(names(classif), function(nm) {
  s <- compute_agreement(classif[[nm]])
  data.frame(
    sample = nm,
    comparison = names(s),
    agreement = as.numeric(s),
    row.names = NULL
  )
},mc.cores=ncores) |> bind_rows()



score_heatmap_un = aggregate(agreement ~ comparison, scores_df2, mean)
score_heatmap_un[, c("method1", "method2")] <- do.call(rbind, strsplit(score_heatmap_un$comparison, " vs ", fixed = TRUE) )
score_heatmap_un$comparison = NULL
score_heatmap_un2=score_heatmap_un

score_heatmap_un2$method1 = score_heatmap_un$method2
score_heatmap_un2$method2 = score_heatmap_un$method1

score_heatmap_un=rbind(score_heatmap_un2, score_heatmap_un)

score_triangle_un <- score_heatmap_un %>%
  bind_rows(
    tibble(
      method1 = unique(c(score_heatmap_un$method1, score_heatmap_un$method2)),
      method2 = unique(c(score_heatmap_un$method1, score_heatmap_un$method2)),
      agreement = 1
    ) %>%
      filter(method1 == method2)
  )


order_labels=c("- Trinity-Blast -\nnt",
  "- Spades-Blast -\nnt",
  "- K2 -\nnt",
  "- KUniq -\nmicrobialDB",
  "- Hybrid-Trinity-KUniq -\nnt + microbialDB",
  "- Hybrid-Spades-KUniq -\nnt + microbialDB",
  "- KUniq-K2 -\nmicrobialDB+nt")

order_methods=c("trinityBlast" , "spadesBlast" , 
  "kraken2" , "kuniq", 
  "trinityKu" , 
  "spadesKu" , 
  "kuk2") 

score_triangle_un=merge(label_df, score_triangle_un, by.x="Method", by.y="method1")
score_triangle_un=merge(label_df, score_triangle_un, by.x="Method", by.y="method2")


score_triangle_un$Method.y = factor(score_triangle_un$Method.y, levels = order_methods)
score_triangle_un$Method = factor(score_triangle_un$Method, levels = order_methods)


score_triangle_un$labels2lines.x = factor(score_triangle_un$labels2lines.x, 
                                       levels = order_labels)
score_triangle_un$labels2lines.y = factor(score_triangle_un$labels2lines.y, 
                                       levels = order_labels)


p=ggplot(score_triangle_un, aes(labels2lines.x, labels2lines.y, fill = agreement)) +
  geom_tile(color = "white") +
  geom_text(aes(
    label = sprintf("%.3f", agreement),
    color = ifelse(agreement > 0.5, "white", "black")  # conditional text color
  ), size = 3, fontface="bold") +
  scale_color_manual(values=c("white","black"))+
  guides(color="none")+
  scale_fill_viridis(option = "F") +
  theme_linedraw() +
  labs(x = NULL,y = NULL, fill = "Agreement") + 
  ggh4x::facet_nested(cols=vars(family_method.x,labels2lines.x), scales="free")



score_triangle_un$metric="with_unclassified"
score_triangle$metric="without_unclassified"


score_triangle_full=rbind(score_triangle,score_triangle_un)
score_triangle_full <- score_triangle_full %>%
  mutate(
    ix = as.numeric(Method),
    iy = as.numeric(Method.y)
  )

upper_df <- score_triangle_full %>%
  filter(metric == "with_unclassified", ix <= iy)

lower_df <- score_triangle_full %>%
  filter(metric == "without_unclassified", ix >= iy)


plot_df <- bind_rows(upper_df, lower_df)

p=ggplot(plot_df, aes(labels2lines.x, labels2lines.y, fill = agreement)) +
  geom_tile(color = "white") +
  geom_text(aes(
    label = sprintf("%.3f", agreement),
    color = ifelse(agreement > 0.5, "white", "black")  # conditional text color
  ), size = 4, fontface="bold") +
  scale_color_manual(values=c("white","black"))+
  guides(color="none")+
  scale_fill_viridis(option = "F") +
  theme_linedraw() +
  labs(x = NULL,y = NULL, fill = "Agreement") + 
  ggh4x::facet_nested(cols=vars(family_method.x,labels2lines.x), scales="free")



##### Merge both:
heatmap_aggreement2=theme_fig6(p) + 
  labs(fill = "Pairwise agreement") +
  guides(fill = guide_colorbar(title.position = "top")) +
  theme(legend.position = "bottom",legend.title = element_text(hjust = 0.5) )

heatmap_aggreement2

figsupp4jpeg="~/results/Simulation/figures_review/figsupp4.jpeg"
figsupp4svg="~/results/Simulation/figures_review/figsupp4.svg"

ggsave(heatmap_aggreement2,
       filename = figsupp4jpeg,
       device = "jpeg",
       width = 29,
       height = 20,
       dpi = 600,
       units = "cm",
       create.dir = TRUE)

magick::image_read(figsupp4jpeg)

ggsave(heatmap_aggreement2,
       filename = figsupp4svg,
       device = "svg",
       width = 29,
       height = 20,
       dpi = 600,
       units = "cm",
       create.dir = TRUE)



################################################################################
## ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Proportions as in the poster


################################################################################ Tab 3

count_read=parallel::mclapply(classif, function(df){
  # remove 'read' column
  df_no_read <- df[, setdiff(names(df), "read")]
  
  # function to count categories per column
  count_categories <- function(col) {
    c(
      unclassified = sum(col %in% c("unclassified","Unassembled"), na.rm = TRUE),
      Unassembled  = sum(col == "Unassembled",  na.rm = TRUE),
      Homo_sapiens = sum(col == "Homo_sapiens", na.rm = TRUE),
      other        = sum(!(col %in% c("unclassified", 
                                      "Unassembled", 
                                      "Homo_sapiens")), 
                         na.rm = TRUE)
    )
  }
  
  # apply to all columns
  result <- as.data.frame(t(sapply(df_no_read, count_categories)))
  
  df_long <- result %>%
    as.data.frame() %>%
    mutate(method = rownames(.)) %>%
    pivot_longer(-method,
                 names_to = "type",
                 values_to = "n_seq")
  
  
  
  return(df_long)
},mc.cores=ncores)
count_read=ldply(count_read)

mean_read=aggregate(n_seq ~ method + type, count_read, mean)
sd_read=aggregate(n_seq ~ method + type, count_read, sd)
colnames(mean_read)[3] = "mean_read"
colnames(sd_read)[3] = "sd_read"
count_read=merge(mean_read,sd_read)

library(scales)

count_read <- count_read %>%
  mutate(
    mean_sd = paste0(
      comma(round(mean_read, 0)),
      " \u00B1 ",
      comma(round(sd_read, 0))))




# write.table(count_read[,c(1,2,5)], file = "~/results/Simulation/figures_review/table3_data_reviews.tsv",
#             quote=FALSE, row.names=FALSE, sep="\t")


################################################################################
### Trying rebuilding quant_lvl item: quantification for each level of classification
ranks=c("superkingdom","phylum", "class", "order", "family", "genus")
ranks=c("superkingdom","phylum")

### dataframe with: method, taxon/phylum, number of read

cols <- setdiff(names(classif$SRR14418854), "read")

classif_smaller_list <- parallel::mclapply(classif, function(x) {
  l <- lapply(cols, function(method) {
    tmp=format_read_table(x, method)
    
    tmp_clean <- tmp %>%
      mutate(
        is_srr = str_detect(taxon_phylum, "^SRR"),
        
        taxon_phylum = if_else(is_srr, "unclassified", taxon_phylum),
        superkingdom = if_else(is_srr, "unclassified", superkingdom),
        phylum = if_else(is_srr, "", phylum)
      ) %>%
      select(-is_srr) %>%
      group_by(taxon_phylum, superkingdom, phylum) %>%
      summarise(NumReads = sum(NumReads), .groups = "drop")
    
    return(tmp_clean)
    
  })

  names(l) <- cols
  l <- ldply(l)
  colnames(l)[1] <- "Method"
  return(l)
},mc.cores=ncores)
names(classif_smaller_list) = names(classif)


##############################################
## replace "unclassified bacteria phylum" 
## by "unclassified"
classif_smaller=ldply(classif_smaller_list)
total_df=aggregate(NumReads ~ .id + Method, classif_smaller, sum)

classif_smaller$phylum=ifelse(str_detect(classif_smaller$phylum, "unclassified") , "unclassified", classif_smaller$phylum)
classif_smaller$taxon_phylum = ifelse(classif_smaller$phylum == "", 
                                      yes=classif_smaller$superkingdom,
                                      no = paste0(classif_smaller$superkingdom, "-",classif_smaller$phylum))
# Collapse read number by the "new" phylum
classif_smaller=aggregate(NumReads ~ .id + Method + taxon_phylum + superkingdom + phylum, classif_smaller, sum)

classif_fullnames=classif_smaller
################################################################################
################################################################################
# 

key_from_data=rename_taxa(classif_smaller, nbact = 5,nvir = 4,neuk = 5)
key_from_data$final_name = ifelse(key_from_data$final_name %in% c("artificial sequences", "cellular organisms","unclassified"), 
                                  "Unclassified" , key_from_data$final_name)
key_from_data$final_name = ifelse(key_from_data$final_name == "Homo_sapiens", 
                                  "Homo sapiens" , key_from_data$final_name)
key_from_data=unique(key_from_data)
# Creating color_palette:
color_dict=Taxonomy_small_color_palette(key_from_data, "final_name")
key_from_data$final_name[!key_from_data$final_name %in% names(color_dict)]


color_dict[["Unassembled"]]="#444444"
color_dict[["Homo sapiens"]]="darkorchid"
color_dict[["Unclassified"]]="#999999"


setdiff( unique(classif_smaller$taxon_phylum), key_from_data$taxon_full_name )


# change classification column 
classif_smaller=merge(classif_smaller, key_from_data, by.x = "taxon_phylum", by.y="taxon_full_name")
classif_smaller=aggregate(NumReads ~ .id + Method + final_name, classif_smaller, sum)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### barplots

classif_smaller <- classif_smaller %>%
  group_by(.id,Method) %>%
  mutate(total = sum(NumReads)) %>%
  ungroup()  # Remove grouping for further operations
classif_smaller$Freq = classif_smaller$NumReads/classif_smaller$total

classif_smaller=unique(classif_smaller)

tax=unique(classif_smaller$final_name)

tax_order=c(
  "artificial sequences" , "Unclassified","cellular organisms", "Unassembled", "Homo sapiens",
  tax[str_detect(tax,"Bacteria") & ! str_detect(tax, "Other")], "Other Bacteria",
  tax[str_detect(tax,"Eukaryota") & ! str_detect(tax, "Other") ],#& ! str_detect(tax,"myco")], 
  "Other Eukaryota", "Archaea",
  tax[str_detect(tax,"Viruses") & ! str_detect(tax, "Other")], "Other Viruses")


setdiff(tax,tax_order)
tax[!tax %in% tax_order]

classif_smaller$final_name <- factor(classif_smaller$final_name, levels=tax_order)


order_patient=rownames(n_reads)[order(n_reads[,1], decreasing = FALSE)]
classif_smaller$.id <- factor(classif_smaller$.id, levels=order_patient)

classif_smaller=merge(label_df, classif_smaller, by="Method")



classif_smaller$labels2lines = factor(classif_smaller$labels2lines, 
                                       levels = order_labels)
classif_smaller$family_method = ifelse(classif_smaller$family_method == "Combined assembly-free",
                                       "Combined\nassembly-free",as.character(classif_smaller$family_method))
classif_smaller$family_method=factor(classif_smaller$family_method , 
                                     levels = c("Assembly-based","Assembly-free","Hybrid","Combined\nassembly-free"))

p=ggplot(classif_smaller, aes(x=.id, y=NumReads, fill = final_name) ) +
  geom_bar(position="stack", stat="identity", width = 0.95)+#color="darkgrey") +
  theme_linedraw() + scale_fill_manual(values = color_dict) + 
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
  guides(fill=guide_legend("Superkingdom - Phylum" , ncol=1)) +
  xlab("") + 
  ylab(paste("Read quantity")) +
  coord_flip() + 
  ggh4x::facet_nested(cols = vars(family_method, labels2lines), 
               scales="free")
  


fig6supp=theme_fig6(p) + theme(strip.text = element_text(size=7))
figsupp3="~/results/Simulation/figures_review/figsupp3.jpeg"
figsupp3svg="~/results/Simulation/figures_review/figsupp3.svg"

ggsave(fig6supp,
       filename = figsupp3,
       device = "jpeg",
       width = 30,
       height = 20,
       dpi = 300,
       units = "cm",
       create.dir = TRUE)

magick::image_read(figsupp3)
ggsave(fig6supp,
       filename = figsupp3svg,
       device = "svg",
       width = 30,
       height = 20,
       dpi = 600,
       units = "cm",
       create.dir = TRUE)


unique(classif_smaller$Method)


classif_smaller$family_method = ifelse(classif_smaller$family_method == "Combined\nassembly-free",
                                       "Combined assembly-free",as.character(classif_smaller$family_method))
classif_smaller$family_method=factor(classif_smaller$family_method , 
                                     levels = c("Assembly-based","Assembly-free","Hybrid","Combined assembly-free"))

p=ggplot(classif_smaller[classif_smaller$Method %in% c("kuk2","trinityBlast","trinityKu"),], aes(x=.id, y=NumReads, fill = final_name) ) +
  geom_bar(position="stack", stat="identity", width = 0.95)+#color="darkgrey") +
  theme_linedraw() + scale_fill_manual(values = color_dict) + 
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
  guides(fill=guide_legend("Superkingdom - Phylum" , ncol=1)) +
  xlab("") + 
  ylab(paste("Read quantity")) +
  coord_flip() + 
  ggh4x::facet_nested(cols = vars(family_method, labels2lines), 
                      scales="free")



fig6=theme_fig6(p)  + 
  theme(strip.text = element_text(size=9),
        legend.text = element_text(size=10),
        axis.text.y = element_text(size=8))

fig7jpg="~/results/Simulation/figures_review/fig7.jpeg"
fig7svg="~/results/Simulation/figures_review/fig7.svg"

ggsave(fig6,
       filename = fig7jpg,
       device = "jpeg",
       width = 30,
       height = 20,
       dpi = 300,
       units = "cm",
       create.dir = TRUE)

magick::image_read(fig6jpg)

ggsave(fig6,
       filename = fig7svg,
       device = "svg",
       width = 30,
       height = 20,
       dpi = 600,
       units = "cm",
       create.dir = TRUE)



save(color_dict,file = paste0("~/results/Douek_Cleveland/color_dict_fig6.rdata"))


p=ggplot(classif_smaller, aes(x=.id, y=Freq, fill = final_name) ) +
  geom_bar(position="stack", stat="identity", width = 0.95)+#color="darkgrey") +
  theme_linedraw() + scale_fill_manual(values = color_dict) + 
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
  guides(fill=guide_legend("Superkingdom - Phylum" , ncol=1)) +
  xlab("") + 
  ylab(paste("Read quantity")) +
  coord_flip() + 
  facet_nested(cols = vars(family_method, labels2lines), 
               scales="free")

fig6_alt=theme_fig6(p)  
fig6_alt

counting=aggregate(NumReads ~ Method+final_name, classif_smaller, sum)
### Mean number of reads remaining after removing unassembled sequences


#### Change in colors: 
unimportant_reads=classif_smaller[classif_smaller$Method=="trinityBlast",]
unimportant_reads=unimportant_reads[unimportant_reads$final_name %in% c("Unassembled","Unclassified"),]

total_reads=unique(unimportant_reads[,c(".id","total"),])
t=aggregate(NumReads ~ .id, unimportant_reads, sum)
unimportant_reads=merge(t,total_reads, by=".id")
unimportant_reads$Freq=unimportant_reads$NumReads / unimportant_reads$total


total_reads=unique(unimportant_reads[,c(".id","total"),])
t=aggregate(NumReads ~ .id, unimportant_reads, sum)
unimportant_reads=merge(t,total_reads, by=".id")
unimportant_reads$Freq=unimportant_reads$NumReads / unimportant_reads$total

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### barplots - Filtered - removing human, unclassfied and unlikely

classif_red=classif_fullnames[!classif_fullnames$taxon_phylum %in% 
                                c("unclassified","artificial sequences","cellular organisms", "Unassembled", "Homo_sapiens",
                                  "Eukaryota-unclassified","Eukaryota-Chordata","Eukaryota-Arthropoda","Eukaryota-Streptophyta",
                                  "Eukaryota-Mollusca", "Eukaryota-Bryozoa") ,]

key_from_data=rename_taxa(classif_red, nbact = 5,nvir = 4,neuk = 5)
#key_from_data$final_name = ifelse(str_detect(key_from_data$taxon_full_name,"myco"), 'Eukaryota-Fungi',key_from_data$final_name)


color_dict=Taxonomy_small_color_palette(key_from_data, "final_name")
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

tax=unique(classif_red$final_name)
tax_order=c(
  "artificial sequences" , "unclassified","cellular organisms", "Unassembled", "Homo_sapiens",
  tax[str_detect(tax,"Bacteria") & ! str_detect(tax, "Other")], "Other Bacteria", "Archaea",
  #tax[str_detect(tax,"myco") & ! str_detect(tax, "Other")], "Other Eukaryota - Fungi",
  tax[str_detect(tax,"Eukaryota") & ! str_detect(tax, "Other") ], "Other Eukaryota",
  tax[str_detect(tax,"Viruses") & ! str_detect(tax, "Other")], "Other Viruses")

classif_red$final_name <- factor(classif_red$final_name, levels=tax_order)

order_patient=rownames(n_reads)[order(n_reads[,1], decreasing = FALSE)]
classif_red$.id <- factor(classif_red$.id, levels=order_patient)


classif_red_blast=classif_red[classif_red$Method == "trinityBlast" & classif_red$final_name == "Viruses-Kitrinoviricota", ]
blast_order=classif_red_blast[order(classif_red_blast$Freq) ,]$.id
classif_red$.id = factor(classif_red$.id, levels = blast_order)

classif_red=merge(label_df, classif_red, by="Method")


p1=ggplot(classif_red, aes(x=.id, y=Freq, fill = final_name) ) +
  geom_bar(position="stack", stat="identity", width = 0.95)+#color="darkgrey") +
  theme_linedraw() + 
  scale_fill_manual(values = color_dict) + 
  guides(fill=guide_legend("Superkingdom - Phylum" , ncol=1)) +
  xlab("Individual") + 
  ylab(paste("Read frequency")) + 
  coord_flip() + 
  facet_nested(cols = vars(family_method, labels2lines), 
               scales="free")



p2=ggplot(classif_red, aes(y=.id, x=NumReads, fill = final_name) ) +
  geom_bar(position="stack", stat="identity", width = 0.95)+#color="darkgrey") +
  theme_linedraw() + 
  scale_fill_manual(values = color_dict) + 
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
  guides(fill=guide_legend("Superkingdom - Phylum" , ncol=1)) +
  ylab("Individual") + 
  xlab(paste("Read Number")) + 
  scale_x_continuous(labels = scales::comma) +
  facet_nested(cols = vars(family_method, labels2lines), 
               scales="free")
#ggbreak::scale_x_break(c(1e6, 1e6), scales = 0.5)


theme_fig6(p1)
theme_fig6(p2)


counting2=aggregate(NumReads ~ Method+.id ,classif_red, sum)
counting3=aggregate(NumReads ~ Method ,counting2, mean)

ggsave(p1,filename=paste0("~/results/Douek_Cleveland/figures/distribution_proportion_Cleveland.png"),
       device = "png",width = 48, height = 27, units = "cm", create.dir = TRUE)
ggsave(p2,filename=paste0("~/results/Douek_Cleveland/figures/distribution_total_Cleveland.png"),
       device = "png",width = 48, height = 27, units = "cm")


ggsave(theme_fig6(p1),
       filename = paste0("~/results/Simulation/figures_review/fig6_filtered_prop.jpeg"),
       device = "jpeg",
       width = 30,
       height = 20,
       dpi = 600,
       units = "cm",
       create.dir = TRUE)

