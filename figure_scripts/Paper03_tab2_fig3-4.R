
path="~/"
sim_set="human_majority"
load(paste0(path, 'data/Simulation/environnement_',sim_set,".RData"))


palette <- c(
  "Blast" = "#1F78B4",      # Blue
  "SpadesBlast" = "#A6CEE3",
  
  "Hybrid - Blast + Kraken2" = "#CAB2D6",  # Lighter Green
  "kraken2" = "#FB9A99",    # Bright Red (unchanged)
  "Hybrid - Blast + KrakenUniq" = "#6A3D9A",  # Darker, richer green
  "kuniq" = "#E31A1C", # Deep, bold red
  "Hybrid - KrakenUniq + Kraken2" = "#33A02C",
  "Hybrid - Kraken2 + KrakenUniq" = "#B2DF8A",
  # "Centrifuge" = "#FDBF6F" ,
  # "Centrifuger" = "#FF7F00" , 
  
  "hybrid_SpadesBlast-kraken2" = "#FDBF6F" ,
  "hybrid_SpadesBlast-kuniq" = "#FF7F00" #, 
)

color <- c(
  "Blast" = "#084269",       # Darker blue to complement Blast
  "SpadesBlast"="#053158",
  
  "Hybrid - Blast + Kraken2" = "#6A0080",  # Darker green to match hybrid
  "kraken2" = "#940C00",     # Darker but distinguishable red
  "Hybrid - Blast + KrakenUniq" = "#4B006E",  # Darkest green for this hybrid
  "kuniq" = "#730F0F",  # Rich deep red for KrakenUniq
  "Hybrid - KrakenUniq + Kraken2" = "#1F4A07", # Deep purple for hybrid
  "Hybrid - Kraken2 + KrakenUniq" = "#275c09",  # Even deeper purple for hybrid
  # "Centrifuge" = "orangered4" ,
  # "Centrifuger" = "saddlebrown" , 
  
  "hybrid_SpadesBlast-kraken2" = "orangered4" ,
  "hybrid_SpadesBlast-kuniq" = "saddlebrown" #, 
)



palette <- c(
  "Blast" = "#E3D7EB",                        
  "SpadesBlast" = "#FDD9A0",
  
  "Hybrid - Blast + Kraken2" = "#AF8BCC",     # medium purple
  "kraken2" = "#FB9A99",    # Bright Red (unchanged)
  "Hybrid - Blast + KrakenUniq" = "#6A3D9A",  # dark purple (kuniq)
  "kuniq" = "#E31A1C", # Deep, bold red
  "Hybrid - KrakenUniq + Kraken2" = "#33A02C",
  "Hybrid - Kraken2 + KrakenUniq" = "#B2DF8A",
  
  "hybrid_SpadesBlast-kraken2" = "#F19737",    # medium orange
  "hybrid_SpadesBlast-kuniq" = "#E66F00"       # dark orange (kuniq)
)

color <- c(
  "Blast" = "#4B006E",       # Darker blue to complement Blast
  "SpadesBlast" = "#5e1e02",
  
  "Hybrid - Blast + Kraken2" = "#4B006E",  # Darker green to match hybrid
  "kraken2" = "#940C00",     # Darker but distinguishable red
  "Hybrid - Blast + KrakenUniq" = "black",  # Darkest green for this hybrid
  "kuniq" = "#730F0F",  # Rich deep red for KrakenUniq
  "Hybrid - KrakenUniq + Kraken2" = "#1F4A07", # Deep purple for hybrid
  "Hybrid - Kraken2 + KrakenUniq" = "#1F4A07",  # Even deeper purple for hybrid
  # "Centrifuge" = "orangered4" ,
  # "Centrifuger" = "saddlebrown" , 
  
  "hybrid_SpadesBlast-kraken2" = "#5e1e02" ,
  "hybrid_SpadesBlast-kuniq" = "black" #, 
)


labels <- c(
  "Blast" = "Contig-based (nt database)",
  "SpadesBlast" = "SpadesBlast",
  
  "Hybrid - Blast + Kraken2" = "Hybrid - Blast + Kraken2",
  "kraken2" = "Kmer-based (Kraken2 - nt database)",     
  
  "Hybrid - Blast + KrakenUniq" = "Hybrid - Blast + KrakenUniq",
  "kuniq" = "Kmer-based (Krakenuniq - microbialDB database)",
  "Hybrid - KrakenUniq + Kraken2" = "Kmer-based (KrakenUniq + Kraken2)",
  "Hybrid - Kraken2 + KrakenUniq" = "Kmer-based (Kraken2 + KrakenUniq)",
  # "Centrifuge" = "- Centrifuge - \n(nt)" ,
  # "Centrifuger" = "- Centrifuger - \n(RefSeq Bacteria-Viruses)" ,
  
  # Hybrid spades
  "hybrid_SpadesBlast-kraken2" =   "- Hybrid-Spades-K2 -\n(nt)",
  "hybrid_SpadesBlast-kuniq"=  "- Hybrid-Spades-KUniq -\n(nt + microbialDB)"
  
)

labels2lines <- c(
  "Blast" = "- Trinity-Blast -\nnt ",
  "SpadesBlast" = "- SPAdes-Blast -\nnt ",
  
  "Hybrid - Blast + Kraken2" = "- Hybrid-Trinity-K2 -\nnt",
  "kraken2" = "- K2 -\nnt", 
  
  "Hybrid - Blast + KrakenUniq" = "- Hybrid-Trinity-KUniq -\nnt + microbialDB",
  "kuniq" = "- KUniq -\nmicrobialDB",
  
  "Hybrid - KrakenUniq + Kraken2" = "- KUniq-K2 -\nmicrobialDB + nt",
  "Hybrid - Kraken2 + KrakenUniq" = "- K2-KUniq -\nnt + microbialDB",
  
  # "Centrifuge" = "- Centrifuge - \n(nt)" ,
  # "Centrifuger" = "- Centrifuger - \n(RefSeq Bacteria-Viruses)",
  
  # Hybrid spades
  "hybrid_SpadesBlast-kraken2" =    "- Hybrid-SPAdes-K2 -\nnt",
  "hybrid_SpadesBlast-kuniq"=  "- Hybrid-SPAdes-KUniq -\nnt + microbialDB"
  # "Fulgor" = "- Fulgor -\n (Refseq Viruses ?)"
)

label_df=as.data.frame(labels2lines)
label_df$method=rownames(label_df)

label_df$family_method = c("Assembly-based","Assembly-based", 
                           "Hybrid" , "Assembly-free",
                           "Hybrid" , "Assembly-free",
                           "Combined assembly-free", "Combined assembly-free",
                           #"Assembly-free", "Assembly-free",
                           "Hybrid","Hybrid")


label_df$family_method = factor(label_df$family_method, levels = c("Assembly-based", "Assembly-free" ,"Hybrid", "Combined assembly-free") )

label_df$method = factor(label_df$method, levels = c("Blast" , "SpadesBlast" , 
                                                     "kraken2" , "kuniq", 
                                                     
                                                     "Hybrid - Blast + Kraken2" , 
                                                     "hybrid_SpadesBlast-kraken2" , 
                                                     
                                                     "Hybrid - Blast + KrakenUniq" , 
                                                     "hybrid_SpadesBlast-kuniq" , 
                                                     
                                                     "Hybrid - Kraken2 + KrakenUniq",
                                                     "Hybrid - KrakenUniq + Kraken2") )



label_df$labels2lines <- factor(
  label_df$labels2lines,
  levels = unique(label_df$labels2lines[order(label_df$method)])
)


levels(label_df$method)
levels(label_df$labels2lines)


#_______________________________________________ Counting unclassified / TP / FP
### Unclassified per superkingdom ###
# Fetch uncassified/unassembled per superkingdom

classifs_pair = classifs_pair[! str_detect(classifs_pair$method, "Centrifuge"),]
classifs_pair$superkingdom_truth=sapply( classifs_pair$taxon_full_name_truth, function(x){str_split(x, " - ")[[1]][1]} )
classifs_pair$superkingdom_truth=ifelse( str_detect(classifs_pair$taxon_full_name_truth, "Eukaryota - Chordata"),
                                         "Homo sapiens" , classifs_pair$superkingdom_truth)

classifs_pair$superkingdom_pred=sapply( classifs_pair$taxon_full_name, function(x){str_split(x, " - ")[[1]][1]} )
classifs_pair$superkingdom_pred=ifelse( str_detect(classifs_pair$taxon_full_name, "Eukaryota - Chordata"),
                                        "Homo sapiens" , classifs_pair$superkingdom_pred)

classifs_pair$taxon_full_name_truth=ifelse( str_detect(classifs_pair$taxon_full_name_truth, "Eukaryota - Chordata"),
                                            "Homo sapiens" , classifs_pair$taxon_full_name_truth)
classifs_pair$taxon_full_name=ifelse( str_detect(classifs_pair$taxon_full_name, "Eukaryota - Chordata"),
                                      "Homo sapiens" , classifs_pair$taxon_full_name)

#unclassified_sk=classifs_pair[classifs_pair$taxon_full_name %in% c("Unclassified","Unassembled") & classifs_pair$prop > 0,]
classifs_pair$prefiltering=ifelse(str_detect(classifs_pair$method,"human_kuniq") ,"KrakenUniq Human Filter", "No filter" )
classifs_pair$method=str_remove(classifs_pair$method, "human_kuniq_")

# Removing filter / no prefilter (not kept for following results)
classifs_pair=classifs_pair[classifs_pair$prefiltering == "No filter",]
unclassified_sk=classifs_pair[classifs_pair$prop > 0,]


# Creating 3 classes: TP / FP / missing
unclassified_sk$classification=ifelse(unclassified_sk$taxon_full_name_truth == unclassified_sk$taxon_full_name , "TP", "FP")
unclassified_sk$classification=ifelse(unclassified_sk$taxon_full_name %in% c("Unclassified","Unassembled") , 
                                      "U", unclassified_sk$classification)


# Aggregating classif at phylum
tmp=unclassified_sk
tmp$label_truth=ifelse(tmp$classification == "FP" & tmp$superkingdom_truth == "Bacteria" , 
                       "Misclassified Bacteria", tmp$taxon_full_name_truth)
tmp$label_truth=ifelse(tmp$classification == "FP" & tmp$superkingdom_truth == "Eukaryota" , 
                       "Misclassified Eukaryota", tmp$label_truth)
tmp$label_truth=ifelse(tmp$classification == "FP" & tmp$superkingdom_truth == "Archaea" , 
                       "Misclassified Archaea", tmp$label_truth)
tmp$label_truth=ifelse(tmp$classification == "FP" & tmp$superkingdom_truth == "Viruses" , 
                       "Misclassified Viruses", tmp$label_truth)
tmp$label_truth=ifelse(tmp$classification == "FP" & tmp$taxon_full_name_truth == "Viruses - Hofneiviricota" , 
                       "Viruses - Hofneiviricota", tmp$label_truth)
unique(tmp[tmp$taxon_full_name == "Bacteria - Pseudomonadota",]$label_truth)


freq=aggregate(Freq ~ .id + method + taxon_full_name + classification + prefiltering + label_truth , 
               tmp , sum)

tmp=unique(unclassified_sk[,c(1,2,3,10,11,13)])
total=aggregate(total ~ .id + taxon_full_name_truth + method + prefiltering , tmp, sum)
unclassified_phylum=freq


# Aggregating classif at Superkingdom level
freq=aggregate(Freq ~ .id + method + superkingdom_truth + classification + prefiltering, unclassified_sk, sum)
tmp=unique(unclassified_sk[,c(1,2,3,10,11,13)])
total=aggregate(total ~ .id + superkingdom_truth + method + prefiltering , tmp, sum)

unclassified_sk=merge(freq,total, by=c("superkingdom_truth","method","prefiltering",".id"))
unclassified_sk$prop=unclassified_sk$Freq / unclassified_sk$total


###___________________________________________________________ Metrics dataframe

metrics_full=do.call(rbind, metrics_list)
metrics_full = metrics_full[! str_detect(metrics_full$method, "Centrifuge"),]

metrics_full$.id <- substr(metrics_full$replicate_id, 1, 1) # Extract the first character (letter)

colnames(metrics_full)[colnames(metrics_full) == "taxon"]="taxon_full_name"
metrics_full$phylum <- sapply(strsplit(metrics_full$taxon_full_name, "_"), function(x) x[2])
metrics_full$superkingdom <- sapply(strsplit(metrics_full$taxon_full_name, "_"), function(x) x[1])

metrics_full=unique(merge(metrics_full[metrics_full$taxlvl == "phylum", colnames(metrics_full) != "taxon_full_name"],
                          indexes, by=c("phylum","replicate_id",'taxlvl')))


## Creating "alpha_global" : mean number of species for each replicate, per condition
metrics_full = metrics_full %>%
  group_by(taxon_full_name, .id, method) %>%
  mutate(alpha_global=mean(alpha_div))

metrics_full$phylum = ifelse(metrics_full$phylum == "Chordata", 
                             "Homo sapiens", metrics_full$phylum)

metrics_full$superkingdom = ifelse(metrics_full$phylum == "Chordata", 
                                   "", metrics_full$superkingdom)

metrics_full$taxon_full_name = ifelse(metrics_full$phylum == "Chordata", 
                                      "Homo sapiens", metrics_full$taxon_full_name)


metrics_full$phylum = ifelse(metrics_full$phylum == "unclassified Viruses phylum", "Caudoviricetes", metrics_full$phylum)
metrics_full$label=paste0(metrics_full$superkingdom,"\n",metrics_full$phylum, "\n mean n_family = ", round(metrics_full$alpha_global,1)  )

metrics_full$label_factor <- factor(metrics_full$label, levels = unique(metrics_full$label[order(metrics_full$alpha_global)]) )
metrics_full$taxlvl <- factor(metrics_full$taxlvl, levels = c("superkingdom", "phylum", "class", "order", "family") )

################################################################################
### Plots
metrics_full$nrt=as.numeric(str_remove_all(metrics_full$rt, "r/t"))
metrics_full$nrt=factor(as.character(metrics_full$nrt),levels = c("5","10","100"))


metrics_full=merge(metrics_full, misclassified_human, by=c(".id","method"))


### Renaming method names
metrics_full$prefiltering=ifelse(str_detect(metrics_full$method,"human_kuniq") , 
                                 "KrakenUniq Human Filter",
                                 "No filter" )
metrics_full$method=str_remove(metrics_full$method, pattern = "human_kuniq_")
metrics_full = metrics_full[metrics_full$prefiltering == "No filter" , ]

metrics_full_kuniq=metrics_full[metrics_full$method == "kuniq",]
metrics_full_kuniq$prefiltering="KrakenUniq Human Filter"
metrics_full=rbind(metrics_full, metrics_full_kuniq)


metrics_full_red = metrics_full
metrics_full_red$nrt_label=paste0(metrics_full_red$nrt, "\n Reads/transcripts")

metrics_full_red$method = ifelse(metrics_full_red$method == "hybrid_Blast-kuniq", "Hybrid - Blast + KrakenUniq", metrics_full_red$method)
metrics_full_red$method = ifelse(metrics_full_red$method == "hybrid_Blast-kraken2", "Hybrid - Blast + Kraken2", metrics_full_red$method)

metrics_full_red$method = ifelse(metrics_full_red$method == "hybrid_kuniq-kraken2", "Hybrid - KrakenUniq + Kraken2", metrics_full_red$method)
metrics_full_red$method = ifelse(metrics_full_red$method == "hybrid_kraken2-kuniq", "Hybrid - Kraken2 + KrakenUniq", metrics_full_red$method)

print(unique(metrics_full_red$method))


interesting_col=unique(metrics_full_red[,c(1,2,12,21,22,23,24)])


unclassified_phylum$method = ifelse(unclassified_phylum$method == "hybrid_Blast-kuniq", "Hybrid - Blast + KrakenUniq", unclassified_phylum$method)
unclassified_phylum$method = ifelse(unclassified_phylum$method == "hybrid_Blast-kraken2", "Hybrid - Blast + Kraken2", unclassified_phylum$method)
unclassified_phylum$method = ifelse(unclassified_phylum$method == "hybrid_kuniq-kraken2", "Hybrid - KrakenUniq + Kraken2", unclassified_phylum$method)
unclassified_phylum$method = ifelse(unclassified_phylum$method == "hybrid_kraken2-kuniq", "Hybrid - Kraken2 + KrakenUniq", unclassified_phylum$method)

unclassified_sk$method = ifelse(unclassified_sk$method == "hybrid_Blast-kuniq", "Hybrid - Blast + KrakenUniq", unclassified_sk$method)
unclassified_sk$method = ifelse(unclassified_sk$method == "hybrid_Blast-kraken2", "Hybrid - Blast + Kraken2", unclassified_sk$method)
unclassified_sk$method = ifelse(unclassified_sk$method == "hybrid_kuniq-kraken2", "Hybrid - KrakenUniq + Kraken2", unclassified_sk$method)
unclassified_sk$method = ifelse(unclassified_sk$method == "hybrid_kraken2-kuniq", "Hybrid - Kraken2 + KrakenUniq", unclassified_sk$method)


unclassified_sk = merge(unclassified_sk, interesting_col, by=c(".id","method","prefiltering"))
unclassified_sk$classification=factor(unclassified_sk$classification, levels = rev(c("TP","FP","U")))
unclassified_sk=merge(unclassified_sk, label_df, by = "method")


### Proportion of unclassified human over every sequence
total_seq=unique(aggregate(Freq ~ .id + method + prefiltering, unclassified_sk, sum)[,c(1,4)])
colnames(total_seq)[2]="total_seq"

unclassified_human=merge(unclassified_sk,total_seq, by=".id")
unclassified_human$freq_total=unclassified_human$Freq / unclassified_human$total_seq
unclassified_human=unclassified_human[unclassified_human$classification != "TP" & unclassified_human$superkingdom_truth == "Homo sapiens",]


####
unclassified_phylum = merge(unclassified_phylum, interesting_col, by=c(".id","method","prefiltering"))
unclassified_phylum$classification=factor(unclassified_phylum$classification, levels = rev(c("TP","FP","U")))
unclassified_phylum=merge(unclassified_phylum, label_df, by = "method")
unclassified_phylum$superkingdom_pred=sapply( unclassified_phylum$taxon_full_name, function(x){str_split(x, " - ")[[1]][1]} )

unclassified_phylum$prefiltering = factor(unclassified_phylum$prefiltering, levels=rev(c("KrakenUniq Human Filter" , "No filter") ))
####



metrics_full_red$method = factor(metrics_full_red$method, levels = names(color))

# Ensure 'superkingdom' is a factor (if not already) with the desired order
metrics_full_red$superkingdom <- factor(metrics_full_red$superkingdom, 
                                        levels = c("Bacteria", "Archaea", "Eukaryota", "Viruses"))


metrics_full_red <- metrics_full_red %>%
  arrange(superkingdom, alpha_div) %>% # Use dplyr for ordering
  mutate(label_factor = factor(label_factor, levels = unique(label_factor))) # Adjust levels

#levels(metrics_full_red$label_factor)

#______________________________________________________________________________#
###### Performances depending on reads simulated per transcripts
### 100rt - some phylum x: method // y= F1

levels(metrics_full_red$method)

metrics_full_red = merge(label_df, metrics_full_red, by="method")

metrics_full_red$strip_label_skmethod=paste0(metrics_full_red$superkingdom, "\n", metrics_full_red$family_method)
metrics_full_red$method_nofactor=as.character(metrics_full_red$method)

metrics_full_red$interaction_var <- interaction(metrics_full_red$family_method, metrics_full_red$method_nofactor)
metrics_full_red$interaction_method <- interaction(metrics_full_red$family_method, metrics_full_red$method_nofactor)


plot=ggplot(metrics_full_red[
  metrics_full_red$superkingdom != "Archaea" &
    metrics_full_red$taxlvl=="phylum" &
    metrics_full_red$prefiltering == "No filter" &
    metrics_full_red$nrt == "100",] ,
  aes(x=labels2lines, y=F1, fill=method_nofactor, colour=method_nofactor))+#, position=phylum  ) +

  geom_boxplot(outlier.size = .5 , width=.7, position = position_dodge(width = 0.75))+ #geom_jitter(color="black", size=0.4, alpha=0.9)+
  theme_linedraw()+
  xlab("Family of approaches")+ ylab("F1-score") +
  coord_cartesian(ylim=c(.5,1))+
  scale_fill_manual(values=palette, labels=labels, name="")+
  scale_colour_manual(values=color, labels=labels, name="")+
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        strip.background = element_rect(fill="#333333") ,
        legend.key.size = unit(1, "cm"),
        legend.key.width = unit(1,"cm"),
        legend.text = element_text(size = 11, face = "bold"),
        strip.text = element_text(size = 10, face = "bold"),
        axis.title.y= element_text(size = 12, face = "bold"),
        axis.title.x= element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),

  )+
  #facet_nested_wrap(~ superkingdom + family_method, scales = "free_x", nrow = 1) +
  facet_grid( vars(nrt_label), vars(superkingdom), scales = "free_x")

plot



########________________________________________________________________________

unclassified_global <- aggregate(cbind(Freq) ~ method + .id + classification + rt + nrt + nrt_label + labels2lines + Human_as_proteobact,
                                 unclassified_sk, sum)
total_df=as.data.frame(list(
  ".id"=c("A","C","E","M","N","O","P","Q","R"),
  "total"=c(380000,760000,7600000,333000,3330000 ,1665000,16650000,185000,1850000 )
))
unclassified_global=merge(unclassified_global,total_df, by=".id")

unclassified_global$prop=unclassified_global$Freq/unclassified_global$total


theme_fig2=function(plot){
  return(plot + theme(legend.direction = "horizontal",legend.position = "bottom",
                      legend.key.size = unit(1, "cm"),
                      legend.key.width = unit(1,"cm"),
                      strip.background = element_rect(fill="#333333") , 
                      strip.text = element_text(size = 7, face = "bold"),
                      axis.title.y= element_text(size = 10, face = "bold"),
                      axis.title.x= element_text(size = 10, face = "bold"),
                      axis.text.x = element_text(size = 10),
                      axis.text.y = element_text(size = 10 )))
}

#______________________________________________________________________________# Table 2 - Mean F1 +/- sd

# Aggregating results per median F1
metrics = metrics_full_red[metrics_full_red$.id == "P" & metrics_full_red$prefiltering == "No filter", ]

##### Mean F1 + standard error scores per phylum
metrics_hs = metrics
metrics_hs$superkingdom = ifelse(metrics_hs$phylum == "Homo sapiens", "Homo sapiens", as.character(metrics_hs$superkingdom))

df_F1_mean=aggregate(F1 ~ method + superkingdom + prefiltering, metrics_hs, mean)
colnames(df_F1_mean)[4]="mean_F1"
df_F1_sd=aggregate(F1 ~ method + superkingdom + prefiltering, metrics_hs, sd)
colnames(df_F1_sd)[4]="sd_F1"
df_F1=merge(df_F1_mean,df_F1_sd)



df_F1_mean_global=aggregate(F1 ~ method + prefiltering, metrics, mean)
colnames(df_F1_mean_global)[3]="mean_F1"
df_F1_sd_global=aggregate(F1 ~ method + prefiltering , metrics, sd)
colnames(df_F1_sd_global)[3]="sd_F1"
df_F1_global=merge(df_F1_mean_global,df_F1_sd_global)
df_F1_global$superkingdom="Global"

df_F1=rbind(df_F1, df_F1_global)
df_F1$mean_F1 = round(df_F1$mean_F1,3)
df_F1$sd_F1 = round(df_F1$sd_F1,2)
df_F1$label = paste0(df_F1$mean_F1, " ( +/- ",df_F1$sd_F1,")" )


df_F1=merge(label_df, df_F1)

df_F1$superkingdom = factor(df_F1$superkingdom, levels = c("Bacteria", "Eukaryota", "Viruses", "Archaea", "Homo sapiens", "Global"))


#______________________________________________________________________________# Table 2 - min max F1
##### Range of F1 scores

# metrics_small=metrics[metrics$prefiltering == "No filter" & metrics$method != "Centrifuge",]
# metrics_small$superkingdom = ifelse(metrics_small$phylum == "Homo sapiens", "Homo sapiens", as.character(metrics_small$superkingdom))

metrics = metrics[metrics$prefiltering == "No filter", ]
metrics_hs = metrics
metrics_hs$superkingdom = ifelse(metrics_hs$phylum == "Homo sapiens", "Homo sapiens", as.character(metrics_hs$superkingdom))

df_F1_global=aggregate(F1 ~ method + replicate_id + method + family_method + labels2lines, metrics_hs ,mean)         
df_F1_global$superkingdom="Global"

df_F1_mean=aggregate(F1 ~ method + replicate_id + superkingdom + method + family_method + labels2lines, metrics_hs ,mean)         
df_F1_mean = rbind(df_F1_mean, df_F1_global)



df_range <- df_F1_mean %>%
  group_by(superkingdom, method, family_method, labels2lines) %>%
  summarise(F1_min = round(min(F1),2), F1_max = round(max(F1),2), .groups = "drop")

df_range$range=paste0("[",df_range$F1_min," ; ", df_range$F1_max,"]" )
df_range$middle = (df_range$F1_min + df_range$F1_max) / 2
df_range=merge(label_df, df_range)

df_range$superkingdom = factor(df_range$superkingdom, levels = c("Bacteria", "Eukaryota", "Viruses", "Archaea", "Homo sapiens", "Global"))

heat=ggplot(df_range, 
            aes(y=superkingdom, x=labels2lines, fill=middle))+
  geom_tile() + geom_text(aes(label=range)) +
  ylab("")+
  scale_fill_gradientn(
    colors = c("red", "white", "green3"),
    values = scales::rescale(c(min(df_range$middle), mean(df_range$middle), max(df_range$middle)))
  ) + labs(fill='Mean Precision') +
  theme_classic()+
  facet_nested(vars(superkingdom), vars(family_method, labels2lines) ,scales="free")

heat = theme_fig2(heat) + 
  theme(legend.position="right", 
        legend.direction = "vertical",
        strip.text = element_text(color="white"),
        axis.text.x = element_blank())

heat

# ggsave(heat,filename=paste0("~/results/Simulation/with_replacement/figures/heatmap_precision_filtering.jpeg"),
#        device = "jpeg",width = 40, height = 32,dpi = 350 , units = "cm", create.dir = TRUE)



pivoted_range <- reshape(df_range[,c("method","superkingdom","range")],
                         idvar = "superkingdom",
                         timevar = "method",
                         direction = "wide")
pivoted_range$superkingdom = factor(pivoted_range$superkingdom, levels=c("Bacteria","Eukaryota","Viruses","Archaea","Homo sapiens","Global"))



write.table(df_range[,c(1,3,4,5,6,7,8)], file = "~/results/Simulation/figures_review/table2_data_reviews.tsv",
           quote=FALSE, row.names=FALSE,sep="\t")


#______________________________________________________________________________# Figure  4


#########################""_____________________________________________________
### O V E R   R E P R E S E N T A T I O N   O F   H U M A N   R E A D S
### 10rt - some phylum x: % of human reads // y= F1 // facet_grid : phyla + method
human_unc=unclassified_sk[unclassified_sk$superkingdom_truth == "Homo sapiens",]

global_compo=as.data.frame(list(
  "percentage_human"=c(10,10,50,50,90,90),
  ".id"=c('Q',"R", "M", "N", "O", "P")) )


unclassified_90=merge(unclassified_sk[unclassified_sk$prefiltering=="No filter",],total_df, by=".id")

t=unclassified_90$total.y  
unclassified_90$total.y=NULL
unclassified_90$total.x=NULL
unclassified_90$total=t
unclassified_90$freq_total=unclassified_90$Freq / unclassified_90$total

unclassified_90=merge(unclassified_90,global_compo, by=".id")

unclassified_90$human_label=ifelse(unclassified_90$superkingdom_truth == "Homo sapiens",
                                   "Human","non Human")

unclassified_90_reduced=unclassified_90[ unclassified_90$nrt == 100 &
                                           unclassified_90$prefiltering == "No filter" &
                                           unclassified_90$.id %in% c("N","P","R") ,] 


unclassified_90_reduced$label=interaction(unclassified_90_reduced$human_label,unclassified_90_reduced$classification)

unclassified_90_reduced$label=factor(unclassified_90_reduced$label, 
                                     levels=c(
                                       "non Human.U","Human.U","non Human.FP","Human.FP","non Human.TP","Human.TP"))


unclassified_90_reduced$sk_category=ifelse(unclassified_90_reduced$superkingdom_truth == "Homo sapiens","Homo sapiens", "non-human")

test=aggregate(freq_total ~ .id + method+sk_category+classification+labels2lines+family_method+total+percentage_human+human_label+label ,
               
               unclassified_90_reduced[unclassified_90_reduced$nrt==100,], sum)

fig4=ggplot(unclassified_90_reduced[unclassified_90_reduced$method != "Centrifuge" , ],
                    aes(x=as.character(percentage_human), y=freq_total, fill=label )) +
  geom_bar(position="stack", stat="identity", width=.75) + 
  
  geom_hline(yintercept = c(0.1,0.5,0.9), linetype="dashed", alpha=.4)+
  theme_linedraw() + 
  xlab("Percentage of human sequences") + 
  ylab("Sequence classification proportion")+
  
  scale_y_continuous(breaks = c(0.1, 0.3 ,0.5, 0.7,0.9)) +
  scale_fill_manual(
    values=c(
      "Human.U"="#D98B00",  
      "Human.FP"="#6E0000",  
      "Human.TP"="#005A8D",  
      "non Human.U"="#FFB462",
      "non Human.FP"="#B12424",
      "non Human.TP"="#7AADD0") ,
    
    labels = c(
      "Human.U"="Unclassified Human", 
      "Human.FP"="Misclassified Human", 
      "Human.TP"="Correctly predicted Human",
      "non Human.U"="Other Unclassified", 
      "non Human.FP"="Other False prediction", 
      "non Human.TP"="Other True Prediction"
    ),
    name="" )+
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(.75, "cm"),
        legend.key.width = unit(.75,"cm"),
        strip.background = element_rect(fill="#333333") , 
        legend.text = element_text(size = 8, face = "bold"),
        strip.text = element_text(size = 6.25, face = "bold"),
        axis.title.y= element_text(size = 9, face = "bold"),
        axis.text.x = element_text(size = 8, angle=0),
        axis.title.x= element_text(size = 9, face = "bold"),
        axis.text.y = element_text(size = 8 ),
        
  ) + 
  facet_nested(cols = vars(family_method, labels2lines), 
               rows=vars(human_label), scales="free")
#fig4

fig4_jpgname="~/results/Simulation/figures_review/fig4.jpeg"
ggsave(fig4,
       filename = fig4_jpgname,
       device = "jpeg",
       width = 30,
       height = 18,
       dpi = 600,
       units = "cm",
       create.dir = TRUE)
magick::image_read(fig4_jpgname)

fig4_svgname="~/results/Simulation/figures_review/fig4.svg"
ggsave(fig4,
  filename = fig4_svgname,
  device = "svg",
  width = 30,   
  height = 18,  
  units = "cm" , 
  dpi=600
)
      


################################################################################
### --------------------------------------------------------------------- Figure 3

global_compo=as.data.frame(list(
  "percentage_human"=c(10,10,50,50,90,90),
  ".id"=c('Q',"R", "M", "N", "O", "P")) )
metrics_full_red=merge(global_compo,metrics_full_red, by=".id")
metrics_full_red$replicate_number=substring(metrics_full_red$replicate_id, 2, 2)

nonente_df=metrics_full_red[metrics_full_red$percentage_human == 90 &
                              metrics_full_red$nrt == 100 & 
                              metrics_full_red$prefiltering == "No filter",]#c(2,3,5:11,13,20,25,29:32)]                           

dix_df=metrics_full_red[metrics_full_red$percentage_human == 10 & 
                          metrics_full_red$nrt == 100 & 
                          metrics_full_red$prefiltering == "No filter",]#c(2,3,5:11,20,13,25,29:32)]




colnames(dix_df)[9:14]=paste0(colnames(dix_df)[9:14],"_10")
colnames(nonente_df)[9:14]=paste0(colnames(nonente_df)[9:14],"_90")



decrease_df=merge(nonente_df, dix_df[,c(3,6,9:14, 32)], by=c("method","phylum","replicate_number") )
decrease_df$decrease_precision=decrease_df$Precision_90 - decrease_df$Precision_10
decrease_df$decrease_F1=decrease_df$F1_90 - decrease_df$F1_10
decrease_df$decrease_recall=decrease_df$Recall_90 - decrease_df$Recall_10


order_phylum=unique(metrics_full_red[#metrics_full_red$superkingdom == "Bacteria",
  ,c("phylum","label_factor","alpha_global")])
order_phylum=arrange(order_phylum, alpha_global)

decrease_df$phylum = factor(decrease_df$phylum, levels= unique(order_phylum$phylum)) 


decrease_df_summary <- aggregate(cbind(decrease_F1, decrease_precision, decrease_recall) ~ labels2lines + method + superkingdom + phylum + label_factor + family_method, 
                                 decrease_df, function(x) c("mean" = mean(x), "sd" = sd(x)))

decrease_df_summary <- do.call(data.frame, decrease_df_summary)
names(decrease_df_summary) <- gsub("\\.", "_", names(decrease_df_summary))  # Clean column names if needed





### plot manips
decrease_df_summary$lowerF1=decrease_df_summary$decrease_F1_mean-decrease_df_summary$decrease_F1_sd
decrease_df_summary$upperF1=decrease_df_summary$decrease_F1_mean+decrease_df_summary$decrease_F1_sd

decrease_df_summary$lowerPre=decrease_df_summary$decrease_precision_mean-decrease_df_summary$decrease_precision_sd
decrease_df_summary$upperPre=decrease_df_summary$decrease_precision_mean+decrease_df_summary$decrease_precision_sd

decrease_df_summary$lowerRec=decrease_df_summary$decrease_recall_mean-decrease_df_summary$decrease_recall_sd
decrease_df_summary$upperRec=decrease_df_summary$decrease_recall_mean+decrease_df_summary$decrease_recall_sd

decrease_df_summary$lowerPre=decrease_df_summary$decrease_precision_mean-decrease_df_summary$decrease_precision_sd
decrease_df_summary$upperPre=decrease_df_summary$decrease_precision_mean+decrease_df_summary$decrease_precision_sd

### Precision confidence interval:
### Standard error of the mean
t_crit=qt(1 - 0.025, df = 9)
CI_pre = decrease_df_summary$decrease_precision_sd / sqrt(9)
ME_pre = CI_pre * t_crit

decrease_df_summary$lowerCI=decrease_df_summary$decrease_precision_mean - ME_pre
decrease_df_summary$upperCI=decrease_df_summary$decrease_precision_mean + ME_pre

decrease_df_summary$full_name=paste0(decrease_df_summary$superkingdom, ' - ', decrease_df_summary$phylum)
decrease_df_summary$full_name = factor(decrease_df_summary$full_name, levels = unique(decrease_df_summary$full_name))


palette=Taxonomy_small_color_palette2(decrease_df_summary, "full_name")

palette_phylum_names=sapply(names(palette), function(x) str_split(x, pattern = " - ")[[1]][2])

decrease_df_summary$full_name=factor(decrease_df_summary$full_name , levels=names(palette))
decrease_df_summary$phylum=factor(decrease_df_summary$phylum , levels=palette_phylum_names)

human_upperpoint_pre = ggplot(
  decrease_df_summary[decrease_df_summary$superkingdom != "Archaea" & 
                        ! decrease_df_summary$method %in% c("Centrifuge") &
                        #! decrease_df_summary$method %in% c("Centrifuge","Hybrid - Blast + Kraken2","Hybrid - Kraken2 + KrakenUniq") &
                        #decrease_df_summary$family_method != "Combined assembly-free" &
                        decrease_df_summary$phylum != "Homo sapiens",],
  aes(x = superkingdom, y = decrease_precision_mean, 
      ymin=lowerCI , 
      ymax=upperCI , 
      colour = full_name)) + 
  geom_point(position = position_dodge(width = 0.75)) +
  geom_errorbar(position = position_dodge(width = 0.75), width = 0) +
  theme_linedraw()+
  ylab("Difference in precision (90% - 10%)")+ 
  xlab("") +  guides(fill="none", color="none")+
  #scale_color_brewer(palette = "Set3") +
  scale_colour_manual(values = palette, labels = names(palette), name = "") +
  facet_nested_wrap(vars(family_method,labels2lines), nrow=1, scales="free_x")

human_upperpoint_pre=theme_fig2(human_upperpoint_pre) + 
  theme(axis.text.x = element_text(size=9 , angle=25,hjust=1, vjust=1))




meandecrease_dotplot = ggplot(
  decrease_df_summary[decrease_df_summary$superkingdom == "Bacteria" & 
                        decrease_df_summary$method %in% c("Blast","kraken2"),],
  aes(x = phylum, y = decrease_precision_mean, 
      ymin=lowerCI , 
      ymax=upperCI , 
      colour = full_name)) + 
  geom_hline(yintercept = 0, linetype="dashed", 
             color = "#454545", linewidth=.75) +
  geom_point(position = position_dodge(width = 0.75),size=3) +
  geom_errorbar(position = position_dodge(width = 0.75), width = .25, linewidth=1) +
  theme_linedraw()+
  ylab("Difference in precision (90% - 10%)")+ 
  xlab("     Bacterial Phylum") +  guides(fill="none", color="none")+
  #scale_color_brewer(palette = "Set3") +
  scale_color_manual(values = palette, labels = names(palette), name = "") + coord_flip()+
  facet_nested_wrap(vars(family_method,labels2lines), nrow=1)

meandecrease_dotplot=theme_fig2(meandecrease_dotplot) + 
  theme(axis.text.y = element_text(size=10 ))
meandecrease_dotplot

fig3 <- cowplot::plot_grid(
  human_upperpoint_pre + theme(strip.text = element_text(size=6.25)),
  meandecrease_dotplot+ theme(strip.text = element_text(size=10)),
  ncol = 1,
  labels = c("A", "B") , 
  rel_heights = c(.55,.45)
)

fig3
fig3_jpgname="~/results/Simulation/figures_review/fig3.jpeg"
ggsave(fig3,
       filename = fig3_jpgname,
       device = "jpeg",
       width = 30,
       height = 20,
       dpi = 600,
       units = "cm",
       create.dir = TRUE)
magick::image_read(fig3_jpgname)


fig3_svgname="~/results/Simulation/figures_review/fig3.svg"
ggsave(fig3,
       filename = fig3_svgname,
       device = "svg",
       width = 30,   
       height = 20,  
       units = "cm" , 
       dpi=600
)
