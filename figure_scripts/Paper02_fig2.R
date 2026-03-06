library(ggh4x)


path="~/"
sim_set="sequencing_depth"
load(paste0(path, 'data/Simulation/environnement_',sim_set,".RData"))




#colorRampPalette(c("#FDBF6F","#E66F00"))(3)
c("#FDBF6F", "#F19737" ,"#E66F00")
colorRampPalette(c("#CAB2D6","#6A3D9A"))(3)
c("#CAB2D6", "#AF8BCC", "#9464C2" ,"#7A3DB8")

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
  "Blast" = "#4B006E",       
  "SpadesBlast" = "#5e1e02",
  
  "Hybrid - Blast + Kraken2" = "#4B006E", 
  "kraken2" = "#730F0F",    
  "Hybrid - Blast + KrakenUniq" = "black",  
  "kuniq" = "#400808",  
  "Hybrid - KrakenUniq + Kraken2" = "#112904",
  "Hybrid - Kraken2 + KrakenUniq" = "#112904",  
  
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
  
  "Hybrid - KrakenUniq + Kraken2" = "- KUniq-K2 -\n(microbialDB + nt)",
  "Hybrid - Kraken2 + KrakenUniq" = "- K2-KUniq -\n(nt + microbialDB)",

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

unclassified_global <- aggregate(cbind(Freq) ~ method + .id + classification + rt + nrt + nrt_label + labels2lines + Human_as_proteobact + family_method,
                                 unclassified_sk, sum)
total_df=as.data.frame(list(
  ".id"=c("A","C","E","M","N","O","P","Q","R"),
  "total"=c(380000,760000,7600000,333000,3330000 ,1665000,16650000,185000,1850000 )
))
unclassified_global=merge(unclassified_global,total_df, by=".id")
unclassified_global$prop=unclassified_global$Freq/unclassified_global$total


# unclassified_global=merge(label_df, unclassified_global, by="method")

# 
# unclassified_global$family_method=ifelse(unclassified_global$method %in% c("Blast","SpadesBlast","kraken2", "kuniq", "Centrifuge","Centrifuger") , 
#                                          "Single Classifier" , "Hybrid")
# 
# unclassified_global$family_method=ifelse(unclassified_global$method %in% c("Hybrid - Kraken2 + KrakenUniq","Hybrid - KrakenUniq + Kraken2") , 
#                                          "Combined assembly-free" , unclassified_global$family_method)
# 
# metrics_full_red=merge(metrics_full_red, label_df, by="method")
# 
# unclassified_global$family_method=ifelse(unclassified_global$method %in% c("Hybrid - Kraken2 + KrakenUniq","Hybrid - KrakenUniq + Kraken2") , 
#                                          "Combined assembly-free" , unclassified_global$family_method)
# 
# metrics_full_red$family_method=factor(metrics_full_red$family_method, levels = c("Single Classifier","Hybrid","Combined assembly-free") )
# unclassified_global$family_method=factor(unclassified_global$family_method, levels = c("Single Classifier","Hybrid","Combined assembly-free") )

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


####___________________________________ Figures pour rep A C E - Figure 2 papier
upper_firstplot=ggplot(
  metrics_full_red[metrics_full_red$taxlvl == "phylum" & 
                     metrics_full_red$family_method != "Combined assembly-free"
                   #& metrics_full_red$method %in% c("Blast","kraken2","kuniq") 
                   , ], 
  aes(x = nrt, y = F1, fill = method, color=method)) +  
  geom_boxplot(outlier.size = .5, width=0.9)+
  theme_linedraw() +
  xlab("") + ylab("F1-score") +   
  coord_cartesian(ylim=c(.25,1))+
  guides(fill="none", color="none")+
  scale_fill_manual(values=palette, labels=labels, name="")+  
  scale_colour_manual(values=color, labels=labels, name="")+
  theme(strip.placement = "outside") +
  facet_nested_wrap(~ family_method + labels2lines , nrow = 1 )  
upper_firstplot=theme_fig2(upper_firstplot)
upper_firstplot



under_firstplotCD=ggplot(unclassified_global[unclassified_global$family_method == "Combined assembly-free",],
                         aes(x=nrt, y=prop, fill=classification))+
  geom_bar(stat = "identity") + 
  xlab("Reads Simulated per transcripts")+ ylab("") +
  scale_fill_manual(values=c(
    "U"="#FFB462", "FP"= "#B12424", "TP"="#7AADD0"),
    labels = c(
      "Unclassified", "False prediction", "True Prediction"),
    name="" )+
  theme_linedraw() +   guides(fill="none", color="none")+
  facet_nested_wrap(vars(family_method, labels2lines), nrow=1)
under_firstplotCD=theme_fig2(under_firstplotCD)
under_firstplotCD

upper_firstplotCD=ggplot(
  metrics_full_red[metrics_full_red$taxlvl == "phylum" & 
                     metrics_full_red$family_method == "Combined assembly-free"
                   #& metrics_full_red$method %in% c("Blast","kraken2","kuniq") 
                   , ], 
  aes(x = nrt, y = F1, fill = method, color=method)) +  
  geom_boxplot(outlier.size = .5, width=0.9)+
  theme_linedraw() +
  xlab("") + ylab("") +   
  coord_cartesian(ylim=c(.25,1))+
  guides(fill="none", color="none")+
  scale_fill_manual(values=palette, labels=labels, name="")+  
  scale_colour_manual(values=color, labels=labels, name="")+
  theme(strip.placement = "outside") +
  facet_nested_wrap(~ family_method + labels2lines , nrow = 1 )  
upper_firstplotCD=theme_fig2(upper_firstplotCD)
upper_firstplotCD

medians=metrics_full_red[metrics_full_red$method == "Blast",]
medians=aggregate(F1 ~ method + .id, metrics_full_red, median)

under_firstplot=ggplot(unclassified_global[unclassified_global$family_method != "Combined assembly-free" ,],
                       aes(x=nrt, y=prop, fill=classification))+
  geom_bar(stat = "identity") + 
  xlab("Reads Simulated per transcripts")+
  ylab("Proportion of sequence") +
  scale_fill_manual(values=c(
    "U"="#FFB462", "FP"= "#B12424", "TP"="#7AADD0"),
    labels = c(
      "Unclassified      ", "False prediction        ", "True Prediction"),
    guide = guide_legend(keywidth = 2, keyheight = 2, size=28, nrow = 1, ),
    name="" )+
  theme_linedraw() +
  facet_nested_wrap(vars(family_method, labels2lines), nrow=1)
under_firstplot=theme_fig2(under_firstplot) + 
  theme(legend.text = element_text(size=12, face="bold"))

under_firstplot



legend_components <- cowplot::get_plot_component(under_firstplot,"guide-box", return_all = TRUE)


final_plot=cowplot::plot_grid(rel_heights = c(12,1),
                              plot_grid(
                                upper_firstplot + theme(legend.position = "none", axis.text.x=element_blank()),
                                upper_firstplotCD + theme(legend.position = "none",axis.text.x=element_blank(),axis.text.y=element_blank()),
                                under_firstplot + theme(legend.position = "none"),
                                under_firstplotCD + theme(legend.position = "none",axis.text.y=element_blank()),
                                ncol = 2,
                                rel_widths = c(8/10, 2/10),
                                labels = c("A", "C", "B", "D")
                              ),legend_components[[3]], ncol = 1)

final_plot




upper_firstplot=ggplot(
  metrics_full_red[metrics_full_red$taxlvl == "phylum" 
                   #& metrics_full_red$method %in% c("Blast","kraken2","kuniq") 
                   , ], 
  aes(x = nrt, y = F1, fill = method, color=method)) +  
  geom_boxplot(outlier.size = .5, width=0.9)+
  theme_linedraw() +
  xlab("") + ylab("F1-score") +   
  coord_cartesian(ylim=c(.25,1))+
  guides(fill="none", color="none")+
  scale_fill_manual(values=palette, labels=labels, name="")+  
  scale_colour_manual(values=color, labels=labels, name="")+
  theme(strip.placement = "outside") +
  facet_nested_wrap(~ family_method + labels2lines , nrow = 1 )  
upper_firstplot=theme_fig2(upper_firstplot)
upper_firstplot


under_firstplot=ggplot(unclassified_global,
                       aes(x=nrt, y=prop, fill=classification))+
  geom_bar(stat = "identity") + 
  xlab("Reads Simulated per transcripts")+
  ylab("Proportion of sequence") +
  scale_fill_manual(values=c(
    "U"="#FFB462", "FP"= "#B12424", "TP"="#7AADD0"),
    labels = c(
      "Unclassified      ", "False prediction        ", "True Prediction"),
    guide = guide_legend(keywidth = 2, keyheight = 2, size=28, nrow = 1, ),
    name="" )+
  theme_linedraw() +
  facet_nested_wrap(vars(family_method, labels2lines), nrow=1)
under_firstplot=theme_fig2(under_firstplot) + 
  theme(legend.text = element_text(size=12, face="bold"))

under_firstplot



final_plot=cowplot::plot_grid(
  upper_firstplot + theme(legend.position = "none", axis.text.x=element_blank(),
                          strip.text = element_text(size=6.25)),
  under_firstplot + theme(strip.text = element_text(size=6.25)),
  ncol = 1,
  labels = c("A", "B") )



final_plot
fig2_jpgname="~/results/Simulation/figures_review/fig2.jpeg"
ggsave(final_plot,
       filename=fig2_jpgname,
       device = "jpeg",width = 30, height = 20,dpi = 350 , units = "cm", create.dir = TRUE)

magick::image_read(fig2_jpgname)



fig2_svgname="~/results/Simulation/figures_review/fig2.svg"
ggsave(final_plot,
       filename = fig2_svgname,
       device = "svg",
       width = 30,
       height = 20,
       dpi = 600,
       units = "cm",
       create.dir = TRUE)


