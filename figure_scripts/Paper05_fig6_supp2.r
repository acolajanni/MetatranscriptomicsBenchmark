library(stringr)
library(plyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(tidyr)
library(patchwork) 
library(cowplot)
library(ggh4x)
library(tidyr)
library(dplyr)
library(tibble)


################################################################################
# F U N C T I O N S
################################################################################
ranks=c("superkingdom","phylum", "class", "order", "family", "genus", "species")
ncores=19

#metrics_full_red$method = ifelse(metrics_full_red$method == "kuniq", "KrakenUniq" , metrics_full_red$method)

# RColorBrewer::brewer.pal(n = 14, name = "Paired")
# [1] "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C"
# [7] "#FDBF6F" "#FF7F00" "#CAB2D6" "#6A3D9A" "#FFFF99" "#B15928"
labels <- c(
  # Contig solo
  "Blast" = "TrinityBlast",
  "SpadesBlast" = "SpadesBlast",
  
  # Kmer solo
  "kuniq"= "kuniq",
  "kraken2"= "kraken2",
  
  # Hybrid Trinity
  "hybrid_Blast-kraken2"="hybrid_TrinityBlast-kraken2",
  "hybrid_Blast-kuniq"="hybrid_TrinityBlast-kuniq",
  
  # Hybrid Spades
  "hybrid_SpadesBlast-kraken2"="hybrid_SpadesBlast-kraken2",
  "hybrid_SpadesBlast-kuniq"="hybrid_SpadesBlast-kuniq",    
  
  ### Combined Kmer
  "hybrid_kraken2-kuniq"="k2-ku",
  "hybrid_kuniq-kraken2"="ku-k2" )         

labels2lines <- c(
  # Contig solo
  "TrinityBlast" = "- Trinity-Blast -\nnt ",
  "SpadesBlast" = "- SPAdes-Blast -\nnt ",
  
  # Kmer solo
  "kuniq" = "- KUniq -\nmicrobialDB",
  "kraken2" = "- K2 -\nnt",     
  
  # Hybrid Trinity
  "hybrid_TrinityBlast-kraken2" = "- Hybrid-Trinity-K2 -\nnt",
  "hybrid_TrinityBlast-kuniq" = "- Hybrid-Trinity-KUniq -\nnt + microbialDB",
  
  # Hybrid spades
  "hybrid_SpadesBlast-kraken2" =   "- Hybrid-SPAdes-K2 -\nnt",
  "hybrid_SpadesBlast-kuniq"=  "- Hybrid-SPAdes-KUniq -\nnt + microbialDB",
  
  ### Combined Kmer
  "ku-k2" = "- KUniq-K2 -\nmicrobialDB + nt",
  "k2-ku" = "- K2-KUniq -\nnt + microbialDB")

palette <- c(
  "TrinityBlast" = "#E3D7EB",      # Blue
  "SpadesBlast" = "#FDD9A0",

  "kraken2" = "#FB9A99",    # Bright Red (unchanged)
  "kuniq" = "#E31A1C", # Deep, bold red
  
  
  "ku-k2" = "#33A02C",
  "k2-ku" = "#B2DF8A",
  
  "hybrid_SpadesBlast-kraken2" = "#F19737",
  "hybrid_SpadesBlast-kuniq"= "#E66F00",
  
  "hybrid_TrinityBlast-kraken2" = "#CAB2D6",  # Lighter Green
  "hybrid_TrinityBlast-kuniq" = "#6A3D9A"  # Darker, richer green

)

color <- c(
  "TrinityBlast" = "#4B006E",      # Blue
  "SpadesBlast" = "#5e1e02",
  
  "kraken2" = "#730F0F",    # Bright Red (unchanged)
  "kuniq" = "#400808", # Deep, bold red
  
  
  "ku-k2" = "#112904",
  "k2-ku" = "#112904",
  
  "hybrid_SpadesBlast-kraken2" = "#5e1e02",
  "hybrid_SpadesBlast-kuniq"= "black",
  
  "hybrid_TrinityBlast-kraken2" = "#4B006E",  
  "hybrid_TrinityBlast-kuniq" = "black"
  
)

label_df=as.data.frame(labels2lines)
label_df$method=rownames(label_df)
colnames(label_df) = c("labels2lines","method_label")

label_key=as.data.frame(labels)
label_key$method=rownames(label_key)
colnames(label_key) = c("method_label","method")

label_df$family_method = c("Assembly-based","Assembly-based",
                            "Assembly-free","Assembly-free",
                            "Hybrid","Hybrid","Hybrid","Hybrid",
                            "Combined assembly-free", "Combined assembly-free")

label_df=merge(label_key,label_df, by="method_label")



label_df$family_method = factor(label_df$family_method, levels = c("Assembly-based", "Assembly-free" ,"Hybrid", "Combined assembly-free") )

label_df$labels2lines = factor(label_df$labels2lines, levels = c("- Trinity-Blast -\nnt " ,
                                                                 "- SPAdes-Blast -\nnt " ,
                                                                 "- K2 -\nnt",
                                                                 "- KUniq -\nmicrobialDB"  ,
                                                                 "- Hybrid-Trinity-K2 -\nnt",
                                                                 "- Hybrid-SPAdes-K2 -\nnt",
                                                                 "- Hybrid-Trinity-KUniq -\nnt + microbialDB",
                                                                 "- Hybrid-SPAdes-KUniq -\nnt + microbialDB" ,
                                                                 "- K2-KUniq -\nnt + microbialDB",
                                                                 "- KUniq-K2 -\nmicrobialDB + nt"))


theme_fig2=function(plot){
  return(plot + theme(legend.direction = "horizontal",legend.position = "bottom",
                      legend.key.size = unit(1, "cm"),
                      legend.key.width = unit(1,"cm"),
                      strip.background = element_rect(fill="#333333") , 
                      strip.text = element_text(size = 12, face = "bold"),
                      axis.title.y= element_text(size = 12, face = "bold"),
                      axis.title.x= element_text(size = 12, face = "bold"),
                      axis.text.x = element_text(size = 12),
                      axis.text.y = element_text(size = 10 )))
}




get_metrics_right_format=function(metrics_list, tax_lvl, method, rt){
  metrics_phylum=lapply(metrics_list, function(x){return(x[[tax_lvl]])})
  metrics_phylum=ldply(metrics_phylum)
  colnames(metrics_phylum)[1] = "replicate_id"
  metrics_phylum$method = method
  metrics_phylum$rt = rt
  
  return(metrics_phylum)  
}

get_metrics_right_format_all_ranks=function(metrics_list, method, rt){
  ranks=c("superkingdom","phylum","class","order","family","genus")
  
  metrics_list=lapply(ranks, function(lvl){get_metrics_right_format(metrics_list, lvl, method, rt)} )
  names(metrics_list)=ranks
  metrics_list=ldply(metrics_list)
  colnames(metrics_list)[1]="taxlvl"
  return(metrics_list)
}


################################################################################
# V A R I A B L E S
################################################################################

path="~/"

#path = path_local

#setwd("/home/acolajanni/Documents/work/")

# simulation_letter = args[1]
# path_results = args[2]
# strategy = args[3]

path_results=paste0(path, "results/Simulation/with_replacement/kraken/nt/Kraken2/")


sra_id=c(0:49)
#sra_id=sra_id[str_detect(sra_id, simulation_letter)]

################################################################################
# A L G O R I T H M 
################################################################################
replicates=c(0:49)

replicate_letters=c("")


# Retrieving read names to get classification of each to get the diversity
path_reads=paste0(path,"data/Simulation/realist/raw_reads/")


classifs=parallel::mclapply(replicates, function(sra){
  print(sra)
  path_r=paste0(path_reads,sra,"/",sra,"_ReadNames.txt")
  r=readLines(path_r)
  tmp=as.data.frame(list("read"=r, "taxid"=sub(".*taxid:([0-9]+)_.*", "\\1", r)))
  tmp=aggregate(read ~ taxid, tmp, length)
  return(tmp)
},mc.cores=ncores)
names(classifs)=replicates

### Merge with truth to retrieve full classif- retrieve dataframe of read origin
load(paste0(path,"/database_clean/selected_genomes/realist/selected_genomes_list.RData"))
selected_genomes=unique(ldply(random_drawing))
selected_genomes$.id=NULL

truth=selected_genomes[,c(1,3,16,17,18,19,20,21,22)]
colnames(truth) = c('taxid',"species_taxid",ranks)
# Add human into the reads
truth=rbind(truth, as.data.frame(list(9606,9606, "Eukaryota", "Chordata", "Mammalia",
                                      "Primates", "Hominidae", "Homo", "Homo sapiens"),
                                 col.names = colnames(truth) ))
truth=unique(truth)
truth$phylum = str_replace_all(truth$phylum, pattern = "-", replacement = " ")

# 3.2 Merge to get actual classification
truth_vir = truth[truth$superkingdom == "Viruses",]  
truth = truth[! truth$superkingdom == "Viruses",]
truth = truth %>% mutate(across(
  c("superkingdom","phylum", "class", "order", "family", "genus", "species"), ~ gsub("^unclassified.*", "unclassified", .)))

truth=rbind(truth, truth_vir)


### Retrieve metrics dataframes
path_results=paste0(path,"results/Simulation/realist/")

letters_replicates=""

metrics_list=list()
conf_mat_list=list()
classif_list=list()
conf_mat_list_genus=list()

ranks=c("superkingdom","phylum","class","order","family","genus")
path_results=paste0(path,"results/Simulation/realist/pipeline_evaluation/")

path_classif_df=paste0(path,"results/Simulation/realist/")

methods=list.files(path=path_results)



for(i in letters_replicates){
  print(i)
  if       (i=="A"){rt="5r/t"
  }else if (i %in% c("C","M","O","Q") ){rt="10r/t"
  }else if (i %in% c("E","N","P","R")){rt="100r/t"
  }else if (i == ""){rt=NA}
  
  
  tmp_metrics = list()
  tmp_matrix  = list()
  for (m in methods){
    
    if (i %in% c("A","C","E") & str_detect(m , "human_kuniq") ){
      next
    }
    
    if        (m == "Blast")        { meth="Contigs"
    } else if (m == "SpadesBlast")  { meth="Contigs_rnaSpades"
    } else if (m == "hybrid_SpadesBlast-kraken2"){ meth = "hybrid_spadesBlast-kraken2"
    } else if (m == "hybrid_SpadesBlast-kuniq"){   meth = "hybrid_spadesBlast-kuniq"
    } else if (m == "kraken2"){ meth="kraken/nt/k2uniq"
    } else if (m == "kuniq"){ meth="kraken/microbialDB/kuniq"
    } else { meth=m }

    
    if (str_detect(m, "hybrid")){
      classif_df_file=paste0(path_classif_df,"hybrid/",meth,"/classif_df_genus.csv")
    }else{
      classif_df_file=paste0(path_classif_df,meth,"/classif_df_genus.csv")
    }
    
    
    
    metrics_file=paste0(path_results, m,"/Simulation_metrics_",m,"_",i,".RData")
    conf_mat_file=paste0(path_results,m,"/Confusion_matrix_", m,"_",i, ".tsv")
    conf_mat_list_genus_file= paste0(path_results,m,"/Confusion_matrix_genus_", m,"_",i, ".tsv")

    if ( file.exists(metrics_file) ) { # if you file.resists prouves que tu file.exists
      if ( file.info(metrics_file)$size >= 1000 ) {
        print(m)
        load(file=metrics_file)
        metrics=get_metrics_right_format_all_ranks(metrics_full, m , rt)
        tmp_metrics[[m]]= metrics
        rm(metrics_full)
        
        ### Load classif matrix
        x=read.table(classif_df_file, header=TRUE, sep = ",")
        classif_list[[m]] = x
        
        ### Load confusion-genus
        x=unique(read.table(conf_mat_list_genus_file, header=TRUE, sep = "\t"))
        #colnames(x) = c("truth","pred","Freq","total","Prop")
        conf_mat_list_genus[[m]] = x
        
        ### Load Confusion matrix
        x=unique(read.table(conf_mat_file, header=TRUE, sep = "\t"))
        x$method=m
        n_true=length(unique(x$taxon_full_name_truth))
        t=as.data.frame(table(x$taxon_full_name))
        problems=as.character(t[t$Freq > n_true,]$Var1)
        ### Handle duplicate entries
        for (p in problems){
          tmp=x[x$taxon_full_name == p, ]
          duplicated_values=tmp[duplicated(tmp$taxon_full_name_truth) , ]$taxon_full_name_truth
          dup=tmp[ tmp$taxon_full_name_truth %in% duplicated_values, ]
          tmp=tmp[!tmp$taxon_full_name_truth %in% duplicated_values, ]
          dup=dup[dup$Freq > 0, ]
          tmp=rbind(tmp,dup)
          
          x=x[x$taxon_full_name != p, ]
          x=rbind(x,tmp)
        }
        tmp_matrix[[m]] = x
      }
    }
  }
  
  metrics_list[[i]]  = do.call(rbind, tmp_metrics)
  conf_mat_list[[i]] = do.call(rbind, tmp_matrix)
  
}
#### Test
truth_total_per_sample=ldply(lapply(classifs, function(x) sum(x$read)) )
colnames(truth_total_per_sample) = c(".id","total_truth")

total_per_sample =read.table (paste0(path,"results/Simulation/realist/Contigs/Contigs_n_readsClassifs.txt"),quote="\"", comment.char="",
                              col.names = c("total_pred",".id"))[0:50,]


total_per_sample$.id =sub("^\\./([^/]+)/.*", "\\1", total_per_sample$.id)
total_per_sample = merge(total_per_sample, truth_total_per_sample, by=".id")
total_per_sample$diff=total_per_sample$total_truth - total_per_sample$total_pred
####
classifs_pair=conf_mat_list[[1]]
classifs_pair$prop=classifs_pair$prop/100


classifs_pair$classification = ifelse(classifs_pair$taxon_full_name_truth == classifs_pair$taxon_full_name, "TP", "FP")
classifs_pair$classification = ifelse(classifs_pair$taxon_full_name %in% c("Unassembled", "Unclassified"), "U", classifs_pair$classification)
classifs_pair$classification = factor(classifs_pair$classification, levels= rev(c("TP","FP","U")))

classifs_pair$taxon_full_name_truth = ifelse(classifs_pair$taxon_full_name_truth == "Eukaryota - Chordata",
                                             'Eukaryota - Homo sapiens', classifs_pair$taxon_full_name_truth)


setdiff(unique(classifs_pair$method) , unique(label_df$method))

unique(classifs_pair$method)
unique(label_df$method)

classifs_pair = merge(label_df, classifs_pair, by="method")
classifs_pair$family_method = factor(classifs_pair$family_method, levels=c("Assembly-based","Assembly-free", "Hybrid", "Combined assembly-free"))



################################################################################
## TP/FP/U per phylum as barplot

theme_fig2=function(plot){
  return(plot + theme(legend.direction = "horizontal",legend.position = "bottom",
                      legend.key.size = unit(1, "cm"),
                      legend.key.width = unit(1,"cm"),
                      legend.text = element_text(size = 8, face = "bold"),
                      strip.background = element_rect(fill="#333333") , 
                      strip.text = element_text(size = 7, face = "bold"),
                      axis.title.y= element_text(size = 9, face = "bold"),
                      axis.title.x= element_text(size = 9, face = "bold"),
                      axis.text.x = element_text(size = 8),
                      axis.text.y = element_text(size = 8 )))
}



p=ggplot(classifs_pair, aes(x=taxon_full_name_truth, y=prop, fill = classification)) + 
  geom_bar(position="stack", stat="identity") + 
  theme_linedraw()+
  xlab("Read simulated per transcripts") + ylab("Proportion")+
  scale_fill_manual(values=c("U"="#FFB462", "FP"= "#B12424","TP"="#7AADD0"),
                    labels = c("Unclassified", "False prediction", "True Prediction"),
                    name="" )+
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(1, "cm"),
        legend.key.width = unit(1,"cm"),
        legend.text = element_text(size = 11, face = "bold"),
        strip.text = element_text(size = 10, face = "bold"),
        axis.title.y= element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 10, angle=75, hjust=1),
        axis.text.y = element_text(size = 10 ),
        
  ) +  facet_nested_wrap( vars(family_method, labels2lines), nrow=2)


p0=theme_fig2(p)
p0



## TP/FP/U aggregated by method: 1method + 1taxon = 3 valeurs: TP/FP/U (= 10 methods, 11taxons = 330 rows )
agg_classif_pair=aggregate(prop ~ classification + method + method_label + family_method + labels2lines + taxon_full_name_truth,  classifs_pair, sum )

p=ggplot(agg_classif_pair, aes(x=labels2lines, y=prop, fill = classification)) + 
  geom_boxplot()+
  theme_linedraw()+
  xlab(" ") + 
  ylab("Proportion False Prediction")+
  scale_fill_manual(values=c("U"="#FFB462", "FP"= "#B12424","TP"="#7AADD0"),
                    labels = c("Unclassified", "False prediction", "True Prediction"),
                    name="" )+
  facet_nested(
    rows = vars(classification),
    cols = vars(family_method, labels2lines),
    scales = "free"
  )
p1=theme_fig2(p)

### peut etre un jitter avec des points en transparence ?
p=ggplot(agg_classif_pair, aes(x=labels2lines, y=prop, fill = classification)) + 
  geom_boxplot()+
  theme_linedraw()+
  xlab(" ") + 
  ylab("Classification proportion at phylum level ")+
  scale_fill_manual(values=c("U"="#FFB462", "FP"= "#B12424","TP"="#7AADD0"),
                    labels = c("Unclassified", "False prediction", "True Prediction"),
                    name="" )+
  facet_nested_wrap( vars(family_method, labels2lines), scales='free_x', nrow=1 ) 

p2=theme_fig2(p) 
p2 = p2+ theme(axis.text.x=element_blank())


p2

ggsave(p2,filename=paste0("~/results/Simulation/realist/figures/FigSimuRealist1.jpeg"),
       device = "jpeg",width = 20, height = 10,dpi = 350 , units = "cm", create.dir = TRUE)



################################################################################
## TP/FP/U per genus as boxplot
################################################################################

genus_pair=ldply(conf_mat_list_genus)
colnames(genus_pair)[1] = "method"

hs="Eukaryota_Chordata_Mammalia_Primates_Hominidae_Homo"

# Set mammalia prediction as Human
genus_pair$taxon_full_name_truth = ifelse(genus_pair$taxon_full_name_truth == hs, "Homo sapiens",genus_pair$taxon_full_name_truth)
genus_pair$taxon_full_name = ifelse(str_detect(genus_pair$taxon_full_name,"Mammalia"), "Homo sapiens",genus_pair$taxon_full_name)

genus_pair$classification = ifelse(genus_pair$taxon_full_name_truth == genus_pair$taxon_full_name, "TP", "FP")
genus_pair$classification = ifelse(genus_pair$taxon_full_name %in% c("unassembled", "unclassified"), "U", genus_pair$classification)

genus_pair$classification = factor(genus_pair$classification, levels= rev(c("TP","FP","U")))
total_df=unique(genus_pair[,c("taxon_full_name_truth", "total")])
genus_pair_agg=aggregate(Freq ~ taxon_full_name_truth + method + classification, genus_pair, sum )
genus_pair_agg = merge(genus_pair_agg, total_df, by="taxon_full_name_truth")
genus_pair_agg$prop = genus_pair_agg$Freq / genus_pair_agg$total


genus_pair_agg = merge(label_df, genus_pair_agg, by="method")
genus_pair_agg$family_method = factor(genus_pair_agg$family_method, levels=c("Assembly-based","Assembly-free", "Hybrid", "Combined assembly-free"))

genus_pair_agg$sk = str_split_fixed(genus_pair_agg$taxon_full_name_truth, fixed("_"), n = 6)[,1]
genus_pair_agg$sk = ifelse(genus_pair_agg$sk != "Homo sapiens", "Microbial", genus_pair_agg$sk)

p=ggplot(genus_pair_agg[!str_detect(genus_pair_agg$method_label, "kraken2") , ], aes(x=labels2lines, y=prop, fill = method_label)) + 
  geom_boxplot()+
  theme_linedraw()+
  xlab(" ") + 
  ylab("Classification proportion at genus level ")+
  # scale_fill_manual(values=c("U"="#FFB462", "FP"= "#B12424","TP"="#7AADD0"),
  #                   labels = c("Unclassified", "False prediction", "True Prediction"),
  #                   name="" )+
  scale_fill_manual(values=palette) +
  scale_color_manual(values=color) +
  facet_nested(cols=vars(family_method, labels2lines),
               rows=vars(classification) ,scales='free' , 
               strip = strip_nested(size = "variable")) 

p2=theme_fig2(p) 
p2 = p2+ theme(axis.text.x=element_blank())


p2




################################################################################
# Classif metrics
################################################################################

metrics=metrics_list[[1]]
metrics=merge(label_df, metrics, by="method")


unique(metrics$method)
unique(label_df$method)


metrics_g=metrics[metrics$taxlvl=="genus" , ]
metrics_g$superkingdom = sapply(metrics_g$taxon, function(x) strsplit(x, "_")[[1]][1])
metrics_g$phylum = sapply(metrics_g$taxon, function(x) strsplit(x, "_")[[1]][2])
metrics_g$genus = sapply(metrics_g$taxon, function(x) strsplit(x, "_")[[1]][6])

metrics_g$human = ifelse(metrics_g$genus == "Homo", "Human", "Microbial") 

metrics_g$family_method = factor(metrics_g$family_method, levels=c("Assembly-based","Assembly-free", "Hybrid", "Combined assembly-free"))



p=ggplot(metrics_g[!str_detect(metrics_g$method_label, "kraken2") , ], aes(x=method_label, y=F1, fill = method_label,color=method_label)) + 
  geom_boxplot()+
  theme_linedraw() + 
  xlab(" ") + 
  ylab("F1-scores at genus level")+
  scale_fill_manual(values=palette) + 
  scale_color_manual(values=color) + 
  
  facet_nested( cols = vars(family_method, labels2lines), 
                rows = vars(superkingdom), scales='free_x') 

p3 = theme_fig2(p) + theme(axis.text.x=element_blank())
p3

################################################################################ matrix truth vs predicted

### Load true composition: transcripts
path_truth=paste0(path,"data/Simulation/realist/")
truth_classif=read.table(paste0(path_truth,"genus_to_select_lineage.txt"),header=FALSE, sep="\t")
colnames(truth_classif)=c("G","taxid","K","P","C","O","F")
truth_classif$genus = paste0("g_",truth_classif$G)

truth_df=read.table(paste0(path_truth,"transcript_to_draw.tsv"),header=TRUE)
truth_df=truth_df*100
row.names(truth_df) = as.character(c(0:49))

long_truth <- truth_df %>%
  mutate(.id = rownames(truth_df)) %>%
  pivot_longer(
    cols = starts_with("g_"),
    names_to = "genus",
    values_to = "count" )

long_truth=merge(truth_classif, long_truth, by="genus")
long_truth$taxid = NULL

### Load predicted composition
classif_dfs=parallel::mclapply(classif_list,function(df) {
  
  df$genus=sapply(df$classification, function(c) paste0("g_",str_split(c, fixed("|") )[[1]][6] ) )
  df$superkingdom=sapply(df$classification, function(c) str_split(c, fixed("|") )[[1]][1]  )
  df <- df %>%
    separate(
      classification,
      into = c("K", "P", "C", "O", "F", "G"),
      sep = "\\|",
      remove = FALSE,
      fill = "right" )
  return(df)
  }, mc.cores = ncores)

### Retrive entire taxonomy
taxonomy_dfs=parallel::mclapply(classif_dfs,function(df) {
  #df$G = paste0("g_",df$G)
  df=df[,c("K","P","C","O","F","G")]
  return(unique(df))
}, mc.cores = ncores)

taxonomy_df=unique(do.call(rbind, taxonomy_dfs))
taxonomy_df=unique(rbind(taxonomy_df, truth_classif[c("K","P","C","O","F","G")]))
taxonomy_df$OTU = paste0("OTU",1:nrow(taxonomy_df))

taxonomy_df$classification <- paste(
  taxonomy_df$K,taxonomy_df$P,
  taxonomy_df$C,taxonomy_df$O,
  taxonomy_df$F,taxonomy_df$G, 
  sep = "|" )

classif_dfs$truth=long_truth
summarized_df=parallel::mclapply(classif_dfs,function(df) {
  df$G = paste0("g_",df$G)
  df <- df %>%
    mutate( across( c(K, P, C, O, F),
        ~ case_when(
          G == "g_Unclassified" ~ "Unclassified",
          G == "g_Unassembled"  ~ "Unassembled",
          TRUE                  ~ . ) ) )
  df=aggregate(count ~ K + P + C + O + F + G + .id, df, sum)
  return(df)
}, mc.cores = ncores)


### Compare number of predicted genus in viruses / bacteria / eukaryota ==> diff entre vrai et prédit

n_genus_list=lapply(summarized_df, function(df) {
  df=df[,c("K","G",".id")]
  colnames(df) = c("K","G","sample_id")
  df = df[!df$G %in% c("g_unassembled","g_Unclassified"),]
  count_genus=aggregate(G ~ K + sample_id, df, length)
  # Complete the data
  count_genus <- count_genus %>%
    complete(sample_id, K = c("Bacteria", "Eukaryota", "Viruses", "Archaea"),
             fill = list(G = 0)) %>%
    arrange(sample_id, K)
  return(count_genus)
  
  })

truth_count=n_genus_list$truth
colnames(truth_count) = c(".id","K","true_G")
n_genus_list$truth = NULL

n_genus=ldply(n_genus_list)
colnames(n_genus) = c("method", ".id", "K","G")

n_genus = merge(n_genus, truth_count, by=c(".id","K"))
n_genus$diff = n_genus$G - n_genus$true_G


n_genus=merge(label_df, n_genus, by="method")
n_genus$family_method = factor(n_genus$family_method, levels=c("Assembly-based","Assembly-free", "Hybrid", "Combined assembly-free"))



p=ggplot(n_genus[!str_detect(n_genus$method_label, "kraken2|k2-ku") , ], aes(x=method_label, y=diff, fill = method_label,color=method_label)) + 
  geom_hline(yintercept=0, linetype="dashed", 
             color = "darkred", size=1)+
  geom_boxplot()+
  theme_linedraw() + 
  xlab(" ") + 
  ylab("Number of genus over/underprediction")+
  scale_fill_manual(values=palette) + 
  scale_color_manual(values=color) + 
  
  facet_nested( cols = vars(family_method, labels2lines), 
                rows = vars(K), scales='free') 

p4 = theme_fig2(p) + theme(axis.text.x=element_blank())
p4


### Compare frobenius after removing homo sapiens and unclassified ==> in proportion

wide_classif=parallel::mclapply(classif_dfs,function(df) {
  df$.id = as.character(df$.id)
  df=df[,c("G",".id","count")]
  wide_df <- as.data.frame.matrix(xtabs(count ~ .id + G, data = df))  # creates a matrix                   
  return(wide_df) }, mc.cores=ncores)

# List all genus 

colnames(truth_df) = sapply(colnames(truth_df), function(x) str_remove(x, 'g_') )
wide_classif$truth = truth_df

all_genus=lapply(wide_classif, function(x) colnames(x))
all_genus$truth = colnames(truth_df)
all_genus = unique(unlist(all_genus))


# Add these empty columns to each dataframe
classif_samecol <- parallel::mclapply(wide_classif, function(df) {
  missing_cols <- setdiff(all_genus, names(df))
  
  for (col in missing_cols) {
    df[[col]] <- 0 }
  return(df[,all_genus]) }, 
  mc.cores=ncores)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### Sample Wise frobenius/MAE 

# Compute both metrics in one go
distances_raw <- do.call(rbind, parallel::mclapply(names(classif_samecol)[names(classif_samecol) != "truth"], function(method) {
  truth <- as.matrix(classif_samecol$truth)
  pred  <- as.matrix(classif_samecol[[method]])
  
  # Frobenius row-wise
  frob <- sqrt(rowSums((truth - pred)^2))
  # MAE row-wise
  mae <- rowMeans(abs(truth - pred))
  
  # Combine into a dataframe
  data.frame(
    method  = rep(method, length(frob)),
    value = c(frob, mae),
    metric = rep(c("frobenius", "MAE"), each = length(frob)) )
}, mc.cores = ncores))


distances_raw=merge(label_df, distances_raw, by="method")
distances_raw$family_method = factor(distances_raw$family_method, levels=c("Assembly-based","Assembly-free", "Hybrid", "Combined assembly-free"))

p=ggplot(distances_raw[!str_detect(distances_raw$method_label, "kraken2|k2-ku") , ] , aes(x=method_label, y=value, fill = method_label,color=method_label)) + 
  geom_boxplot()+
  theme_linedraw() + 
  xlab(" ") + 
  ylab("distance")+
  scale_fill_manual(values=palette) + 
  scale_color_manual(values=color) + 
  
  facet_nested( cols = vars(family_method, labels2lines), 
                rows = vars(metric), scales='free') 

p5 = theme_fig2(p) + theme(axis.text.x=element_blank())
p5


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### Sample Wise frobenius/MAE After proportion + remove unclassified + human
distances <- do.call(rbind, parallel::mclapply(names(classif_samecol)[names(classif_samecol) != "truth"], function(method) {
  
  truth <- classif_samecol$truth
  pred  <- classif_samecol[[method]]
  
  pred[,c("Unclassified","unassembled","Homo")] = NULL
  truth[,c("Unclassified","unassembled","Homo")] = NULL
  
  pred=prop.table(as.matrix(pred), margin=1)
  truth=prop.table(as.matrix(truth), margin=1)
  
  # Frobenius row-wise
  frob <- sqrt(rowSums((truth - pred)^2))
  # MAE row-wise
  mae <- rowMeans(abs(truth - pred))
  
  # Combine into a dataframe
  data.frame(
    method  = rep(method, length(frob)),
    value = c(frob, mae),
    metric = rep(c("frobenius", "MAE"), each = length(frob)) )
}, mc.cores = ncores))


distances=merge(label_df, distances, by="method")
distances$family_method = factor(distances$family_method, levels=c("Assembly-based","Assembly-free", "Hybrid", "Combined assembly-free"))

p=ggplot(distances[!str_detect(distances$method_label, "kraken2|k2-ku") , ], aes(x=method_label, y=value, fill = method_label,color=method_label)) + 
  geom_boxplot()+
  theme_linedraw() + 
  xlab(" ") + 
  ylab("Number of genus over/underprediction")+
  scale_fill_manual(values=palette) + 
  scale_color_manual(values=color) + 
  
  facet_nested( cols = vars(family_method, labels2lines), 
                rows = vars(metric), scales='free') 

p6 = theme_fig2(p) + theme(axis.text.x=element_blank())
p6




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### Sample Wise frobenius/MAE After proportion + remove unclassified + human 
### ### ### ### ### ### ! Filtrating rare taxons !
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
method="kuniq"

filter_df=function(df, nzero_threshold=45, count_threshold=100){
  
  # Keep if at least 5 samples have this prediction
  count_zeroes=apply(df, MARGIN = 2, function(x){length(x[x==0]) })
  df=df[,names(count_zeroes[count_zeroes <= nzero_threshold])]
  
  # Keep if at least 100 count across the dataset for each column
  count_total=colSums(df)
  df=df[,names(count_total[count_total >= count_threshold])]
  return(df)
}

distances <- parallel::mclapply(names(classif_samecol)[names(classif_samecol) != "truth"], function(method) {
  
  truth <- classif_samecol$truth
  pred  <- classif_samecol[[method]]
  
  pred[,c("Unclassified","unassembled","Homo")] = NULL
  truth[,c("Unclassified","unassembled","Homo")] = NULL
  
  pred=filter_df(pred, 49 , 100)
  truth=filter_df(truth, 49 , 100)

  #---column homogeneisation
  cols=unique(c(colnames(pred),colnames(truth)))
  
  # predicted column not in truth ### false positives
  FP=colnames(pred)[! colnames(pred) %in% colnames(truth)]
  
  # true column not in predicted ### missing column // FN
  missing=colnames(truth)[! colnames(truth) %in% colnames(pred)]
  
  # true column (TP)
  true_col=colnames(truth)[ colnames(truth) %in% colnames(pred)]
  
  truth[,FP]=0
  pred[,missing]=0
  
  pred=prop.table(as.matrix(pred[,cols]), margin=1)
  truth=prop.table(as.matrix(truth[,cols]), margin=1)
  
  frobenius_dist <- norm(truth - pred, type = "F")
  
  # FP_list[[method]]=FP
  # missing_list[[method]]=missing
  
  # Frobenius row-wise
  frob1 = rowMeans(sqrt((truth - pred)^2))
  
  frob2 <- sqrt(rowSums((truth - pred)^2))
  # MAE row-wise
  mae <- rowMeans(abs(truth - pred))
  
  # Combine into a dataframe
  res=data.frame(
    method  = rep(method, length(frob1)),
    value = c(frob1, mae),
    metric = rep(c("Mean euclidean distance", "MAE"), each = length(frob1)) )
  # return(res)
  return(list("res"=res, "FP"=FP, "FN"=missing, "TP"=true_col, "frobenius"=frobenius_dist))
  
}, mc.cores = ncores)


names(distances)=names(classif_samecol)[names(classif_samecol) != "truth"]

dist_list=lapply(distances, function(x) x$res)
FP_list=lapply(distances, function(x) length(x$FP) )
missing_list=lapply(distances, function(x) length(x$FN))
true_list=lapply(distances, function(x) length(x$TP))
frobenius_list=lapply(distances, function(x) x$frobenius)


distances=do.call(rbind, dist_list)

distances=merge(label_df, distances, by="method")
distances$family_method = factor(distances$family_method, levels=c("Assembly-based","Assembly-free", "Hybrid", "Combined assembly-free"))

p=ggplot(distances[!str_detect(distances$method_label, "kraken2|k2-ku") , ], 
         aes(x=method_label, y=value, fill = method_label,color=method_label)) + 
  geom_boxplot()+
  theme_linedraw() + 
  xlab(" ") + 
  ylab("Number of genus over/underprediction")+
  scale_fill_manual(values=palette) + 
  scale_color_manual(values=color) + 
  
  facet_nested( cols = vars(family_method, labels2lines), 
                rows = vars(metric), scales='free') 

p7 = theme_fig2(p) + theme(axis.text.x=element_blank())
p7


FP_list=ldply(FP_list)
missing_list=ldply(missing_list)
TP_list=ldply(true_list)

colnames(FP_list) = c("method","value")
colnames(missing_list) = c("method","value")
colnames(TP_list) = c("method","value")

FP_list$metric = "FP"
missing_list$metric = "Missing"
TP_list$metric = "TP"

count_errors=do.call(rbind, list(TP_list,FP_list,missing_list))
count_errors$value = ifelse(count_errors$metric=="Missing", -count_errors$value, count_errors$value)




count_errors=merge(label_df, count_errors, by="method")
count_errors$family_method = factor(count_errors$family_method, levels=c("Assembly-based","Assembly-free", "Hybrid", "Combined assembly-free"))

count_errors2=count_errors[!str_detect(count_errors$method_label, "kraken2|k2-ku") , ]

true_taxon=colnames(truth_df)[colnames(truth_df) != "Homo"]

p=ggplot(count_errors2, aes(x = method, y = value, fill = metric)) +
  geom_col(position = "stack") +
  geom_hline(yintercept = length(true_taxon), 
             linetype = "dashed", 
             linewidth = 0.6)+
  theme_linedraw() + 
  xlab(" ") +
  geom_hline(yintercept = 0, color = "black") +
  ylab("Number of missing and falsely predicted genus")+
  scale_fill_manual(values=c(
    "TP"="#7AADD0",
    "FP"= "#B12424",
    "Missing"="#FFB462"),
    labels = c("Spurious genera (FP)" ,
               "Missing genera (FN)",
               "Correctly predicted genera (TP)"),
    name="" )+
  geom_text(
    data = subset(count_errors2, metric != "Missing"),
    aes(label = value),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  scale_y_continuous(
    breaks = function(x) sort(unique(c(pretty(x), 79))),
    labels = abs
  )+
  # Labels for Missing (placed below bar)
  geom_text(
    data = subset(count_errors2, metric == "Missing"),
    aes(label = abs(value)),
    vjust = 1.3,   # push downward
    size = 4
  ) +
  

  #scale_fill_manual(values=palette) + 
  #scale_color_manual(values=color) + 
  facet_nested( cols = vars(family_method, labels2lines), 
                #rows = vars(metric), 
                scales='free') 
#p
p8 = theme_fig2(p) + theme(axis.text.x=element_blank())
p8





### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### diff in prediction numbers, sample per samples
### ### ### ### ### ### ! Filtrating rare taxons - based on minimal count 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

get_pred_taxon=function(truth, pred, threshold=1 ){
  test=lapply(as.character(c(0:49)), function(index_sample) {
    
    tmp=as.data.frame(t(truth))
    true_tax=row.names(tmp[tmp[[index_sample]] >= threshold , ])
    
    tmp=as.data.frame(t(pred))
    pred_tax=row.names(tmp[tmp[[index_sample]] >= threshold , ])
    
    #---column homogeneisation
    cols=unique(c(true_tax,true_tax))
    # predicted column not in truth ### false positives
    FP=pred_tax[! pred_tax %in% true_tax]
    # true column not in predicted ### missing column // FN
    missing=true_tax[! true_tax %in% pred_tax]
    # true column (TP)
    true_col=true_tax[ true_tax %in% pred_tax]
    
    return(list("FP"=FP, "FN"=missing, "TP"=true_col))
  })
  return(test)
}


thresh_list=list()
thresh_vect=c(1,5,10,25,50,75,100)
for (t in thresh_vect){
  diff_genus <- parallel::mclapply(names(classif_samecol), function(method) {
    
    truth <- classif_samecol$truth
    pred  <- classif_samecol[[method]]
    
    pred[,c("Unclassified","unassembled","Homo")] = NULL
    truth[,c("Unclassified","unassembled","Homo")] = NULL
    
    tmp=get_pred_taxon(truth, pred, t)
    
    
    FP_list=as.data.frame(ldply( lapply(tmp, function(x) length(x$FP) )  ))
    FN_list=as.data.frame(ldply( lapply(tmp, function(x) length(x$FN) )  ))
    TP_list=as.data.frame(ldply( lapply(tmp, function(x) length(x$TP) )  ))
    
    FP_list$metric="FP"
    TP_list$metric="TP"
    FN_list$metric="FN"
    
    
    res=do.call(rbind, list(FP_list,TP_list,FN_list))
    
    return(res)
    
  }, mc.cores = 20)
  
  names(diff_genus) = names(classif_samecol)
  
  n_genus=ldply(diff_genus)
  colnames(n_genus) = c("method" , "value", "metric")
  
  n_genus$threshold=t
  
  thresh_list[[as.character(t)]] = n_genus
}

genus_count=do.call(rbind, thresh_list)


label_df <- label_df %>%
  add_row(
    method_label = "ground_truth",
    method = "truth",
    family_method = " ",
    labels2lines = "- Ground Truth -"
  )
label_df=unique(label_df)

genus_count=merge(label_df, genus_count, by="method")
genus_count$family_method = factor(genus_count$family_method, levels=c(" ","Assembly-based","Assembly-free", "Hybrid", "Combined assembly-free"))


genus_count <- genus_count[genus_count$method != "truth" | genus_count$threshold == 1, ]

genus_count$threshold = ifelse(genus_count$method == "truth", "No filter", as.character(genus_count$threshold))
genus_count$threshold = factor(as.character(genus_count$threshold),levels=c('No filter' , as.character(thresh_vect) ))


genus_count2=genus_count[!str_detect(genus_count$method_label, "kraken2|k2-ku") , ]



palette["ground_truth"] = "white"
color["ground_truth"] = "black"


genus_count2$metric_label=ifelse(genus_count2$metric == "FN", "Missing genus" , genus_count2$method)
genus_count2$metric_label=ifelse(genus_count2$metric == "FP", "False positive genus" , genus_count2$metric_label)
genus_count2$metric_label=ifelse(genus_count2$metric == "TP", "True positive genus" , genus_count2$metric_label)
genus_count2$metric_label = factor(genus_count2$metric_label, levels=c("True positive genus","False positive genus","Missing genus") )


genus_count2$labels2lines = factor(genus_count2$labels2lines, levels = c(
  "- Ground Truth -"  ,
  "- Trinity-Blast -\nnt " ,
  "- SPAdes-Blast -\nnt " ,
  "- KUniq -\nmicrobialDB"  ,
  "- Hybrid-Trinity-KUniq -\nnt + microbialDB",
  "- Hybrid-SPAdes-KUniq -\nnt + microbialDB" ,
  "- KUniq-K2 -\nmicrobialDB + nt"))


p=ggplot(genus_count2, aes(x=threshold, y=value, fill = method_label,color=method_label)) + 
  geom_boxplot()+
  theme_linedraw() + 
  xlab("Read threshold to keep a taxon") + 
  ylab("Number of predicted genus")+
  scale_fill_manual(values=palette) +
  scale_color_manual(values=color) +
  guides(fill="none" , color="none")+
  # facet_wrap(vars(metric), scales="free")  
  facet_nested( cols = vars(family_method, labels2lines),
                rows = vars(metric_label), scales='free')



p9 = theme_fig2(p)
p9

genus_count_summary=aggregate(value ~ method + metric + threshold + family_method, genus_count2, median)




p=ggplot(genus_count2[!genus_count2$threshold %in% c("1","5"),], aes(x=threshold, y=value, fill = method_label,color=method_label)) + 
  geom_boxplot()+
  theme_linedraw() + 
  xlab("Read threshold to keep a taxon") + 
  ylab("Number of predicted genus")+
  scale_fill_manual(values=palette) +
  scale_color_manual(values=color) +
  guides(fill="none", color="none") +
  # facet_wrap(vars(metric), scales="free")  
  facet_nested( cols = vars(family_method, labels2lines),
                rows = vars(metric_label), scales='free')

p10 = theme_fig2(p)
p10

supp2jpeg="~/results/Simulation/figures_review/figsupp2.jpeg"
ggsave(p9,filename=supp2jpeg,
       device = "jpeg",width = 30, height = 20,dpi = 350 , units = "cm", create.dir = TRUE)

magick::image_read(supp2jpeg)

fig6jpeg="~/results/Simulation/figures_review/fig6.jpeg"
ggsave(p10,filename=fig6jpeg,
       device = "jpeg",width = 30, height = 20,dpi = 350 , units = "cm", create.dir = TRUE)


magick::image_read(fig6jpeg)


supp2jsvg="~/results/Simulation/figures_review/figsupp2.svg"
ggsave(p9,filename=supp2jsvg,
       device = "svg",width = 30, height = 20,dpi = 600 , units = "cm", create.dir = TRUE)

fig6svg="~/results/Simulation/figures_review/fig6.svg"
ggsave(p10,filename=fig6svg,
       device = "svg",width = 30, height = 20,dpi = 600 , units = "cm", create.dir = TRUE)











