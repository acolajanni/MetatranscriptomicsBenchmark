library(stringr)
library(plyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(tidyr)
library(patchwork) 
library(ggh4x)
library(ggtext)

################################################################################
# F U N C T I O N S
################################################################################
ranks=c("superkingdom","phylum", "class", "order", "family", "genus", "species")


Get_alpha_div_per_samples=function(truth_df, classif_df){
  # Merging with real annotation
  alphadiv = merge(classif_df, truth_df[c("taxid", "superkingdom","phylum","class","order","family","genus","species")], by="taxid")
  # Alpha div per phylum
  taxonomic_levels=c("superkingdom","phylum","class","order","family","genus")
  
  results <- list()
  for (level in taxonomic_levels) {
    columns_to_aggregate=taxonomic_levels[1:which(level==taxonomic_levels)]
    # Aggregate species count for the current level
    agg_result <- aggregate(species ~ ., data = alphadiv[, c(level, "species")], length)

    if (!level=="superkingdom"){
      agg_result =  unique(merge(agg_result, truth_df[,columns_to_aggregate], by=level ))
      agg_result$taxon_full_name <- apply(agg_result[,columns_to_aggregate], 1, function(row) paste(row, collapse = "_"))
      agg_result=agg_result[,c(level,"species","taxon_full_name")]
    } else{
      colnames(agg_result)[1]="taxon"
      agg_result$taxon_full_name=agg_result$taxon
    }
    colnames(agg_result) = c("taxon","species","taxon_full_name")
    agg_result$taxlvl=level
    # Save the result in a list
    results[[level]] <- agg_result
  }
  alphadiv=do.call(rbind,results)
  return(alphadiv)
}



# Shanon: For each species i, sum( pi * log2(pi)
# Where pi : proportion of one specie in a sample // 
# Need: dataframe with n columns: taxid + phylum n_reads / taxid // need to calculate the total per phylum 
Get_indexes_phylum=function(truth_df, classif_list){
  
  shanon = parallel::mclapply(classif_list, function(x){
    indexes=merge(x, truth_df[,c("taxid","superkingdom","phylum")], by="taxid")
    indexes$taxon_full_name = paste0(indexes$superkingdom, " - ", indexes$phylum)
    
    indexes=indexes %>%
      group_by(taxon_full_name) %>%
      mutate(total_phylum = sum(read),
             alpha_div = length(taxid)) %>%
      select(- superkingdom)
    
    indexes$pi = indexes$read / indexes$total_phylum
    indexes$pilog2pi = indexes$pi * log2(indexes$pi)
    indexes$pisquare = indexes$pi * indexes$pi
    
    indexes=indexes %>%
      group_by(taxon_full_name) %>%
      mutate(shanon_index = -sum(pilog2pi),
             simpson_index = sum(pisquare))
    indexes=indexes[,c("shanon_index","simpson_index","alpha_div",'taxon_full_name',"phylum")]
    indexes$taxlvl="phylum"
    return(indexes)
  }, mc.cores=20)
  
  return(shanon) }


# Need: dataframe with n columns: taxid + phylum n_reads / taxid // need to calculate the total per phylum 
Get_diversity_indexes=function(truth_df, classif_list, tax_lvl="family"){
  
  shanon = parallel::mclapply(classif_list, function(x){
    indexes=merge(x, truth_df[,c("taxid","superkingdom", "phylum" , tax_lvl)], by="taxid")
    indexes$taxon_full_name = paste0(indexes$superkingdom, " - ", indexes[[tax_lvl]])
    
    indexes=aggregate(read ~ superkingdom + phylum + family + taxon_full_name, indexes, sum )
    indexes=indexes %>%
      group_by(phylum) %>%
      mutate(total_phylum = sum(read),
             alpha_div = length(family)) 
    
    indexes$pi = indexes$read / indexes$total_phylum
    indexes$pilog2pi = indexes$pi * log2(indexes$pi)
    indexes$pisquare = indexes$pi * indexes$pi
    
    indexes=indexes %>%
      group_by(phylum) %>%
      mutate(shanon_index = -sum(pilog2pi),
             simpson_index = sum(pisquare))
    
    indexes$taxon_full_name = paste0(indexes$superkingdom, " - ", indexes$phylum)

    indexes=indexes[,c("shanon_index","simpson_index","alpha_div",'taxon_full_name',"phylum")]
    indexes$taxlvl="phylum"
    return(indexes)
  }, mc.cores=20)
  
  return(shanon) }

### Merge across multiple samples to get all the existent taxid
### Then 
Get_alpha_div=function(truth_df, classif_list){
  # Get occurrences of each taxid
  classif_df=do.call(rbind, classif_list)
  classif_df=aggregate(read ~ taxid, classif_df, sum)
  
  alphadiv=Get_alpha_div_per_samples(truth_df, classif_df)
  return(alphadiv)
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


path_results=paste0(path, "results/Simulation/family/kraken/nt/Kraken2/")


sra_id=readLines(paste0(path, "data/Simulation/family/patient_name.txt"))
#sra_id=sra_id[str_detect(sra_id, simulation_letter)]

################################################################################
# A L G O R I T H M 
################################################################################
replicates=c(
  paste0("A",0:9),
  paste0("B",0:9),
  paste0("C",0:9),
  paste0("D",0:9))


replicate_letters=c("A","B","C","D")


# Retrieving read names to get classification of each to get the diversity
path_reads=paste0(path,"data/Simulation/family/raw_reads_10/")

#path_reads=paste0(path_local,"data/Simulation/family/")
classifs=parallel::mclapply(replicates, function(sra){
  print(sra)
  path_r=paste0(path_reads,sra,"/",sra,"_ReadNames.txt")
  r=readLines(path_r)
  tmp=as.data.frame(list("read"=r, "taxid"=sub(".*taxid:([0-9]+)_.*", "\\1", r)))
  tmp=aggregate(read ~ taxid, tmp, length)
  return(tmp)
},mc.cores=10)
names(classifs)=replicates
#save(classifs,file=paste0(path,"/data/Simulation/family/Read_taxid_count_MNOPQR.RData") )
#load(paste0(path,"/data/Simulation/family/Read_taxid_count_MNOPQR.RData"))

### Merge with truth to retrieve full classif- retrieve dataframe of read origin
load(paste0(path,"/database_clean/scripts/selected_genomes_list.RData"))
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

### Get read classification
alpha_div=lapply(classifs, function(x){Get_alpha_div_per_samples(truth, x)})
alpha_div=ldply(alpha_div)
alpha_div <- alpha_div %>%
  separate(taxon_full_name, into = c("superkingdom", "phylum"), sep = "_", remove = FALSE)
alpha_div$taxon=paste0(alpha_div$superkingdom,"_",alpha_div$phylum)
alpha_div$taxon <- sub("_NA$", "", alpha_div$taxon)

colnames(alpha_div)[3]="alpha_sample"
colnames(alpha_div)[1]="replicate_id"


# Alpha div global: For A, C and E
t=list(
  "A" = Get_alpha_div(truth, classifs[1:10]),
  "B" = Get_alpha_div(truth, classifs[11:20]),
  "C" = Get_alpha_div(truth, classifs[21:30]),
  "D" = Get_alpha_div(truth, classifs[31:40])
  )

alpha_div_global=ldply(t)
alpha_div <- alpha_div %>%
  separate(taxon_full_name, into = c("superkingdom", "phylum"), sep = "_", remove = FALSE)
colnames(alpha_div_global)[3]="alpha_global"


#### Retrive index
indexes=Get_diversity_indexes(truth, classifs, "family")
indexes=ldply(indexes)
colnames(indexes)[1]="replicate_id"

### Retrieve metrics dataframes
path_results=paste0(path,"results/Simulation/family/")

letters_replicates=c("A","B","C","D")

metrics_list=list()
metrics_phylum_list=list()
conf_mat_list=list()
ranks=c("superkingdom","phylum","class","order","family","genus")
path_results=paste0(path,"results/Simulation/family/pipeline_evaluation/")
methods=list.files(path=path_results)

for(i in letters_replicates){
  
  if       (i %in% c("A","C") ){rt="10r/t"
  }else if (i %in% c("B","D")){rt="100r/t"
  }

  
  tmp_metrics = list()
  tmp_matrix  = list()
  tmp_metrics_phylum = list()
  for (m in methods){
    
    if (i %in% c("A","C","E") & str_detect(m , "human_kuniq") ){
      next
    }
    
    metrics_file=paste0(path_results, m,"/Simulation_metrics_",m,"_",i,".RData")
    metrics_phylum_file=paste0(path_results, m,"/Simulation_metrics_phylum_",m,"_",i,".RData")
    conf_mat_file=paste0(path_results,m,"/Confusion_matrix_", m,"_",i, ".tsv")
    
    
    if ( file.exists(metrics_phylum_file) ) { # if you file.resists prouves que tu existes
      if ( file.info(metrics_phylum_file)$size >= 1000 ) {
        print(m)
        load(file=metrics_file)
        metrics=get_metrics_right_format_all_ranks(metrics_full, m , rt)
        tmp_metrics[[m]]= metrics
        rm(metrics_full)
        
        ### Load metrics aggregated at phylum level
        load(file=metrics_phylum_file)
        metrics_phylum=ldply(metrics_phylum)
        colnames(metrics_phylum)[1] = "replicate_id"
        colnames(metrics_phylum)[3] = "taxlvl"
        metrics_phylum$method = m
        metrics_phylum$rt=rt
        tmp_metrics_phylum[[m]]= metrics_phylum
        
        
        
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
  
  metrics_list[[i]]        = do.call(rbind, tmp_metrics)
  conf_mat_list[[i]]       = do.call(rbind, tmp_matrix)
  metrics_phylum_list[[i]] = do.call(rbind, tmp_metrics_phylum)
  
}

metrics_df=ldply(metrics_phylum_list)
metrics_df$taxlvl=factor(metrics_df$taxlvl, levels = ranks)


metrics_df$percentage_human=ifelse(!metrics_df$.id %in% c("A","B"), "90", "1/120" )



unique(metrics_df$method)




label_df <- as.data.frame(
  list(
    method = c(
      # Contig solo
      "Blast","SpadesBlast" ,
      # Hybrid Trinity
      "hybrid_Blast-kraken2", "hybrid_Blast-kuniq", 
      # Hybrid spades
      "hybrid_SpadesBlast-kraken2","hybrid_SpadesBlast-kuniq",
      # double kmer
      "hybrid_kraken2-kuniq", "hybrid_kuniq-kraken2", 
      # kmer solo
      "kraken2", "kuniq" ),
    
    labels2lines = c(
      # Contig solo
      "Trinity-Blast\nnt",
      "SPAdes-Blast\nnt",
      # Hybrid Trinity
      "Hybrid-Trinity-K2\nnt",
      "Hybrid-Trinity-KUniq\nnt + microbialDB",
      # Hybrid spades
      "Hybrid-SPAdes-K2\nnt",
      "Hybrid-SPAdes-KUniq\nnt + microbialDB",
      # double kmer
      "K2-KUniq\nnt+microbialDB",
      "KUniq-K2\nmicrobialDB+nt",
      # kmer solo
      "K2\nnt",
      "KUniq\nmicrobialDB" ),
    
    family_method = c(
      "Assembly-based","Assembly-based",
      "Hybrid", "Hybrid",
      "Hybrid", "Hybrid",
      "Combined assembly-free","Combined assembly-free",
      "Assembly-free","Assembly-free" ), 
    
    
    color=c("#053158","#084269", 
            "#6A0080","#4B006E",
            "saddlebrown","orangered4",
            "#275c09","#1F4A07",
            "#940C00","#730F0F" ) , 
    
    
    
    fill = c(
      "#0F4F8A", "#6FA9C9",
      "#9E84BE", "#4A1F78",
      "#D98E2F", "#C45500",
      "#7FB85C", "#1F6F1A",
      "#D95C5A", "#9E0E12" )
    ))


label_df <- as.data.frame(
  list(
    method = c(
      # Contig solo
      "Blast","SpadesBlast" ,
      # Hybrid Trinity
      "hybrid_Blast-kraken2", "hybrid_Blast-kuniq", 
      # Hybrid spades
      "hybrid_SpadesBlast-kraken2","hybrid_SpadesBlast-kuniq",
      # double kmer
      "hybrid_kraken2-kuniq", "hybrid_kuniq-kraken2", 
      # kmer solo
      "kraken2", "kuniq" ),
    
    labels2lines = c(
      # Contig solo
      "Trinity-Blast\nnt",
      "SPAdes-Blast\nnt",
      # Hybrid Trinity
      "Hybrid-Trinity-K2\nnt",
      "Hybrid-Trinity-KUniq\nnt + microbialDB",
      # Hybrid spades
      "Hybrid-SPAdes-K2\nnt",
      "Hybrid-SPAdes-KUniq\nnt + microbialDB",
      # double kmer
      "K2-KUniq\nnt+microbialDB",
      "KUniq-K2\nmicrobialDB+nt",
      # kmer solo
      "K2\nnt",
      "KUniq\nmicrobialDB" ),
    
    family_method = c(
      "Assembly-based","Assembly-based",
      "Hybrid", "Hybrid",
      "Hybrid", "Hybrid",
      "Combined assembly-free","Combined assembly-free",
      "Assembly-free","Assembly-free" ), 
    
    
    color=c("#053158","#084269", 
            "#6A0080","#4B006E",
            "saddlebrown","orangered4",
            "#275c09","#1F4A07",
            "#940C00","#730F0F" ) , 
    


    # fill = c(
    #   "#0F4F8A", "#6FA9C9",
    #   "#9E84BE", "#4A1F78",
    #   "#D98E2F", "#C45500",
    #   "#7FB85C", "#1F6F1A",
    #   "#D95C5A", "#9E0E12" ),
    # fill = c(
    #   "#9E84BE", "#F19737",
    #   "black", "#6d3b9c",
    #   "black", "#ff7e03",
    #   "black", "#2ea227",
    #   "black", "#e61a1b" ),
    fill = c(
      "#9E84BE", "goldenrod1",
      "black", "#5F2E8F",
      "black", "#E66F00",
      "black", "#2ea227",
      "black", "#C91516"
    )
    
  ))



metrics_df=merge(metrics_df, label_df)

metrics_df$superkingdom=sapply(metrics_df$taxon_phylum, function(x){
  x=strsplit(x, "_")[[1]][1]
  return(x)
})
metrics_df$phylum=sapply(metrics_df$taxon_phylum, function(x){
  x=strsplit(x, "_")[[1]][2]
  return(x)
})

std <- function(x) sd(x)/sqrt(length(x))

summary_metrics_df <- aggregate(cbind(Precision, Recall, F1) ~ superkingdom + phylum + labels2lines + method + family_method + .id + taxlvl + rt + percentage_human,
                                metrics_df, function(x) c("mean" = mean(x), "std" = std(x)))
summary_metrics_df <- do.call(data.frame, summary_metrics_df)
names(summary_metrics_df) <- gsub("\\.", "_", names(summary_metrics_df))  # Clean column names if needed



summary_metrics_df$lowerF1=summary_metrics_df$F1_mean-summary_metrics_df$F1_std
summary_metrics_df$upperF1=summary_metrics_df$F1_mean+summary_metrics_df$F1_std

summary_metrics_df$lowerPre=summary_metrics_df$Precision_mean-summary_metrics_df$Precision_std
summary_metrics_df$upperPre=summary_metrics_df$Precision_mean+summary_metrics_df$Precision_std

summary_metrics_df$lowerRec=summary_metrics_df$Recall_mean-summary_metrics_df$Recall_std
summary_metrics_df$upperRec=summary_metrics_df$Recall_mean+summary_metrics_df$Recall_std


summary_metrics_df$phylum = ifelse(summary_metrics_df$phylum == "unclassified Viruses phylum", "Incertae sedis", summary_metrics_df$phylum)


# summary_metrics_df$labels2lines = factor(summary_metrics_df$labels2lines, 
#                                          levels=c( "- Contig-based -\nnt",   
#                                                    "- K2 -\nnt",
#                                                    "- KUniq -\nmicrobialDB",
#                                                    "- Hybrid-K2 -\nnt",   
#                                                    "- Hybrid-KUniq -\nnt+microbialDB",
#                                                    "- K2-KUniq -\nnt+microbialDB"  ,
#                                                    "- KUniq-K2 -\nmicrobialDB+nt"))

### Jeter un oeil phylum par phylum ?

### Supplementary figure: All phylum
for(sk in c("Bacteria","Viruses", "Eukaryota")){}


theme_fig2 <- function(plot) {
  plot + theme(
    legend.direction = "horizontal",
    legend.position = "bottom",
    legend.key.size = unit(0.3, "cm"),
    legend.key.width = unit(0.3, "cm"),
    
    strip.background = element_rect(fill = "#333333"),
    strip.text = element_text(size = 7),  
    
    axis.title.y = element_text(size = 8, face = "bold"),
    axis.title.x = element_text(size = 8, face = "bold"),
    
    axis.text.x = element_text(size = 7),   # ~9 pt
    axis.text.y = element_text(size = 7)    # ~9 pt
  )
}
#______________________________________________________________________________#
###### Big general overview of boxplots

fill=label_df$fill
names(fill)=label_df$labels2lines

########### change size:

plots <- lapply(c("Bacteria", "Viruses","Eukaryota"), function(sk) {
  add_xlabel = FALSE
  nrow = 1
  if (sk == "Eukaryota") {
    add_xlabel = TRUE  
  } else { ylim = c(.5,1) }
  
  if (sk == "Viruses") {
    nrow=2
  }
  
  if (sk == "Bacteria") {
    sk = c("Bacteria","Archaea") 
  }
  

  
  plot_F1 = ggplot(summary_metrics_df[
    summary_metrics_df$superkingdom %in% sk & 
      summary_metrics_df$rt == "100r/t" & 
      (!str_detect(summary_metrics_df$method, "kraken2") | summary_metrics_df$method == "hybrid_kuniq-kraken2") &
      summary_metrics_df$percentage_human == "90", ],
    
    aes(x = taxlvl, y = F1_mean, color = labels2lines,
        ymin = lowerF1, 
        ymax = upperF1)) +
    geom_line(aes(group = labels2lines), linewidth = 0.65) +
    #geom_point(size = 0.4) +
    xlab("") + ylab("F1-score") +
    geom_errorbar(width = 0.1) +
    scale_color_manual(values = fill, name = "") +
    theme_linedraw() +
    theme(
      legend.position = "bottom", 
      legend.direction = "horizontal",
      legend.key.size = unit(1, "cm"),
      legend.key.width = unit(1,"cm"),
      strip.background = element_rect(fill = "#333333"),
      panel.background = element_rect(fill = "#FFFFFF"),
      legend.text = element_text(size = 9, face = "bold"),
      strip.text = element_text(size = 7, face = "bold"),
      axis.title.y = element_text(size = 7, face = "bold"),
      axis.title.x = element_text(size = 7, face = "bold"),
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size=0, color = "transparent")
    ) +
    facet_nested_wrap(vars(superkingdom, phylum), nrow = nrow)
  
  if (add_xlabel){
    plot_F1 = plot_F1+theme(
      axis.text.x = element_text(
      size = 7, 
      angle = 45, 
      hjust = 1, 
      vjust = 0.9,
      color = "black" ))
  }
  
  
  if ("Bacteria" %in% sk) {
    plot_F1 = plot_F1 + coord_cartesian(ylim = c(0.25,1))
  }
  
  return(plot_F1)
})

plots[[1]] <- plots[[1]] + theme(plot.margin = unit(c(0, .2, 0 , .2), "cm"))  # top, right, bottom, left
plots[[2]] <- plots[[2]] + theme(plot.margin = unit(c(-0.25, 0.2, -0.25, .2), "cm"))
plots[[3]] <- plots[[3]] + theme(plot.margin = unit(c(0, .2, 0, .2), "cm"))

final_plot_allphy <- wrap_plots(plots, ncol = 1, heights = c(1,2,1)) + 
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "",
    tag_suffix = ""
  ) & 
  theme(
    plot.tag = element_text(size = 11, face = "bold"),
    panel.spacing = unit(0.5, "lines"),      # smaller spacing between facets
    strip.text.y = element_text(margin = margin(0,0,0,0)),
    legend.position = "bottom",
    legend.key.size = unit(1, "cm"),
    legend.key.width = unit(1, "cm"),
  ) & 
  guides(
    color = guide_legend(override.aes = list(size = 1.5, linewidth = 1))
  )


print(final_plot_allphy)


figsupp1="~/results/Simulation/figures_review/figsupp1.jpeg"
ggsave(final_plot_allphy,
       filename = figsupp1,
       device = "jpeg",
       width = 30,
       height = 20,
       dpi = 600,
       units = "cm",
       create.dir = TRUE)


magick::image_read(figsupp1)

figsupp1svg="~/results/Simulation/figures_review/figsupp1.svg"
ggsave(final_plot_allphy,
       filename = figsupp1svg,
       device = "svg",
       width = 30,
       height = 20,
       dpi = 600,
       units = "cm",
       create.dir = TRUE)

tmp=summary_metrics_df[summary_metrics_df$taxlvl == "phylum" & summary_metrics_df$rt == "100r/t" & summary_metrics_df$percentage_human == "90",]
result <- tmp %>%
  group_by(phylum) %>%
  filter(F1_mean == max(F1_mean, na.rm = TRUE)) %>%  # handle ties
  ungroup() %>%
  count(method, name = "n_best_phylum")

result

winners <- tmp %>%
  group_by(phylum) %>%
  slice_max(order_by = F1_mean, n = 1, with_ties = TRUE) %>%
  ungroup()

winners


### Find winners at all taxonomic rank
tmp=summary_metrics_df[summary_metrics_df$rt == "100r/t" & summary_metrics_df$percentage_human == "90",]

best_by_lvl <- tmp %>%
  filter(taxlvl %in% c("phylum", "genus")) %>%
  group_by(phylum, taxlvl) %>%
  slice_max(order_by = F1_mean, with_ties = TRUE) %>%
  summarise(methods = list(method), .groups = "drop")

comparison <- best_by_lvl %>%
  pivot_wider(
    names_from = taxlvl,
    values_from = methods,
    names_prefix = "best_"
  ) %>%
  mutate(
    same_method = purrr::map2_lgl(
      best_phylum,
      best_genus,
      ~ length(intersect(.x, .y)) > 0
    )
  )

comparison


summary_count <- comparison %>%
  summarise(
    n_total = n(),
    n_same = sum(same_method, na.rm = TRUE)
  )

summary_count

##########


### Selecting phylums:

bact=c("Bacillota","Pseudomonadota","Actinomycetota","Bacteroidota", "Thermodesulfobacteriota")
vir=c("Artverviricota","Hofneiviricota","Incertae sedis","Uroviricota","Pisuviricota")

summary_metrics_df$phylum = ifelse(summary_metrics_df$phylum == "Thermodesulfobacteriota", 
                                   "Thermodesulfo-\nbacteriota", summary_metrics_df$phylum)

bact=c("Bacillota","Pseudomonadota","Actinomycetota","Bacteroidota", "Thermodesulfo-\nbacteriota")


superkingdoms=c("Bacteria","Viruses")
# Adjust sizes proportionally for smaller figure (20x18 cm)
plots_reduced <- lapply(superkingdoms, function(sk) {
  add_xlabel = FALSE
  if (sk == "Eukaryota") add_xlabel = TRUE
  if (sk == "Bacteria") {
    phylums = bact
  } else {
    phylums = vir
  }
  
  plot_F1 = ggplot(summary_metrics_df[
    summary_metrics_df$superkingdom %in% sk & 
      summary_metrics_df$rt == "100r/t" & 
      summary_metrics_df$phylum %in% phylums &
      (!str_detect(summary_metrics_df$method, "kraken2") | summary_metrics_df$method == "hybrid_kuniq-kraken2") &
      summary_metrics_df$percentage_human == "90", ],
    
    aes(x = taxlvl, y = F1_mean, color = labels2lines,
        ymin = lowerF1, ymax = upperF1)) +
    geom_line(aes(group = labels2lines), linewidth = 0.85) +  # smaller
    #geom_point(size = 0.5) +  # smaller
    xlab("") + ylab("F1-score") +
    geom_errorbar(width = 0.12) +  # smaller
    scale_color_manual(values = fill, name = "") +
    theme_linedraw() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.size = unit(0.8, "cm"),   # smaller
      legend.key.width = unit(0.8, "cm"),
      strip.background = element_rect(fill = "#333333"),
      panel.background = element_rect(fill = "#FFFFFF"),
      legend.text = element_text(size = 9, face = "bold"),
      strip.text = element_text(size = 11, face = "bold"),
      axis.title.y = element_text(size = 8, face = "bold"),
      axis.title.x = element_text(size = 8, face = "bold"),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 0.9)
    ) +
    facet_nested_wrap(vars(superkingdom, phylum), nrow = 1)
  
  if ("Bacteria" %in% sk) {
    plot_F1 = plot_F1 + coord_cartesian(ylim = c(0.25, 1))
  }
  
  return(plot_F1)
})

final_plot <- wrap_plots(plots_reduced, ncol = 1) + 
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "",
    tag_suffix = ""
  ) & 
  theme(
    plot.tag = element_text(size = 12, face = "bold"),  # smaller
    legend.position = "bottom",
    legend.key.size = unit(0.8, "cm"),
    legend.key.width = unit(0.8, "cm")
  ) & 
  guides(
    color = guide_legend(override.aes = list(size = 1, linewidth = 1))
  )

print(final_plot)

fig5jpg="~/results/Simulation/figures_review/fig5.jpg"

ggsave(final_plot,
       filename = fig5jpg,
       device = "jpeg",
       width = 29.7,
       height = 21,
       dpi = 600,
       units = "cm",
       create.dir = TRUE)

magick::image_read(fig5jpg)

fig5svg="~/results/Simulation/figures_review/fig5.svg"

ggsave(final_plot,
  filename=fig5svg,
  device = "svg",
  width = 29.7,   
  height = 21,  
  units = "cm"
)




