library(stringr)
library(plyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(tidyr)
library(patchwork) 
library(cowplot)
library(ggh4x)


sim_set="sequencing_depth"
sim_set="human_majority"


################################################################################
# F U N C T I O N S
################################################################################
ranks=c("superkingdom","phylum", "class", "order", "family", "genus", "species")


Taxonomy_small_color_palette2 = function(df, reference_colum){
  
  df=df[order(df[[reference_colum]]), ]
  df=df[!str_detect(df[[reference_colum]],"Other"), ]
  
  color_df_1 = data.frame("variable" = unique(df[str_detect(df[[reference_colum]], "Bacteria"), reference_colum  ] ) )
  color_df_1$col = colorRampPalette(c("lightblue","royalblue4" ))(nrow(color_df_1))
  
  color_df_2 = data.frame("variable" = unique(df[str_detect(df[[reference_colum]], "Eukaryota") , reference_colum ] ) )
  #color_df_2$col = colorRampPalette(c("thistle1","#960018"))(nrow(color_df_2))
  color_df_2$col = colorRampPalette(rev(c("#026e02","darkolivegreen3")))(nrow(color_df_2))
  
  color_df_3 = data.frame("variable" = unique(df[str_detect(df[[reference_colum]], "Viruses"), reference_colum ] ) )
  color_df_3$col = colorRampPalette(c("yellow2","darkorange2"))(nrow(color_df_3))
  
  
  
  color_df = rbind(color_df_1,color_df_2,color_df_3)
  
  color_dict = color_df$col
  names(color_dict) = color_df$variable
  return(color_dict) }

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
  ranks=c("superkingdom","phylum","class","order","family")
  
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
path_results=paste0(path, "results/Simulation/with_replacement/kraken/nt/Kraken2/")


sra_id=readLines(paste0(path, "data/Simulation/with_replacement/patient_name_ORhuman.txt"))
#sra_id=sra_id[str_detect(sra_id, simulation_letter)]

################################################################################
# A L G O R I T H M 
################################################################################

if (sim_set == "human_majority") {
  replicates=c(
    paste0("M",0:9),
    paste0("N",0:9),
    paste0("O",0:9),
    paste0("P",0:9),
    paste0("Q",0:9),
    paste0("R",0:9))
  
  replicate_letters=c("M","N","O","P","Q","R")
  path_reads=paste0(path,"data/Simulation/with_replacement/raw_reads_45/")
  
} else {
  
  replicates=c(
    paste0("A",0:9),
    paste0("C",0:9),
    paste0("E",0:9))
  replicate_letters=c("A","C","E")
  path_reads=paste0(path,"data/Simulation/with_replacement/raw_reads_200/")
  
}


classifs=parallel::mclapply(replicates, function(sra){
  print(sra)
  path_r=paste0(path_reads,sra,"/",sra,"_ReadNames.txt")
  print(path_r)
  r=readLines(path_r)
  tmp=as.data.frame(list("read"=r, "taxid"=sub(".*taxid:([0-9]+)_.*", "\\1", r)))
  tmp=aggregate(read ~ taxid, tmp, length)
  return(tmp)
},mc.cores=20)
names(classifs)=replicates


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

# Merge to get actual classification
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



# Alpha div global (Je suis pas sûr que les réultats soient correct)

if (sim_set == "human_majority") {
  t=list(
    "M" = Get_alpha_div(truth, classifs[1:10]),
    "N" = Get_alpha_div(truth, classifs[11:20]),
    "O" = Get_alpha_div(truth, classifs[21:30]),
    "P" = Get_alpha_div(truth, classifs[31:40]),
    "Q" = Get_alpha_div(truth, classifs[41:50]),
    "R" = Get_alpha_div(truth, classifs[51:60]) )
  letters_replicates=c("M","N","O","P",'Q',"R")[1:(length(replicates)/10)]
  
} else {
  
  t=list(
    "A" = Get_alpha_div(truth, classifs[1:10]),
    "E" = Get_alpha_div(truth, classifs[21:30]),
    "C" = Get_alpha_div(truth, classifs[11:20]) )
  letters_replicates=c("A","C","E")[1:(length(replicates)/10)]
}

alpha_div_global=ldply(t)
alpha_div <- alpha_div %>%
  separate(taxon_full_name, into = c("superkingdom", "phylum"), sep = "_", remove = FALSE)
colnames(alpha_div_global)[3]="alpha_global"


#### Retrive index
indexes=Get_diversity_indexes(truth, classifs, "family")
indexes=ldply(indexes)
colnames(indexes)[1]="replicate_id"


################################################################################

#____________________________________________ Importing results

################################################################################


### Retrieve metrics dataframes ==> All the results should be in the same folder: path_results (transfère les résultats de fulgor et centrifuger dedans) 
### Besoin d'avoir les "confusion matrix" file
path_results=paste0(path,"results/Simulation/with_replacement/")

metrics_list=list()
conf_mat_list=list()
ranks=c("superkingdom","phylum","class","order","family")
path_results=paste0(path,"results/Simulation/with_replacement/pipeline_evaluation/")
methods=list.files(path=path_results)

for(i in letters_replicates){
  if       (i=="A"){rt="5r/t"
  }else if (i %in% c("C","M","O","Q") ){rt="10r/t"
  }else if (i %in% c("E","N","P","R")){rt="100r/t"
  }
  
  
  tmp_metrics = list()
  tmp_matrix  = list()
  for (m in methods){
    
    if (i %in% c("A","C","E") & str_detect(m , "human_kuniq") ){
      next
    }
    
    metrics_file=paste0(path_results, m,"/Simulation_metrics_",m,"_",i,".RData")
    conf_mat_file=paste0(path_results,m,"/Confusion_matrix_", m,"_",i, ".tsv")
    
    
    if ( file.exists(metrics_file) ) { # if you file.resists prouves que tu existes
      if ( file.info(metrics_file)$size >= 1000 ) {
        print(m)
        load(file=metrics_file)
        metrics=get_metrics_right_format_all_ranks(metrics_full, m , rt)
        tmp_metrics[[m]]= metrics
        rm(metrics_full)
        
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

# Formating loaded results
classifs_pair=ldply(conf_mat_list)
misclassified_human=classifs_pair[classifs_pair$taxon_full_name_truth == "Eukaryota - Chordata" & classifs_pair$taxon_full_name == "Bacteria - Pseudomonadota",]
misclassified_human=misclassified_human[,c(1,10,5)]
colnames(misclassified_human)=c(".id",'method',"Human_as_proteobact")


save.image(file=paste0(path, 'results/Simulation/with_replacement/environnement_',sim_set,".RData"))


