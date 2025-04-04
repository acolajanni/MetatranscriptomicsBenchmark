library(stringr)
library(plyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(caret)
library(ggtext)

args = commandArgs(trailingOnly=TRUE)

################################################################################
# F U N C T I O N S
################################################################################
ranks=c("superkingdom","phylum", "class", "order", "family", "genus", "species")

extract_tax_levels <- function(tax_string, levels = c("superkingdom", "phylum")) {
  # Split the taxonomy string by underscores
  tax_levels <- unlist(strsplit(tax_string, "_"))
  
  # Define the standard taxonomic ranks
  standard_ranks <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
  
  # Get the indices of the requested levels
  start_index <- match(levels[1], standard_ranks)
  end_index <- match(levels[2], standard_ranks)
  
  # Extract the relevant levels if indices are valid
  if (!is.na(start_index) && !is.na(end_index) && start_index <= end_index && end_index <= length(tax_levels)) {
    return(paste(tax_levels[start_index:end_index], collapse = "_"))
  } else {
    return(NA)  # Return NA if invalid levels are provided
  }
}

per_phylum_classif=function(compare_df){
  # Counting per phylum, association ==> TP/FN/FP
  compare_df$Freq=1
  counting_pred=aggregate(Freq ~ tax_lvl + taxon + phylum_pred + taxon_truth + phylum_truth , compare_df , sum)
  
  # returns a list: 
  phylum_metrics=list()
  for (p in unique(counting_pred$phylum_truth) ){
    tmp=counting_pred[counting_pred$phylum_pred == p | counting_pred$phylum_truth == p , ]
    
    tmp$classification=ifelse( (tmp$phylum_truth == tmp$phylum_pred &  tmp$taxon_truth == tmp$taxon), 
                               "TP","ELSE")
    
    tmp$classification=ifelse( (tmp$phylum_truth == p & tmp$taxon_truth != tmp$taxon), 
                               "FN",tmp$classification)
    
    tmp$classification=ifelse( (tmp$phylum_truth != p), 
                               "FP",tmp$classification)
    phylum_metrics[[p]] = tmp
  }
  
  return(phylum_metrics)
}

# Function to compute metrics for each tax_lvl and .id
compute_metrics <- function(df) {
  df %>%
    group_by(.id, tax_lvl) %>%
    summarise(
      TP = sum(Freq[classification == "TP"], na.rm = TRUE),
      FP = sum(Freq[classification == "FP"], na.rm = TRUE),
      FN = sum(Freq[classification == "FN"], na.rm = TRUE),
      Precision = ifelse((TP + FP) > 0, TP / (TP + FP), 0),
      Recall = ifelse((TP + FN) > 0, TP / (TP + FN), 0),
      F1 = ifelse((Precision + Recall) > 0, 2 * (Precision * Recall) / (Precision + Recall), 0)
    ) %>%
    ungroup()
}

metrics_per_phylum=function(compare_df){
  
  ### Get phylum truth and predicted
  phylum_truth=ldply(unique(lapply(unique(compare_df$taxon_full_name_truth),function(x){
    phylum=sub("^(([^_]+)_([^_]+)).*", "\\1", x)
    df=data.frame(list(
      "taxon_full_name_truth"= x, 
      "phylum_truth" = phylum ))
    return(df) })))
  
  phylum_pred=ldply(unique(lapply(unique(compare_df$taxon_full_name),function(x){
    phylum=sub("^(([^_]+)_([^_]+)).*", "\\1", x)
    
    df=data.frame(list(
      "taxon_full_name"= x, 
      "phylum_pred" = phylum ))
    return(df) })))
  
  compare_df=merge(compare_df, phylum_truth, by="taxon_full_name_truth")
  compare_df=merge(compare_df, phylum_pred, by="taxon_full_name")
  
  
  phylum_metrics=list()
  for (lvl in c("genus","family","order","class","phylum")){
    
    print(lvl)
    
    compare_df$taxon=sapply(compare_df$taxon_full_name, function(x) {
      return(extract_tax_levels(x, c("superkingdom",lvl) )) })
    
    compare_df$taxon_truth=sapply(compare_df$taxon_full_name_truth, function(x) {
      return(extract_tax_levels(x, c("superkingdom",lvl) )) })
    
    compare_df$tax_lvl=lvl
    
    
    phylum_metrics[[lvl]] = ldply(per_phylum_classif(compare_df))
    
  }
  
  phylum_metrics=do.call(rbind, phylum_metrics)
  
  test=aggregate(Freq ~ classification + tax_lvl + .id, phylum_metrics, sum)
  
  # Define all possible classifications
  all_classifications <- c("FN", "TP", "FP")
  
  # Create all possible combinations
  df_complete <- test %>%
    distinct(.id, tax_lvl ) %>%
    tidyr::expand(.id, tax_lvl, classification = all_classifications) %>%
    left_join(test, by = c(".id", "tax_lvl", "classification")) %>%
    mutate(Freq = tidyr::replace_na(Freq, 0)) # Assign 0 for missing frequencies
  
  # Apply function to your dataset
  df_metrics <- compute_metrics(df_complete)
  
  colnames(df_metrics)[1]="taxon_phylum"
  
  return(df_metrics)
}

std <- function(x) sd(x)/sqrt(length(x))


get_classif_at_lvl=function(classif_df, lvl){
  if (lvl == "superkingdom"){
    wanted_ranks=lvl
    classif_df$taxon_full_name = classif_df$superkingdom
  } else{
    wanted_ranks=ranks[1:which(ranks == lvl)]
    classif_df$taxon_full_name <- apply(classif_df[, wanted_ranks], 1, paste, collapse = "_")
  }
  classif_df$taxon = classif_df[[lvl]]
  classif_df[,ranks]=NULL
  classif_df$tax_lvl = lvl
  return(classif_df)
}

### Function: 
# truth_df columns: taxid + all ranks
# classif_df: readname + true taxid + all ranks
get_classif_trueVSexp=function(truth_df, classif_df, tax_lvl){
  
  classif_lvl= get_classif_at_lvl(classif_df, tax_lvl)
  truth_lvl  = get_classif_at_lvl(truth_df, tax_lvl)
  
  colnames(truth_lvl)[colnames(truth_lvl)=="taxon_full_name"] ="taxon_full_name_truth"
  colnames(truth_lvl)[colnames(truth_lvl)=="taxon"] ="taxon_truth"
  truth_lvl$tax_lvl=NULL
  truth_lvl=unique(truth_lvl)
  
  classif_compare=merge(classif_lvl, truth_lvl, by="taxid")
  return(classif_compare)
}

### Function: take classification prediction vs truth - For one taxonomic level
### input:
# classif_df: 1 line = 1 read ---  taxon_full_name + taxon_full_name_truth 
get_metrics_from_classif=function(classif_df){
  
  # Using caret's confusionMatrix function
  all_levels <- union(levels(factor(classif_df$taxon_full_name_truth)), levels(factor(classif_df$taxon_full_name)))
  conf_matrix <- confusionMatrix(factor(classif_df$taxon_full_name, levels = all_levels), factor(classif_df$taxon_full_name_truth, levels = all_levels) )
  metrics=as.data.frame(conf_matrix$byClass)
  metrics=metrics[ !is.na(metrics$F1), c("Precision","Recall","F1","Balanced Accuracy","Sensitivity","Specificity")]
  metrics$taxon= str_remove(row.names(metrics), "Class: ")
  return(metrics)
}


### 1 dataframe with taxon_full_name(_truth) + taxon(_truth)
get_metrics_from_classif_df=function(classif_df, ranks){
  
  classif_metrics=lapply(ranks, function(lvl){
    print(lvl)
    
    # If dataframe already in good format, no need to change it
    if ( unique(classif_df$tax_lvl) != lvl ){ 
      df=extract_taxon_dataframe(classif_df, lvl)
      
    } else{ 
      df = classif_df } 
    
    return( get_metrics_from_classif(df) )  } )
  
  
  names(classif_metrics) = ranks
  return(classif_metrics)
  
}


### Function: take classification prediction vs truth - For all taxonomic level
get_metrics_from_classif_list=function(classif_list, ranks){
  
  classif_metrics=lapply(ranks, function(lvl){
    print(lvl)
    return( get_metrics_from_classif(classif_list[[lvl]])  ) })
  names(classif_metrics) = ranks
  return(classif_metrics)
}


extract_taxon_dataframe <- function(df, rank) {
  ranks <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
  idx <- match(rank, ranks)
  
  
  # Concatenate taxon_full_name and taxon_full_name_truth for unique processing
  combined <- unique(c(df$taxon_full_name, df$taxon_full_name_truth))
  result <- data.frame(combined_key = combined, taxon = NA, stringsAsFactors = FALSE)
  
  for (i in seq_along(combined)) {
    words <- unlist(strsplit(combined[i], "_"))
    
    result$taxon_full_name[i] <- paste(words[1:idx], collapse = "_")
    result$taxon[i] <- words[idx]
  }
  
  # Rename for merging later
  colnames(result) <- c("combined_key","taxon","taxon_full_name" )
  result_truth <- setNames(result, c("combined_key", "taxon_truth","taxon_full_name_truth"))
  
  
  # Merge results back to original dataframe based on taxon_full_name then on taxon_full_name_truth 
  colnames(df)[colnames(df) == "taxon_full_name"] = "combined_key"
  df <- merge(df[, ! colnames(df) %in% c("taxon", "taxon_full_name")], result, by = "combined_key", all.x = TRUE)
  df$combined_key=NULL
  
  colnames(df)[colnames(df) == "taxon_full_name_truth"] = "combined_key"
  df <- merge(df[, ! colnames(df) %in% c("taxon_truth", "taxon_full_name_truth")], result_truth, by = "combined_key", all.x = TRUE)
  df$combined_key=NULL
  
  df$tax_lvl=rank
  
  return(df)
}




Get_alpha_div_per_samples=function(truth_df, classif_df, return_table=FALSE){
  classif_df$taxid <- sub(".*taxid:([0-9]+)_.*", "\\1", classif_df$read)
  classif_df=classif_df[,c("read","taxid")]
  # Counting number of read per species
  classif_df=aggregate(read ~ taxid, classif_df, length)
  
  if(return_table==TRUE){return(classif_df)}
  
  # Merging with real annotation
  alphadiv = merge(classif_df, truth_df[c("taxid", "superkingdom","phylum","family","species")], by="taxid")
  # Alpha div per phylum
  alphadiv = aggregate(species ~ superkingdom + phylum , alphadiv ,length)
  return(alphadiv)
}


Get_alpha_div=function(truth_df, classif_list){
  # Get occurrences of each taxid
  classif_list=lapply(classifs, function(x){ Get_alpha_div_per_samples(truth_df, x, TRUE) })
  classif_df=do.call(rbind, classif_list)
  classif_df=aggregate(read ~ taxid, classif_df, sum)
  
  # Merging with real annotation
  alphadiv = merge(classif_df, truth_df[c("taxid", "superkingdom","phylum","family","species")], by="taxid")
  # Alpha div per phylum
  alphadiv = aggregate(species ~ superkingdom + phylum , alphadiv ,length)
  
}



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

################################################################################
# V A R I A B L E S
################################################################################

# path="/home/acolajanni/Documents/work/"
# path="/run/user/1001/gvfs/sftp:host=core.cluster.france-bioinformatique.fr,user=acolajanni/shared/projects/microbiome_translocation/"
path="/shared/projects/microbiome_translocation/"


#setwd("/home/acolajanni/Documents/work/")

simulation_letter <- ifelse(length(args) >= 1, args[1], "M")
result_dir <- ifelse(length(args) >= 2, args[2], paste0(path, "results/Simulation/with_replacement/"))
strategy <- ifelse(length(args) >= 3, args[3], "normal") # normal / hybrid

method1 <- ifelse(length(args) >= 4, args[4], "Centrifuger") # c("Blast","kraken2","kuniq")
method2 <- ifelse(length(args) >= 5, args[5], "none") # c("Blast","kraken2","kuniq","none")
unmapped <- ifelse(length(args) >= 6, args[6], "aligned") # c("unclassified_kuniq","unclassified_kraken2","human_kuniq","aligned")
compare_assembled_only <- ifelse(length(args) >= 7, args[7], "FALSE") # FALSE

print(simulation_letter) 
print(strategy) 
print(method1) 
print(method2) 
print(unmapped) 
print(compare_assembled_only) 

print('--------')
print(result_dir)
print('--------')


sra_id=paste0( simulation_letter,c(0:9) )

kuniq_path   = paste0(result_dir,"/kraken/microbialDB/kuniq/")
kraken2_path = paste0(result_dir,"/kraken/nt/k2uniq/")
centrifuger_path = paste0(result_dir,"/kraken/RefSeq_noEukaryota/Centrifuger/")

################################################################################
# A L G O R I T H M 
################################################################################

strategy_name=method1

if(method2  != "none"){    strategy_name=paste0(strategy_name,"-",method2)  }
if(strategy != "normal"){  strategy_name=paste0("hybrid_",strategy_name)  }
if(unmapped != "aligned"){ strategy_name=paste0(unmapped,"_",strategy_name)  }

print(strategy_name)

dir.create( paste0(result_dir,"/pipeline_evaluation/") , showWarnings = FALSE)
dir.create( paste0(result_dir,"/pipeline_evaluation/",strategy_name) , showWarnings = FALSE)

output_dir=paste0(result_dir,"/pipeline_evaluation/",strategy_name,"/")




if(method2=="Blast"){ order="reversed" 
}   else{             order="forward"}


if (method1=="kuniq"){ 
  path_classif1=kuniq_path
  data1="kraken"
  
}else if (method1=="kraken2"){
  path_classif1=kraken2_path
  data1="kraken"
  
} else if (method1=="Blast"){
  data1="salmon"
  if (unmapped == "aligned"){ path_classif1=paste0(result_dir,"/Contigs/") 
  } else if(unmapped == "human_kuniq"){ path_classif1=paste0(result_dir,"/hybrid/human/unmapped_kuniq/Contigs/") }

}else if (method1=="Centrifuger"){
  path_classif1=centrifuger_path
  data1="kraken"
}


if (method2=="kuniq"){ 
  path_classif2=kuniq_path
  data2="kraken"
}else if (method2=="kraken2"){
  path_classif2=kraken2_path
  data2="kraken"
} else if (method2=="Blast"){
  data2="salmon"
  if(unmapped == "unclassified_kuniq"){          path_classif2=paste0(result_dir,"/hybrid/unclassified/unmapped_kuniq/") 
  } else if(unmapped == "unclassified_kraken2"){ path_classif2=paste0(result_dir,"/hybrid/unclassified/unmapped_kraken2/") }
}

### Load classification - First one
print("loading files")
start_time <- proc.time()

classifs_1=parallel::mclapply(sra_id, function(sra){load_ReadsClassif(sra, path_classif1, data1)},mc.cores=10)
names(classifs_1)=sra_id


### Load KrakenUniq if we want to filter human first with krakenuniq
if(unmapped == "human_kuniq") {
  
  print("use krakenuniq to filter human reads")
  
  # load
  path_kraken=paste0(result_dir,"/kraken/microbialDB/kuniq/")
  human_classifs=parallel::mclapply(sra_id, function(sra){load_ReadsClassif(sra, path_kraken, "kraken")} ,mc.cores=10)
  names(human_classifs)=sra_id
  
  # Removing human classified reads with KrakenUniq ; other classification by the first classifier of the list
  classifs_1=parallel::mclapply(sra_id, function(sra){
    ku =human_classifs[[sra]]
    A  =classifs_1[[sra]]
    # human reads: krakenuniq classification as mammals
    human_reads=ku[ku$class == "Mammalia",]
    # other reads: keep the classification for the other reads
    other_reads=A[! A$read %in% human_reads$read ,]
    hybrid=rbind(human_reads,other_reads)
    return(hybrid)
  },mc.cores=10)
  names(classifs_1)=sra_id
}


if(strategy == "normal"){
  
  # Straightforward, nothing more is needed 
  classifs=classifs_1
  rm(classifs_1)
  
}else{
  
  classifs_2=parallel::mclapply(sra_id, function(sra){load_ReadsClassif(sra, path_classif2, data2)},mc.cores=10)
  names(classifs_2)=sra_id
  
  ###
  # à vérifier : Marcherait avec reversed ?
  ###
  # Removing unassembled and unclassified read from first tool, taking classification from the second
  classifs=parallel::mclapply(sra_id, function(sra){
    A=classifs_1[[sra]]
    B=classifs_2[[sra]]
    
    # Identify unclassified/unassembled reads
    misclassified=A[ A$superkingdom %in% c("unassembled", "unclassified"), ]$read
    classified_A=A[! A$superkingdom %in% c("unassembled", "unclassified"),] 
    
    # Take the classification of these reads with the other classifier
    classified_B=B[B$read %in% misclassified,] 
    hybrid=rbind(classified_A,classified_B)
    return(hybrid)
  },mc.cores=10)
  
  
}
names(classifs)=sra_id

end_time <- proc.time() - start_time
print(paste0("Files loaded, took ", round(end_time[3]), " seconds" ))

### Looking only at reads that are assembled into contigs
if (compare_assembled_only == "TRUE"){
  
  if(unmapped == "unclassified_kuniq"){          
    path_unassembled=paste0(result_dir,"/hybrid/unclassified/unmapped_kuniq/Contigs/") 
  } else if(unmapped == "unclassified_kraken2"){ 
    path_unassembled=paste0(result_dir,"/hybrid/unclassified/unmapped_kraken2/Contigs/") 
  } else if (unmapped == "human_kuniq"){
    path_unassembled=paste0(result_dir,"/hybrid/human/unmapped_kuniq/Contigs/") 
  } else{
    path_unassembled=paste0(result_dir,"/Contigs/") }
  
  
  classifs=parallel::mclapply(sra_id, function(sra){
    
    p=paste0(path_unassembled,sra,"/ContigsToReads/missing_ids.txt")
    unassembled=readLines(p) 
    c=classifs[[sra]]
    c=c[!c$read %in% unassembled,]
    return(c) },mc.cores=10)
  
  names(classifs)=sra_id
  
  output_dir=paste0(result_dir,"/pipeline_evaluation/",strategy_name,"/")
  
  if (compare_assembled_only == "TRUE"){
    
    dir.create( paste0(result_dir,"/pipeline_evaluation/assembled_reads/") , showWarnings = FALSE)
    dir.create( paste0(result_dir,"/pipeline_evaluation/assembled_reads/",strategy_name) , showWarnings = FALSE)
    output_dir= paste0(result_dir,"/pipeline_evaluation/assembled_reads/",strategy_name,"/")
    
  }
}




### Rajouter les reads manquants pour l'humain ### ( G H I J K L )
if (str_detect(sra_id[1], "I|L|H|K")){
  for (sra in sra_id){ 
    
    human_reads_paths=paste0(result_dir,"/aligned_ReadsNames/HumanReads/",sra,"_human.txt")
    x=read.table(human_reads_paths, header=FALSE, fill = TRUE, sep = "\t", quote = "")
    x[,c("superkingdom","phylum", "class", "order", "family", "genus", "species", "strain")] =  
      ldply(rep(list(c("Eukaryota", "Chordata", "Mammalia", "Primates", 
                       "Hominidae", "Homo", "Homo sapiens", "Homo sapiens")), 
                each = nrow(x)))
    x[,2]=NULL
    colnames(x)=c("read", "superkingdom","phylum", "class", "order", "family", "genus", "species", "strain")
    x$strain=NULL
    x$taxid <- sub(".*taxid:([0-9]+)_.*", "\\1", x$read)
    
    classifs[[sra]] = rbind(classifs[[sra]], x)
    
  }
}


############
# 3: Build Results vs Truth

# 3.1: retrieve dataframe of read origin

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

# 3.2.2 get dataframe in right format (the most precise that we want = family)
compare_df=parallel::mclapply(classifs, function(x){
  return(get_classif_trueVSexp(truth, x, "genus") )  },mc.cores=10)


metrics_full=parallel::mclapply(compare_df, function(tmp){
  #tmp=get_classif_trueVSexp(truth, x, "family")  
  return(get_metrics_from_classif_df(tmp, ranks[1:5] ) )
},mc.cores = 10)


names(compare_df) = sra_id

##### Metrics for each phylum
metrics_phylum=parallel::mclapply(compare_df, function(x){metrics_per_phylum(x)}, mc.cores = 10)
names(metrics_phylum) = sra_id

#metrics_phylum=ldply(metrics_phylum)




####
save(metrics_phylum, file=paste0(output_dir,"Simulation_metrics_phylum_",strategy_name,"_",simulation_letter,".RData" ) )
save(metrics_full, file=paste0(output_dir,"Simulation_metrics_",strategy_name,"_",simulation_letter,".RData" ) )
#load(paste0(output_dir,"Simulation_metrics_",strategy_name,"_",simulation_letter,".RData" ))
####

stop()

# Save current environment to restart at this step and debug: find where the NA comes
# save.image('/home/acolajanni/Documents/work/results/Simulation/debug_env.RData')
# load('/home/acolajanni/Documents/work/results/Simulation/debug_env.RData')

############
### Get the relative frequency for each true label. 
# 1 - Merge incorrect predicted taxa
# Retrieve misclassified taxa
test=ldply(parallel::mclapply(compare_df, function(x) {
  return(extract_taxon_dataframe(x, "phylum"))  },mc.cores=10))

incorrect_pred=unique(test[!test$taxon_full_name %in% test$taxon_full_name_truth , ]$taxon_full_name)
incorrect_pred=incorrect_pred[! str_detect(incorrect_pred,"unclassified") ]

incorrect_pred=as.data.frame(incorrect_pred)
incorrect_pred[,c("sk","p")] = str_split_fixed(incorrect_pred$incorrect_pred, fixed("_"),n=2 )
incorrect_pred$p = NULL
incorrect_pred$sk = ifelse(incorrect_pred$sk == "unassembled", incorrect_pred$sk, paste0("Misclassified ",incorrect_pred$sk) )


un=data.frame("incorrect_pred" = "unclassified_unclassified", "sk" = "unclassified")

incorrect_pred=rbind(incorrect_pred, un)
colnames(incorrect_pred) = c("taxon_full_name","predicted")

already_present=incorrect_pred[NULL,]
present=unique(test[!test$taxon_full_name %in% incorrect_pred$taxon_full_name , ]$taxon_full_name)
already_present=data.frame("taxon_full_name"=present, "predicted"=present)

dict=rbind(already_present,incorrect_pred)


# 1bis- get the total seq for each true value 
test=merge(test, dict, by="taxon_full_name")

conf_mat=as.data.frame(table("taxon_full_name"=test$predicted, "taxon_full_name_truth"=test$taxon_full_name_truth))

total_seq=as.data.frame(table("taxon_full_name_truth"=test$taxon_full_name_truth))
colnames(total_seq) = c("taxon_full_name_truth", "total")

conf_mat=merge(total_seq, conf_mat, by="taxon_full_name_truth")
conf_mat$prop=conf_mat$Freq / conf_mat$total

conf_mat$taxon_full_name_truth = factor(conf_mat$taxon_full_name_truth)
conf_mat$taxon_full_name = factor(conf_mat$taxon_full_name, levels = rev(c(
  unique(conf_mat$taxon_full_name[!conf_mat$taxon_full_name %in% conf_mat$taxon_full_name_truth]),
  unique(conf_mat$taxon_full_name_truth)
)))


conf_mat$color=NULL
conf_mat$label=NULL

# Wanted order
test_factor = rev(c(
  "unassembled", "unclassified" , "Archaea_unclassified", "Other Archaea",
  "Bacteria_unclassified", "Other Bacteria",  
  "Eukaryota_unclassified", "Other Eukaryota",
  "Viruses_unclassified", "Other Viruses",
  unique(as.character(conf_mat$taxon_full_name_truth))
))

test_factor = rev(c(
  "unassembled", "unclassified" , 
  "Misclassified Archaea", 
  "Misclassified Bacteria",  
  "Misclassified Eukaryota",
  "Misclassified Viruses", 
  unique(as.character(conf_mat$taxon_full_name_truth))
))

values_confmat=unique(conf_mat$taxon_full_name)
missing_id1=as.character(values_confmat[! values_confmat %in% test_factor])

test_factor=str_to_title(str_replace(test_factor, pattern = "_", replacement = " - "))

conf_mat$taxon_full_name=str_to_title(str_replace(conf_mat$taxon_full_name, pattern = "_", replacement = " - "))
conf_mat$taxon_full_name_truth=str_to_title(str_replace(conf_mat$taxon_full_name_truth, pattern = "_", replacement = " - "))



missing_ids= c(
  as.character(unique(conf_mat[ !conf_mat$taxon_full_name_truth  %in% unique(conf_mat$taxon_full_name) , ]$taxon_full_name_truth)),
  test_factor[!test_factor %in% unique(as.character(conf_mat$taxon_full_name))]
)
print(missing_ids)

true_labels=unique(conf_mat$taxon_full_name_truth)

### Need to add the missing true label to keep the diagonal in the figure \\ penser à rajouter le label '0' dans la diagonale
for (missing in missing_ids){
  print(missing)
  
  missing_df=data.frame(
    "taxon_full_name_truth"=true_labels ,
    "total"=rep(unique(conf_mat$total), length(true_labels)),
    "taxon_full_name"=rep(missing,length(true_labels)),
    "Freq"=rep(0,length(true_labels) ),
    "prop"=rep(0,length(true_labels)) )
  
  conf_mat=unique(rbind(missing_df, conf_mat))
}

#
final_factor = factor(conf_mat$taxon_full_name, levels = test_factor)
print(length(unique(final_factor)))
print(length(unique(conf_mat$taxon_full_name)))

missing=unique(conf_mat$taxon_full_name)[!unique(conf_mat$taxon_full_name) %in% levels(final_factor)]
print("missing_id:")
print(missing)
for(m in missing){
  
  conf_mat$pred = ifelse(conf_mat$taxon_full_name == m & str_detect(m, "Viruses") , 
                         "Misclassified Viruses", conf_mat$taxon_full_name )
  
  conf_mat$pred = ifelse(conf_mat$taxon_full_name == m & str_detect(m, "Archaea") , 
                         "Misclassified Archaea", conf_mat$pred )
  
  conf_mat$pred = ifelse(conf_mat$taxon_full_name == m & str_detect(m, "Eukaryota") , 
                         "Misclassified Eukaryota", conf_mat$pred )
  
  conf_mat$pred = ifelse(conf_mat$taxon_full_name == m & str_detect(m, "Bacteria") , 
                         "Misclassified Bacteria", conf_mat$pred )
  
  conf_mat$taxon_full_name = conf_mat$pred 
  conf_mat$pred = NULL
}

final_factor = factor(conf_mat$taxon_full_name, levels = test_factor)

levels(final_factor)[!levels(final_factor) %in% unique(conf_mat$taxon_full_name)]
unique(conf_mat$taxon_full_name)[!unique(conf_mat$taxon_full_name) %in% levels(final_factor)]

conf_mat$taxon_full_name=final_factor 


conf_mat <- conf_mat %>%
  group_by(taxon_full_name_truth, total, taxon_full_name) %>%
  summarise(
    Freq = sum(Freq),
    prop = sum(prop),
    .groups = 'drop'  # Avoids unnecessary grouping
  )


conf_mat$label = ifelse(conf_mat$prop != 0, conf_mat$Freq, "" )
conf_mat$label = ifelse(conf_mat$taxon_full_name_truth == conf_mat$taxon_full_name & conf_mat$Freq == 0, 
                        "0", conf_mat$label )



conf_mat$prop = conf_mat$prop * 100

conf_mat$label_freq = ifelse(conf_mat$prop != 0, round(conf_mat$prop*10, 1), "" )
conf_mat$label_freq = ifelse(conf_mat$taxon_full_name_truth == conf_mat$taxon_full_name & conf_mat$Freq == 0, 
                             "0", conf_mat$label_freq )
conf_mat$color = ifelse(conf_mat$prop >= 50, "white", "black" )


################################################################################
### Coloring taxa based on their names

# Create colors based on the taxon
colors <- ifelse(grepl("Bacteria", conf_mat$taxon_full_name), "#27408B",       # Dark Blue for Bacteria
                 ifelse(grepl("Viruses", conf_mat$taxon_full_name), "#FFB300",   # Orange for Viruses
                        ifelse(grepl("Eukaryota", conf_mat$taxon_full_name), "#FF3333",  # Lighter Red for Eukaryota
                               ifelse(grepl("Archaea", conf_mat$taxon_full_name), "#5F7DAD",  # Light Blue for Archaea
                                      "black"))))  # Default to black for others


colors=ifelse(grepl("Misclassified", conf_mat$taxon_full_name), "black", colors)

# Generate a named vector with the taxon name as names and the color as values
colored_labels <- setNames(colors, conf_mat$taxon_full_name)

# Create a custom ggplot
plot=ggplot(conf_mat, aes(x=taxon_full_name, y=taxon_full_name_truth, fill = prop)) +
  geom_tile(width = .95 , height = .95) + 
  geom_text(aes(label=label, color = color), size = 7.5/.pt) + 
  scale_color_identity() +
  scale_fill_viridis_c(limits = c(0, 100), oob = scales::squish, option="rocket", direction = -1) +
  theme_dark() + 
  ylab("True label") + xlab("Predicted label") +
  labs(fill = "(%)") +
  theme(
    axis.text.x = element_markdown(angle = 30, size = 12, vjust = 1, hjust = 1, face = "bold"),
    axis.text.y = element_markdown(size = 12, face = "bold"),
    panel.background = element_rect(fill = "lightgrey"),
    panel.grid.major.y = element_line(color = "darkgrey"),
    panel.grid.major.x = element_line(color = "black"),
    legend.key.size = unit(0.75, 'cm'),
    legend.key.height = unit(.75, 'cm'), 
    legend.key.width = unit(.75, 'cm'),
    legend.title = element_text(size=14, face="bold"), 
    legend.text = element_text(size=10)) + 
  scale_x_discrete(labels = function(x) {
    # Apply the color to the x-axis labels (taxon names) using HTML <span>
    sapply(x, function(taxon) { paste0("<span style='color:", colored_labels[taxon], "'>", taxon, "</span>")
    })
  }) + # Apply custom x-axis labels with colored text
  scale_y_discrete(labels = function(x) {
    sapply(x, function(taxon) { paste0("<span style='color:", colored_labels[taxon], "'>", taxon, "</span>") })
  })  # Apply custom x-axis labels with colored text

plot
#save.image(file = paste0(output_dir,strategy_name,"_",simulation_letter,"_simulation_letter_env.RData" ))
#save(metrics_full, file=paste0(output_dir,"Simulation_metrics_",simulation_letter,".RData" ) )

write.table(conf_mat, file = paste0(output_dir,"Confusion_matrix_",strategy_name,"_",simulation_letter,".tsv" ),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")

ggsave(plot,filename= paste0(output_dir,"Confusion_heatmap_",strategy_name,"_",simulation_letter,".png" ) , 
       device = "png",height = 34, width = 60, units = "cm")

stop()
