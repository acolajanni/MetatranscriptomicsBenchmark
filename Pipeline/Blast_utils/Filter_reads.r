#!/shared/ifbstor1/software/miniconda/envs/r-4.2.3/bin/Rscript

library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
#library(ggalluvial)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

################################################################################
# F U N C T I O N S
################################################################################

get_exaequo_rank = function(sorted_dataframe, col_to_rank){
  # We use a sorted dataframe, so, the first line is the rank first
  rank_vector = c(1)
  # we then start at the second line, and get the current and previous row
  for (i in 2:nrow(sorted_dataframe)) { 
    previous = sorted_dataframe[i-1,col_to_rank]
    current = sorted_dataframe[i,col_to_rank]
    # If different, not ex aequo, so different rank
    if (! identical(as.vector(previous) , as.vector(current))  ){ rank = i }
    # If identical, same rank as previous
    else{ rank = rank_vector[i-1] } 
    rank_vector = append(rank_vector, rank) }
  sorted_dataframe$rank = rank_vector
  return(sorted_dataframe) }

get_best_taxon = function(df, col){
  agg = aggregate(list(df$rank), by=list(df[[col]]), FUN=mean)
  best_taxon = agg[agg[,2] == max(agg[,2]),1]
  return(best_taxon) }

get_frequent_taxon = function(df, col = 'taxon'){
  count=as.data.frame(table(df[[col]]))
  frequent_taxon=as.vector(count[count$Freq == max(count$Freq), ]$Var1)
  return(frequent_taxon) }

LCA_detection = function(classif_df, tax_levels=c("subject_species","strain","species",
                                                  "genus","family","order",
                                                  "class","phylum","superkingdom")){
  # if one row and species different than uncultured bacteria or synthetic construct
  # return the classification as is
  # else, find the most precise taxonomic level 
  # ex: uncultered pseudomonas sp.: most precise is phylum pseudomonas
  if (nrow(classif_df) == 1 ) { 
  
    if(str_detect(classif_df$species, "uncultured | synthetic")){
      count=4 #=genus column index
      LCA = FALSE
      while(LCA == FALSE){
        
        clade = classif_df[[ tax_levels[count] ]]
        if(clade != 0){ LCA = TRUE}
        else{count=count+1}
        
        if(count == length(tax_levels)){LCA=TRUE}
      }
      tmp = classif_df
      tmp[, tax_levels[1:count-1] ] = 0  
      return(tmp) 
    }
    else{return(classif_df)}
  } 
    
  # While there is more than one species, genus, phylum... that is not an uncultured bacteria
  # the loop continues until we find the most precise common ancestry of all queries 
  len_unique_tax = nrow(classif_df)
  count=0
  while(len_unique_tax > 1){
    count = count + 1
    unique_tax = unique(classif_df[[ tax_levels[count] ]])
    len_unique_tax = length( unique_tax )
    if(length(unique_tax[str_detect(unique_tax, 'uncultured | synthetic')]) != 0  ){
      len_unique_tax = 2
    }
  }
  
  tmp = classif_df
  tmp[, tax_levels[1:count-1] ] = 0
  return(tmp[nrow(tmp),]) }


################################################################################
# V A R I A B L E S
################################################################################
#path = "/run/user/1001/gvfs/sftp:host=core.cluster.france-bioinformatique.fr,user=acolajanni/shared/projects/microbiome_translocation"
path = "/shared/projects/microbiome_translocation/"
# path = "/home/acolajanni/Documents/work/"

setwd(path)

SRA_id = args[1]
result_dir = args[2]


#path = "/home/acolajanni/Documents/work/"
#result_dir = file.path(path, "results/Douek_Montreal/DNA/Contigs/")
# SRA_id="SRR14419357"
#setwd(path)
#SRA_in_manuscript = readLines(file.path(data_dir, "sra_subset_in_manuscript.txt"))

time=Sys.time()
print(c("Load packages - ", as.character(time)))


# query_result = read.csv(paste0(result_dir,SRA_id,"/",SRA_id,"_Blast_query.formated.csv"), row.names = 1, )

query_result = read.table(paste0(result_dir,SRA_id,"/",SRA_id,"_Blast_query_classif.formated.tsv"),
                          sep = "\t",fill = TRUE, 
                          col.names = c("contig_ID","evalue","pident","qcov_hsp","align_len","staxids","classif") )

query_result[c("superkingdom","kingdom","phylum",
               "class","order","family","genus",
               "species","strain")] = str_split_fixed(query_result$classif, fixed('|'), 9)



time=Sys.time()
print(c("Blast query loaded - ", as.character(time)))
################################################################################
# A L G O R I T H M 
################################################################################
# Remove uncomplete rows
query_result = query_result[!(is.na(as.numeric(query_result$evalue))),]
query_result$evalue = as.numeric(query_result$evalue)

# Select top 20 hit
trinity_list = list()

# time=Sys.time()
# for (contig_id in unique(query_result$contig_ID)){
#   tmp = query_result[query_result$contig_ID == contig_id , ]
#   tmp = unique(tmp)
#   tmp$log_evalue = -log10(tmp$evalue)
#   tmp = arrange(tmp, qcov_hsp,pident,log_evalue )
#   if (nrow(tmp) > 20) {trinity_list[[contig_id]] = tmp[ 1:20,  ] }
#   else{trinity_list[[contig_id]] = tmp} }
# Sys.time()-time


doParallel::registerDoParallel(cores = 12)

trinity_list = llply(.data=unique(query_result$contig_ID),
                           .fun=function(contig_id){
                              tmp = query_result[query_result$contig_ID == contig_id , ]
                              tmp = unique(tmp)
                              tmp$log_evalue = -log10(tmp$evalue + 1e-1000) 
                              tmp = plyr::arrange(tmp, qcov_hsp,pident,log_evalue )
                              if (nrow(tmp) > 20) {trinity_list[[contig_id]] = tmp[ 1:20,  ] }
                              else{trinity_list[[contig_id]] = tmp} },
                           .parallel = TRUE, .paropts = list(.export=ls())
                           )

names(trinity_list) = unique(query_result$contig_ID)

query_result = ldply(trinity_list)
query_result$.id = NULL

time=Sys.time()
print(c("Algo select 20 first hit done - ", as.character(time)))
################################################################################
# Select the taxon among best ranked, / most frequent
query_result$taxon = paste0(query_result$superkingdom,"_",query_result$phylum,"_", query_result$genus )

A_list = list() # If contig has one match only
B_list = list() # multiple taxon with one clear best among most frequent
C_list = list() # multiple taxon return the one with the best mean rank
D_list = list() # other
for (contig_id in unique(query_result$contig_ID)){
  tmp = query_result[query_result$contig_ID == contig_id , ]
  # count=as.data.frame(table(tmp$taxon))
  # frequent_taxon=as.vector(count[count$Freq == max(count$Freq), ]$Var1)
  frequent_taxon = get_frequent_taxon(tmp)
  # 1. Rank the value based on evalue, coverage, pident with ex aequo
  tmp = arrange(tmp, qcov_hsp,pident,log_evalue )
  if (nrow(tmp) > 1) { tmp = get_exaequo_rank(tmp, c("qcov_hsp","log_evalue","pident")) }
  else { 
    tmp$rank = 1
    A_list[[contig_id]] = tmp
    next
  }
  # 2. Aggregate by their mean rank - 3. select the best ranked taxon
  best_taxon = get_best_taxon(tmp, 'taxon')
  
  # First select the most frequent hit 
  # If there is equally represented taxon, use the mean ranked to select the best
  # Else, just select the best possible hit among the same taxa
  if(length(best_taxon) > 1 & length(frequent_taxon)>1){
    C_list[[contig_id]] = tmp[tmp$taxon %in% frequent_taxon & tmp$taxon %in% best_taxon,]
    
  }else if (length(best_taxon) > 1){
    
    tmp = tmp[tmp$taxon %in% best_taxon , ]
    frequent_taxon = get_frequent_taxon(tmp)
    B_list[[contig_id]] = tmp[tmp$taxon %in% frequent_taxon,]
    
  }else{ D_list[[contig_id]] = tmp[tmp$taxon == best_taxon,] }
}



time=Sys.time()
print(c("Algo select best hit done - ", as.character(time)))
################################################################################
# Filter consensual reads

filtered_trinity = c(A_list, B_list, C_list, D_list)
superkingdom_contig = lapply(filtered_trinity, function(x){
  classification = unique(x$superkingdom)
  return(classification)
})

kingdom_contig = lapply(filtered_trinity, function(x){
  classification = unique(x$kingdom)
  return(classification)
})

class_contig = lapply(filtered_trinity, function(x){
  classification = unique(x$class)
  return(classification)
})

phylum_contig = lapply(filtered_trinity, function(x){
  classification = unique(x$phylum)
  return(classification)
})

t = as.vector(superkingdom_contig)
bact = names(t[t %in% c("Bacteria","Viruses") ])

t = as.vector(kingdom_contig)
fungi = names(t[t %in% c("Fungi") ])

t = as.vector(phylum_contig)
parasites = names(t[t %in% c("Annelida","Platyhelminthes","Nematoda","Apicomplexa") ])


t_alt = as.vector(class_contig)
no_mammal = names(t_alt[! t_alt %in% c("Mammalia") ])

no_mammal_contig = filtered_trinity[no_mammal]

bacteria_contig = filtered_trinity[unique(c(bact,fungi,parasites))]
other_contig = filtered_trinity[! names(filtered_trinity) %in% names(bacteria_contig)]


## Détection de Last Common Ancestry
bact_contig_classification = ldply( lapply(bacteria_contig, function(x) LCA_detection(x)) )

bact_contig_classification$kingdom = ifelse(bact_contig_classification$kingdom == "0" & 
                                              bact_contig_classification$superkingdom == "Bacteria", 
                                            yes= bact_contig_classification$superkingdom,
                                            no = bact_contig_classification$kingdom)




print(paste0(result_dir,SRA_id,"/Contig_classification.csv"))
write.csv(bact_contig_classification, file = paste0(result_dir,SRA_id,"/Contig_classification.csv"), row.names = FALSE )


time=Sys.time()
print(c("Bacterial contig written: done - ", as.character(time) ))



# no_mammal_contig_classification = ldply( lapply(no_mammal_contig, function(x) LCA_detection(x)) )
# 
# no_mammal_contig_classification$kingdom = ifelse(no_mammal_contig_classification$kingdom == "0" & no_mammal_contig_classification$superkingdom == "Bacteria", 
#                                             yes= no_mammal_contig_classification$superkingdom,
#                                             no = no_mammal_contig_classification$kingdom)
# 
# to_see = no_mammal_contig_classification[ ! no_mammal_contig_classification$contig_ID %in% bact_contig_classification$contig_ID,]

# what's left: eukaryotic synthetic construct / plantae / arthropode

################################################################################
# regroup contig by whether or not it can be classified
# As in the VRC code

bact_contig_classification=read.csv2(file = paste0(result_dir,SRA_id,"/Contig_classification.csv"),sep=',')

tax_levels=c("strain","species", "genus","family","order",
             "class","phylum","superkingdom")




for (lvl in rev(tax_levels)){
  print(lvl)
  classif_tmp = bact_contig_classification[bact_contig_classification[[lvl]] != 0 , ]
  
  # Uncultered bacteria have not a 0 in species column 
  if (lvl == "species"){
    classif_tmp = classif_tmp[classif_tmp$family != 0, ] }
  
  path=paste0(result_dir,SRA_id,"/")
  #dir.create(paste0(path,"Quantification/"),)
  path=paste0(path,"Quantification/")

  writeLines(unique(classif_tmp$contig_ID), paste0(path,"contig_ID_",lvl,".txt") )

}


time=Sys.time()
print(c("R script: done - ", as.character(time)))



stop()

################################################################################
# TESTING
#
# t=filtered_trinity[["TRINITY_DN358_c0_g1_i1"]]
# len_unique_tax = nrow(t)
# count=0
# while(len_unique_tax > 1){
#   count = count + 1
#   print(tax_level[count])
#   len_unique_tax = length(unique(t[[ tax_level[count] ]]) )
# }
# t2 = t
# t2[, tax_level[1:count-1] ] = 0
# t2 = t2[1,]


################################################################################
# Bar plot 

tax_levels=rev(c("species","genus","family","order","class","phylum","kingdom","superkingdom"))
tax_levels = c("kingdom", "class")

count_list = list()
for (lvl in tax_levels){
  tmp = bact_contig_classification[bact_contig_classification[,lvl] != "0" , ]
  #if (! lvl %in% c("kingdom","superkingdom") ){ tmp$taxon = paste0(tmp$kingdom,"_",tmp[,lvl]) }
  #else{tmp$taxon = tmp[,lvl]}
  tmp$taxon = tmp[,lvl]
  count_df=as.data.frame(table(tmp$taxon))
  colnames(count_df) = c("taxon","freq")
  count_df$freq_percent = count_df$freq / sum(count_df$freq)
  
  count_df$lvl = lvl
  count_df = arrange(count_df, freq)
  count_list[[lvl]] = count_df
  
}

test = ldply(count_list)
test$lvl = as.factor(test$lvl)
# tmp = count_list$class
# tmp$lvl = "class"

ggplot(data = test, aes(x=lvl, y=freq_percent, fill=taxon))+
  geom_col(position = position_stack(vjust = 0.5),  color="darkgrey", size=0.05) + 
  xlab("Number of contigs mapping to a taxon") + ylab("Individual")+ #theme(legend.position = "none" )
  geom_text(aes(label = paste0(freq, " - ", taxon)),
            position = position_stack(vjust = 0.5)) +
  ggtitle("Number of contig mapping to various taxa") + labs(fill="Taxa" ) + 
  theme_linedraw() + theme(plot.title = element_text(hjust=0.5),text = element_text(size = 16))    +
  theme(legend.position = "none") #+ facet_wrap(~ .id)


tmp=reorder(tmp$taxon, tmp$freq)

################################################################################
# Alluvial Plot
tmp = bact_contig_classification[bact_contig_classification$phylum != "0" ,]
tmp_test = tmp[,c("superkingdom","kingdom","phylum","class","family")]
tmp_test = as.data.frame(table(tmp_test))
tmp_test = tmp_test[tmp_test$Freq != 0,]


tmp_test=arrange(tmp_test, kingdom,phylum)
for (lvl in colnames(tmp_test)){
  tmp_test[,lvl] = factor(tmp_test[,lvl], levels = unique(tmp_test[,lvl])) 
}


ggplot(tmp_test,
       aes(y = Freq, axis2 = phylum, axis3 = class, axis1= kingdom, axis4=family)) +
  #geom_alluvium(aes(fill = superkingdom),
  #              aes.bind = TRUE) +
  geom_flow(aes(fill = superkingdom),alpha=0.5,color="grey")+
  geom_stratum(aes(fill = superkingdom),size=.05,fill='grey') +  
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("kingdom", "Phylum","Class","Order"), expand = c(.05, .05)) +
  #scale_fill_brewer(type = "qual", palette = "Set1") +
  theme_linedraw()+
  ggtitle(" ")

################################################################################
################################################################################
# Alluvial Plot (doesn't work)
# a = melt(tmp_test)
# 
# test = test[,!colnames(test) %in% c(".id","freq")]
# 
# 
# test$freq_percent = round(test$freq_percent*100)
# 
# 
# ggplot(test,
#        aes(x = lvl, stratum = taxon, alluvium = freq_percent,
#            fill = taxon, label = taxon)) +
#   scale_fill_brewer(type = "qual", palette = "Set2") +
#   geom_flow(stat = "alluvium", lode.guidance = "frontback",
#             color = "darkgray") +
#   geom_stratum() +
#   theme(legend.position = "bottom") +
#   ggtitle("student curricula across several semesters")
# 
# is_alluvia_form( test )
# 
# test_2 = is_lodes_form(test, 
#                          key="lvl", 
#                          value="taxon", 
#                          id = "freq_percent")



