#!/shared/ifbstor1/software/miniconda/envs/r-4.2.3/bin/Rscript
library(plyr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

################################################################################
# V A R I A B L E S
################################################################################
path = "/run/user/1001/gvfs/sftp:host=core.cluster.france-bioinformatique.fr,user=acolajanni/shared/projects/microbiome_translocation/"
# path = "/shared/projects/microbiome_translocation/"
# path = "/home/acolajanni/Documents/work/"

letter_replicate = args[1]
# result_dir = args[2]


path_data=paste0(path,"data/Simulation/with_replacement/")
sra_id=readLines(paste0(path_data, "patient_name_ORhuman.txt"))

sra_id=sra_id[str_detect(sra_id,letter_replicate)]

setwd("/home/acolajanni/Documents/work/")


path_readNames=paste0(path_data,"aligned_ReadsNames/HumanReads/")

# AlignedReads=parallel::mclapply(sra_id, function(sra){
#   print(sra)
#   x=read.delim(paste0(path_readNames,sra,"_human.txt"), header = FALSE)
#   colnames(x) = c("readname","mapper")
#   x$taxid <- sub(".*taxid:([0-9]+)_.*", "\\1", x$readname)
#   return(x) },
#   mc.cores=10)
# names(AlignedReads)=sra_id

### Global parameters df:
total_reads = data.frame(
  ID = c("I", "L", "H", "K", "G", "J"),
  transcripts = c(7400 , 66600, 7400, 66600, 7400, 66600) , 
  human_perc = c( 0.5 , 0.9, 0.5, 0.9, 0.5, 0.9) ,
  rpert = c(100, 100, 10, 10, 5, 5) )
total_reads$human_reads= total_reads$transcripts * total_reads$rpert
total_reads$total=total_reads$human_reads/total_reads$human_perc
total_reads$non_human=total_reads$total - total_reads$human_reads



### 1: Nombre de faux positif + faux négatif (= on connait le total)
mapper=c("Bowtie2_mapping_hg38" , "STAR_mapping_hg38" ,
         "Bowtie2_mapping_chm13", "STAR_mapping_chm13" ,
         "Bowtie2_mapping_hg19" , "STAR_mapping_hg19")
complete_mappers <- data.frame(mapper = rep(mapper,2), readname = 0, taxid=c(rep("TP",6),rep("FP",6) ))

complete_mappers <- expand.grid(mapper = mapper, taxid = c("TP","FP"),
                                     stringsAsFactors = FALSE)
complete_mappers$readname=0

# count_reads=lapply(names(AlignedReads), function(x){
#   t=AlignedReads[[x]]
#   t=aggregate(readname ~ mapper + taxid, t, length)
#   t$taxid=ifelse(t$taxid == 9606, "TP", "FP")
#   t=aggregate(readname ~ mapper + taxid, t, sum)
#   
#   t <- merge(complete_mappers, t, by = c("mapper","taxid"), all.x = TRUE)
#   # Replace NA in the readname column with the value from `complete_mappers`
#   t$readname <- ifelse(is.na(t$readname.y), 
#                        t$readname.x, 
#                        t$readname.y)
#   t$mapper = factor(t$mapper, levels = mapper, ordered = TRUE)
#   
#   # Cumulative sum at each step
#   t = t %>% 
#     arrange(mapper) %>%
#     group_by(taxid) %>%
#     mutate(cumsum = cumsum(readname))
#   
#   # Get letter serie of simulation(H,I,L, ...)
#   letter = substring(x, 1, 1)
#   t$serie = letter
#   t[,c("human_reads","total","non_human")] = total_reads[total_reads$ID == letter, c("human_reads","total","non_human")]
#   
#   # Get remaining number of reads nonhuman/human after each step
#   t <- t %>%
#     mutate(remaining = if_else(
#       taxid == "TP", 
#       human_reads - cumsum,  # Subtract readname from human_reads when taxid == "TP" ==> Faux Négatifs
#       non_human - cumsum     # Subtract readname from non_human otherwise            ==> Faux 
#     ))
#   
# 
#   return(t[ , ! colnames(t) %in% c("readname.x", "readname.y")])
# })
# names(count_reads) = sra_id
# 
# # Reformat dataframes to obtain 4 rows for a mapper: TP / FP / TN / FN with cumulative values
# count_reads=lapply(names(count_reads), function(x){
#   t=count_reads[[x]]
#   
#   t_neg=t
#   t_neg$cumsum = t_neg$remaining
#   t_neg$taxid = ifelse(t_neg$taxid == "FP", "TN", "FN")
#   t_neg$remaining=NULL
#   t$remaining=NULL
#   
#   t= rbind(t, t_neg)
#   
#   t$mapper = factor(t$mapper, levels = mapper, ordered = TRUE)
#   
#   return(t)
# })
# names(count_reads) = sra_id


save(count_reads, file = paste0(path,"/results/Simulation/with_replacement/rdata/",letter_replicate,"_count_reads.RData") )
### Above for IFB
### Bellow exploratory
exit()
readcount=list()
for( i in c("I","L","H","K") ){
  
  path_rdata = "/home/acolajanni/Documents/work/results/Simulation/with_replacement/rdata/"
  
  load(paste0(path_rdata,i,"_count_reads.RData"))
  readcount=append(readcount,count_reads)
}


count_reads=ldply(readcount)


count_reads$taxid = factor(count_reads$taxid, levels = rev(c("TP","FN", "FP", "TN")))


result_df <- do.call(rbind, lapply(seq_len(nrow(total_reads)), function(i) {
  row <- total_reads[i, ]
  data.frame(
    cumsum = c(row$human_reads, row$non_human, 0, 0),
    readname = c(row$human_reads, row$non_human, 0, 0),
    taxid = c("TP", "TN", "FP", "FN"),
    serie = row$ID,
    human_reads = row$human_reads,
    total = row$total,
    non_human = row$non_human )
}))
result_df$cumsum = result_df$cumsum * 10
result_df$readname = result_df$readname * 10
result_df$mapper = "Truth"
result_df$.id = "Truth"
count_reads_merge=rbind(count_reads, result_df)

custom_colors <- c(
  "TP" = "darkred",
  "FN" = "red3",
  "TN" = "darkblue",
  "FP" = "turquoise"
)

count_reads_merge$mapper = factor(count_reads_merge$mapper, levels = rev(c("Truth",mapper)), ordered = TRUE)


count_reads_merge$legend=ifelse(count_reads_merge$serie == "H", 
                                "H: 74K Human reads per replicate - (50%)" , count_reads_merge$serie)
count_reads_merge$legend=ifelse(count_reads_merge$legend == "K", 
                                "K: 666K Human reads per replicate - (90%)" , count_reads_merge$legend)
count_reads_merge$legend=ifelse(count_reads_merge$legend == "I", 
                                "I: 740K Human reads per replicate - (50%)" , count_reads_merge$legend)
count_reads_merge$legend=ifelse(count_reads_merge$legend == "L", 
                                "L: 6660K Human reads per replicate - (90%)" , count_reads_merge$legend)

ggplot(count_reads_merge[!count_reads_merge$serie %in% c("G","J"),], aes(y=mapper, x=cumsum, fill = taxid))+
  geom_bar(position="stack", stat="identity") + 
  #ylab("# of Correctly Mapped reads against Human Genomes") + 
  ylab("Mapping steps") + xlab("Number of mapped read") +
  geom_vline(aes(xintercept = human_reads*10),count_reads_merge[!count_reads_merge$serie %in% c("G","J"),]) +
  scale_fill_manual(values = custom_colors) +
  theme_linedraw()+
  theme(
    axis.text.x = element_text(angle = 30, size = 10, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    panel.grid.major.x = element_line(color = "black"),
    strip.text = element_text(size = 12, face = "bold"),
    legend.key.size = unit(0.75, 'cm'),
    legend.key.height = unit(.75, 'cm'), 
    legend.key.width = unit(.75, 'cm'),
    legend.title = element_text(size=14, face="bold"), 
    legend.text = element_text(size=10)) + 
  facet_wrap(~ legend,nrow = 2, scales = "free_x")





### 2: calculer le nombre de reads restants à chaque étape // ==> recup seqkit stats ?
### 2bis: Récup faux négatifs et vrai négatifs

#last_step=count_reads_merge[count_reads_merge$mapper == "STAR_mapping_hg19" , ]
reshaped_df <- count_reads_merge[count_reads_merge$mapper != "Truth" , ] %>%
  pivot_wider(
    names_from = taxid,    # Use the taxid column to create new columns
    values_from = cumsum,  # Use the cumsum column as values
    values_fill = 0        # Fill missing values with 0
  ) %>% 
  group_by(.id, mapper, serie, human_reads, total, non_human, legend) %>%
  summarise(
    FP = sum(FP, na.rm = TRUE),
    TP = sum(TP, na.rm = TRUE),
    TN = sum(TN, na.rm = TRUE),
    FN = sum(FN, na.rm = TRUE),
    .groups = "drop"
  )
reshaped_df$precision = reshaped_df$TP / (reshaped_df$TP + reshaped_df$FP)
reshaped_df$recall = reshaped_df$TP / (reshaped_df$TP + reshaped_df$FN)
reshaped_df$F1 = reshaped_df$TP / (reshaped_df$TP + 0.5*(reshaped_df$FN + reshaped_df$FP))


reshaped_df$mapper = factor(reshaped_df$mapper, levels = c("Truth",mapper), ordered = TRUE)


library(ggbreak)

ggplot(reshaped_df, aes(x=mapper, y=recall, fill = mapper))+
  geom_boxplot()+
  #ylab("# of Correctly Mapped reads against Human Genomes") + 
  xlab("Mapping steps") + ylab("Mapping Recall score at each step of classification") +
  #scale_fill_manual(values = custom_colors) +
  theme_linedraw()+
  labs(fill = "") +
  #coord_cartesian(ylim = c(.97,.975))+
  theme(
    axis.text.x = element_text(angle = 30, size = 10, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    panel.grid.major.x = element_line(color = "black"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.key.size = unit(0.75, 'cm'),
    legend.key.height = unit(.75, 'cm'), 
    legend.key.width = unit(.75, 'cm'),
    legend.title = element_text(size=14, face="bold"), 
    legend.text = element_text(size=10)) + 
  facet_wrap(~ legend,nrow = 2, scales = "free_x")











library(ggbreak)

ggplot(count_reads[count_reads$taxid=="TP",], aes(x=mapper, y=readname, fill = mapper))+
  geom_boxplot() + 
  ylab("# of Correctly Mapped reads against Human Genomes")+
  scale_y_break(c(100, 2600), scales = c(5,1)) +
  facet_wrap(~ taxid, scales = 'free_y', nrow = 2)

ggplot(count_reads[count_reads$taxid=="TP",], aes(x=mapper, y=cumsum, fill = mapper))+
  geom_boxplot() + 
  ylab("# of Correctly Mapped reads against Human Genomes")+
  coord_cartesian(ylim = c(0,75000))
  #scale_y_break(c(100, 2600), scales = c(5,1)) +
  facet_wrap(~ taxid, scales = 'free_y', nrow = 2)

ggplot(count_reads[count_reads$taxid=="FP",], aes(x=mapper, y=readname, fill = mapper))+
  geom_boxplot() + 
  ylab("# of Incorrectly Mapped reads against Human Genomes")+
  facet_wrap(~ taxid, scales = 'free_y', nrow = 2)


# Pivot taxid values into separate columns (TP / FP)
count_reads_wide <- count_reads %>%
  pivot_wider(names_from = taxid, values_from = readname, values_fill = 0)

count_reads_wide$total_mapped=count_reads_wide$FP + count_reads_wide$TP
count_reads_wide$precision = count_reads_wide$TP / count_reads_wide$total_mapped


ggplot(count_reads_wide, aes(x=mapper, y=precision, fill = mapper))+
  geom_boxplot() + 
  ylab("# of Incorrectly Mapped reads against Human Genomes")
  #facet_wrap(~ taxid, scales = 'free_y', nrow = 2)

 



