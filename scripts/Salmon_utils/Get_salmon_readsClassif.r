#!/shared/ifbstor1/software/miniconda/envs/r-4.2.3/bin/Rscript

library(stringr)
library(plyr)
library(ggsankey)
library(ggplot2)
library(data.table)


### Functions
# For each contig gives the full lineage annotation for human and unclassified
annotate_missing_contig=function(contig_ids, annotation="human"){
  if (annotation=="human"){
    lineage=c(9606 , "Eukaryota",	"Chordata",	"Mammalia",	"Primates",	"Hominidae",	
              "Homo",	"Homo sapiens", "unclassified Homo sapiens subspecies/strain")
  }else if (annotation=="unclassified"){ lineage=c(0 , rep("unclassified",8)) }
  
  columns=c("taxid","D","P","C","O","F","G","S",'T')
  
  lineage_df=data.frame(matrix(lineage, nrow = 1))
  colnames(lineage_df)=columns
  contig_list=lapply(unique(contig_ids), function(x){data.frame("contig_ID"=x)})
  
  contig_list=lapply(contig_list, function(x) {cbind(x,lineage_df) })
  return(do.call(rbind, contig_list))
  
}

### 


args = commandArgs(trailingOnly=TRUE)


path = "/run/user/1001/gvfs/sftp:host=core.cluster.france-bioinformatique.fr,user=acolajanni/shared/projects/microbiome_translocation/"
path = "/shared/projects/microbiome_translocation/"

args = commandArgs(trailingOnly=TRUE)

path_result = args[1]
SRA_id  = args[2]

# path_result=paste0(path,'results/Douek_Cleveland/Contigs/')
# SRA_id="SRR14418867"

### Salmon - Blast classif ###
# 1- Load Blast results: 

result_dir_salmon=paste0(path_result,"/",SRA_id,"/")

# Blast classification
contig_blast=read.csv2(paste0(result_dir_salmon,"Contig_classification.csv"), sep=",")
contig_blast[,c(".id","evalue","pident","qcov_hsp","align_len","classif","log_evalue","rank","taxon","kingdom")]=NULL
colnames(contig_blast)[colnames(contig_blast) == "staxids"] = "taxid"
colnames(contig_blast)[3:10]=c("D","P","C","O","F","G","S",'T')
contig_blast$subject_species = NULL

human_readstocontigs=read.table(paste0(result_dir_salmon,"ContigsToReads/human_readsToContigs.txt"), 
                        header = FALSE, sep="\t",col.names = c("read_name","contig_ID") )
unclassified_readstocontigs=read.table(paste0(result_dir_salmon,"ContigsToReads/unclassified_readsToContigs.txt"), 
                               header = FALSE, sep="\t",col.names = c("read_name","contig_ID") )

# 2: Merge them

human_contig=annotate_missing_contig(human_readstocontigs$contig_ID, annotation = "human")
unclassified_contig=annotate_missing_contig(unclassified_readstocontigs$contig_ID, annotation = "unclassified")

contig_classif=do.call(rbind, list(human_contig, unclassified_contig, contig_blast))

####
# 3- Load Salmon reads ==> Contigs

reads_to_contig=read.table(paste0(result_dir_salmon,"ContigsToReads/reads_to_contig.txt"), 
                           header = FALSE, sep="\t",col.names = c("read_name","contig_ID") )                  

reads_to_contig=do.call(rbind, list(reads_to_contig,human_readstocontigs,unclassified_readstocontigs))
reads_to_contig=unique(reads_to_contig)


# 4: LCA when multimapping
reads_to_contig=merge(reads_to_contig, contig_classif, by="contig_ID")

t=as.data.frame(table(reads_to_contig$read_name))
multimapped_reads=t[t$Freq > 1,]$Var1
multimapped_classif=reads_to_contig[reads_to_contig$read_name %in% multimapped_reads,]
reads_to_contig=reads_to_contig[! reads_to_contig$read %in% multimapped_reads,]

reads_to_contig[,c("contig_ID", "taxid")]=NULL
c("D","P","C","O","F","G","S",'T')

# What we want: for multiplampped reads, they will have the exact same classification across different contigs (contig column will be removed later)
# Check optimisation ChatGPT
# LCA <- parallel::mclapply(unique(multimapped_classif$read_name),
#                           function(r){
#                             tmp=multimapped_classif[multimapped_classif$read_name == r, ]
#                             if (length(unique(tmp$T)) != 1) {tmp$T = "unclassified" }
#                             if (length(unique(tmp$S)) != 1) {tmp$S = "unclassified" }
#                             if (length(unique(tmp$G)) != 1) {tmp$G = "unclassified" }
#                             if (length(unique(tmp$F)) != 1) {tmp$F = "unclassified" }
#                             if (length(unique(tmp$O)) != 1) {tmp$O = "unclassified" }
#                             if (length(unique(tmp$C)) != 1) {tmp$C = "unclassified" }
#                             if (length(unique(tmp$P)) != 1) {tmp$P = "unclassified" }
#                             if (length(unique(tmp$D)) != 1) {tmp$D = "unclassified" }
#                             return(unique(tmp)) },
#                           mc.cores = 20)
# LCA=do.call(rbind,LCA)
# LCA[,c("contig_ID", "taxid")]=NULL
# LCA=unique(LCA)


# Convert your data.frame to data.table
dt <- as.data.table(multimapped_classif)

# Define taxonomic levels in hierarchical order
levels <- c("D", "P", "C", "O", "F", "G", "S", "T")

# Group by 'read_name' and propagate "unclassified"
dt[, (levels) := {
  # Create a vector to track which columns need to be set to "unclassified"
  unclassified_levels <- rep(FALSE, length(levels))
  
  # Loop through taxonomic levels
  for (i in seq_along(levels)) {
    level <- levels[i]
    if (uniqueN(.SD[[level]]) > 1) {  # Check if there are inconsistencies
      unclassified_levels[i:length(levels)] <- TRUE  # Mark lower levels for update
      break
    }
  }
  
  # Apply updates only to marked levels
  lapply(seq_along(levels), function(i) {
    if (unclassified_levels[i]) "unclassified" else .SD[[levels[i]]]
  })
}, by = read_name, .SDcols = levels]

# Drop unnecessary columns
dt[, c("contig_ID", "taxid")] = NULL

# Remove duplicate rows
LCA <- unique(dt)




final=rbind(reads_to_contig, LCA)

write.table(final,paste0(result_dir_salmon,"ContigsToReads/ReadsClassif.txt") ,
            quote=FALSE, row.names = FALSE, sep="\t")

stop()

#### TEST chatGPT:
# 
# library(future)
# library(dplyr)
# plan(multisession, workers = 12)
# levels <- c("D", "P", "C", "O", "F", "G", "S", "T")
# 
# # Function to propagate "unclassified" hierarchically
# propagate_unclassified <- function(data) {
#   for (level in levels) {
#     if (n_distinct(data[[level]]) != 1) {
#       # Set the current level and all lower levels to "unclassified"
#       lower_levels <- levels[which(levels == level):length(levels)]
#       data[lower_levels] <- "unclassified"
#       break  # Stop once we set a level to "unclassified"
#     }
#   }
#   # Return only the first row of the group
#   return(data[1, ])
# }
# 
# # Group by read_name and apply the hierarchical unclassified rule
# result <- multimapped_classif %>%
#   select(-contig_ID, -taxid) %>%  # Drop unwanted columns
#   group_by(read_name) %>%
#   group_modify(~ propagate_unclassified(.x)) %>%
#   slice(1) %>%  # Ensure only the first row of each group is kept
#   ungroup()

###





final=rbind(reads_to_contig, LCA)


write.table(final,paste0(result_dir_salmon,"ContigsToReads/ReadsClassif.txt") ,
            quote=FALSE, row.names = FALSE, sep="\t")




library(microbenchmark)
library(dplyr)
library(future)
plan(multisession, workers = 20)  # Set up parallel backend
# Benchmark mclapply approach
benchmark_mclapply <- function() {
  LCA <- parallel::mclapply(unique(multimapped_classif$read_name),
                            function(r){  
                              tmp=multimapped_classif[multimapped_classif$read_name == r, ]
                              if (length(unique(tmp$T)) != 1) {tmp$T = "unclassified" }
                              if (length(unique(tmp$S)) != 1) {tmp$S = "unclassified" }
                              if (length(unique(tmp$G)) != 1) {tmp$G = "unclassified" }
                              if (length(unique(tmp$F)) != 1) {tmp$F = "unclassified" }
                              if (length(unique(tmp$O)) != 1) {tmp$O = "unclassified" }
                              if (length(unique(tmp$C)) != 1) {tmp$C = "unclassified" }
                              if (length(unique(tmp$P)) != 1) {tmp$P = "unclassified" }
                              if (length(unique(tmp$D)) != 1) {tmp$D = "unclassified" }
                              return(unique(tmp)) },
                            mc.cores = 20)
  LCA <- do.call(rbind, LCA)
  LCA[,c("contig_ID", "taxid")] <- NULL
  LCA <- unique(LCA)
}

# Benchmark dplyr approach
benchmark_dplyr <- function() {
  levels <- c("D", "P", "C", "O", "F", "G", "S", "T")
  propagate_unclassified <- function(data) {
    for (level in levels) {
      if (n_distinct(data[[level]]) != 1) {
        lower_levels <- levels[which(levels == level):length(levels)]
        data[lower_levels] <- "unclassified"
        break
      }
    }
    return(data[1, ])
  }
  result <- multimapped_classif %>%
    select(-contig_ID, -taxid) %>%
    group_by(read_name) %>%
    group_modify(~ propagate_unclassified(.x)) %>%
    slice(1) %>%
    ungroup()
}

# Compare the two methods
microbenchmark(
  mclapply = benchmark_mclapply(),
  dplyr = benchmark_dplyr(),
  times = 10
)
  



### data.table:

library(data.table)


# Convert your data.frame to data.table
dt <- as.data.table(multimapped_classif)

# Define taxonomic levels in hierarchical order
levels <- c("D", "P", "C", "O", "F", "G", "S", "T")

# Group by 'read_name' and propagate "unclassified"
dt[, (levels) := {
  # Create a vector to track which columns need to be set to "unclassified"
  unclassified_levels <- rep(FALSE, length(levels))
  
  # Loop through taxonomic levels
  for (i in seq_along(levels)) {
    level <- levels[i]
    if (uniqueN(.SD[[level]]) > 1) {  # Check if there are inconsistencies
      unclassified_levels[i:length(levels)] <- TRUE  # Mark lower levels for update
      break
    }
  }
  
  # Apply updates only to marked levels
  lapply(seq_along(levels), function(i) {
    if (unclassified_levels[i]) "unclassified" else .SD[[levels[i]]]
  })
}, by = read_name, .SDcols = levels]

# Drop unnecessary columns
dt[, c("contig_ID", "taxid")] = NULL

# Remove duplicate rows
result <- unique(dt)




microbenchmark(
  mclapply = benchmark_mclapply(),
  dplyr = benchmark_dplyr(),
  data.table = benchmark_datatable(),
  
  times = 3
)





