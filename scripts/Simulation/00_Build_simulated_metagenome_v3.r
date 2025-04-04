
library(stringr)

################################################################################

#                              F U N C T I O N S                               #

################################################################################

################################################################################

#                                   C O D E                                    #

################################################################################
args = commandArgs(trailingOnly=TRUE)

n_transcript = args[1]
n_replicate  = args[2]
replicate_number = args[3]
project_name = args[4]
replacement = args[5]
onlyhuman = args[6]
family= args[7]

# n_transcript = 10
# n_replicate  = 10
# replicate_number = 0
# project_name = "Simulation/family/"
# replacement = "TRUE"
# onlyhuman = "FALSE"
# family = "TRUE"

n_transcript=as.numeric(n_transcript)
n_replicate=as.numeric(n_replicate)
replicate_number=as.numeric(replicate_number)

path = "/shared/projects/microbiome_translocation/"

load( "/shared/projects/microbiome_translocation/database_clean/selected_genomes/selected_genomes_list.RData")
random_drawing=lapply(random_drawing, function(x) {return(unique( x[,c("taxid","F","G")] ) ) })
taxid_taxo_key=unique(do.call(rbind , random_drawing) )
taxid_taxo_key=rbind(
  as.data.frame(list("taxid"=9606,
             "F" = "Hominidae" , 
             "G" = "Homo" )) , 
  taxid_taxo_key)

path_database_short=paste0(path,'database_clean/transcript_database/')
                     
if (family == "TRUE"){    path_database_short=paste0(path_database_short,"family/") }

path_database=paste0(path_database_short,replicate_number,'/dedup_transcripts_ID/')


transcripts_database=list.files(path_database, full.names = FALSE, pattern = '.tsv')

if (onlyhuman == TRUE ){
  transcripts_database=transcripts_database[str_detect(transcripts_database, "Homo_sapiens_rna")]
}

phylum=str_remove(transcripts_database,'_boolaray.tsv')
print(phylum)

load_df=function(x){return(read.delim(paste0(path_database,x))) }
transcripts_database <- parallel::mclapply(transcripts_database, load_df, mc.cores = 2)
names(transcripts_database) = phylum

transcripts_database[["Homo_sapiens_rna"]] = read.delim(paste0(path_database_short,"Homo_sapiens_rna_boolaray.tsv"))



dir.create(paste0(path, 'data/',project_name,'transcripts/'), showWarnings = FALSE)
data_dir=paste0(path, 'data/',project_name)
        
        
# Set of 10 replicates: select randomly n transcripts for every phylum: 
# To change: for each replicate load the corresponding "selected sequences" to draw from
set.seed(replicate_number)
print(replicate_number)

final_list=list()
for (p in phylum) {
  
  print(p)
  df=transcripts_database[[p]]
  df$taxid <- sub(".*taxid:([0-9]+)_.*", "\\1", df$TranscriptID)
  
  df=merge(df, taxid_taxo_key , by="taxid")
    
  dir.create(paste0(data_dir, 'transcripts/',replicate_number ), showWarnings = FALSE)
  current_dir=paste0(data_dir, 'transcripts/',replicate_number ,"/")
    
  
  
  if (family == "TRUE") {
    print("Family")
    
    family_file=paste0(current_dir,p,"_",n_transcript,".txt")

    df=df[df$boolarray == TRUE , ]
    if(nrow(df) == 0){
      print("lost family")
      next
    }else if(nrow(df) == n_transcript){
      writeLines(df$TranscriptID, paste0(current_dir,p,"_",n_transcript,".txt") )
      final_list[[p]]=df
      next 
    }
    
    
    n_genera=length(unique(df$G))
    n_transcript_per_G=n_transcript / n_genera
    
    # if more genus than transcript, randomly select genus
    if(n_transcript < n_genera){
      print(p)
      print("more genus than transcripts")
      df=df[ df$G %in% sample(unique(df$G), n_transcript) ,]
      n_genera=length(unique(df$G))
      n_transcript_per_G=n_transcript / n_genera
    }
    
    n_transcript_per_genus=aggregate(F ~ G,  df, length)
    colnames(n_transcript_per_genus) = c("G","n_example")
    
    # If n transcript per genus is not round, set value
    if (n_transcript %% n_transcript_per_G != 0){

      v <- rep(floor(n_transcript_per_G), n_genera)  # Base integer values
      remainder <- n_transcript - sum(v)  # Compute remainder
      if (remainder > 0) {
        v[sample(1:n_genera, remainder)] <- v[sample(1:n_genera, remainder)] + 1  # Distribute remainder
      }
    } else { v = rep(n_transcript_per_G, n_genera) }

    
    #------ Another solution is to force the repartition based on the number of transcripts there is for each genus
    
    available_counts <- table(df$G)
    v <- pmin(floor(n_transcript_per_G), available_counts)  # Ensure we don't exceed available rows
    remainder <- n_transcript - sum(v)  # Compute the remainder
    names(v) <- names(available_counts)

    # Distribute the remainder while respecting available rows
    # floor : 1.75 ==> 1 // 1.33 ==> 1
    while (remainder > 0) {
      possible_indices <- which(v < table(df$G))  # Find genera that can take more samples
      if (length(possible_indices) == 0) break  # Stop if no more room to distribute

      chosen <- sample(possible_indices, min(remainder, length(possible_indices)))  # Select random genera
      v[chosen] <- v[chosen] + 1  # Increase count for selected genera
      remainder <- n_transcript - sum(v)  # Recalculate remainder
    }
    v
    
    fam_list=list()
    for (g in unique(df$G)){
      df_G=df[df$G == g , ]
      nrow_to_sample=v[[g]]
      df_G = df_G[  sample(nrow_to_sample) , ]
      fam_list[[g]] = df_G
    }
    df_G = do.call(rbind, fam_list)
    df_G = df_G[ ! is.na(df_G$taxid), ]
    
    # If not enough rows, complete to n_transcripts
    if (nrow(df_G) < n_transcript) {
      extra_rows <- df_G[sample(1:nrow(df_G), n_transcript - nrow(df_G), replace = TRUE), ]
      df_G <- rbind(df_G, extra_rows) 
    }
    
    writeLines(df_G$TranscriptID, paste0(current_dir,p,"_",n_transcript,".txt") )
    
    final_list[[p]]=df_G
    next 
  }
  
  # If more transcripts wanted than in "transcripts.fasta"; sample with or without replacement
  if (n_transcript > nrow(df[df$boolarray==TRUE,])) {
    
    if (replacement == "TRUE") {
      print("random sampling with replacement")
      # Count the number of rows where boolarray is TRUE
      n_rows <- sum(df$boolarray == TRUE)
      # Sample with replacement based on the count
      writeLines(
        df[ sample(which(df$boolarray == TRUE), n_transcript, replace = n_rows < n_transcript) , ]$TranscriptID , 
        paste0(current_dir,p,"_",n_transcript,".txt")) }
    
    else if (replacement == "FALSE") {
      print("random sampling without replacement")
      writeLines(df[df$boolarray==TRUE,]$TranscriptID, paste0(current_dir,p,"_",n_transcript,".txt") ) }
    
  # If enough transcript, just pick random
  }else {
    df2 = df[  sample(which(df$boolarray == TRUE) , n_transcript) , ]
    writeLines(df2$TranscriptID, paste0(current_dir,p,"_",n_transcript,".txt") )
  }
 
}







