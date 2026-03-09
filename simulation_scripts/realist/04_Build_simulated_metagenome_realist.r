
library(stringr)

################################################################################

#                              F U N C T I O N S                               #

################################################################################


################################################################################

#                                   C O D E                                    #

################################################################################
args = commandArgs(trailingOnly=TRUE)


project_name = "Simulation/realist/"
dir.create(paste0(path, 'data/',project_name,'transcripts/'), showWarnings = FALSE)
data_dir=paste0(path, 'data/',project_name)

load( "~//database_clean/selected_genomes/selected_genomes_list.RData")

path_parameters = "~//database_clean/selected_genomes/realist/"
transcripts_to_draw=read.delim(paste0(path_parameters,"transcript_to_draw.tsv"), sep="\t")
row.names(transcripts_to_draw) = sapply(c(0:49), toString)
  
path = "~//"
path_database_short=paste0(path,'database_clean/transcript_database/realist/')
n_replicate=as.numeric(nrow(transcripts_to_draw))
colnames(transcripts_to_draw) = str_remove(colnames(transcripts_to_draw), "g_")

for (replicate_number in c(0: (n_replicate-1) ) ) {
  set.seed(replicate_number)
  print(replicate_number)
  
  
  
  path_database=paste0(path_database_short,replicate_number,'/dedup_transcripts_ID/')
  
  for (taxon in colnames(transcripts_to_draw)) {
    if ( taxon == "Homo"){
      transcripts_database=read.delim(paste0(path_database_short,"Homo_sapiens_rna_boolaray.tsv"))
      avail_transcript=transcripts_database[transcripts_database$boolarray , ]
    }else{
      transcripts_database=read.delim(paste0(path_database,taxon,"_boolaray.tsv"))
      avail_transcript=transcripts_database[transcripts_database$boolarray , ]
    }
    

    
    n_transcript=as.numeric(transcripts_to_draw[toString(replicate_number) , taxon])
    
    if(n_transcript == 0){ next }
    
    # sample with replacement if we want more transcript than available
    if (n_transcript > nrow( avail_transcript )) {
      
      #print("random sampling with replacement")
      # Count the number of rows where boolarray is TRUE
      n_rows <- nrow(avail_transcript)
      selected_transcripts = avail_transcript[ sample(nrow(avail_transcript), n_transcript, replace = n_rows < n_transcript) , ]
      
    }else {
      #print("random sampling")
      selected_transcripts = avail_transcript[ sample(nrow(avail_transcript), n_transcript) , ]
    }
    
    #writeLines(selected_transcripts$TranscriptID, paste0(current_dir,p,"_",n_transcript,".txt") )
    current_dir=paste0(data_dir, 'transcripts/',replicate_number ,"/")
    dir.create(current_dir, showWarnings = FALSE)
    
    writeLines(selected_transcripts$TranscriptID, paste0(current_dir,taxon,"_all.txt") )
    
    
    
    selected_transcripts_replicates=as.data.frame(table(selected_transcripts$TranscriptID))
    selected_transcripts_replicates$Var1=selected_transcripts_replicates$Var1 <- sub("_[^_]+$", "", selected_transcripts_replicates$Var1)
    colnames(selected_transcripts_replicates) = c("TranscriptID","Freq")
    write.table(selected_transcripts_replicates, file=paste0(current_dir,taxon,"_selected_transcripts.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
  }
  

}











