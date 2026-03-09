#!/shared/ifbstor1/software/miniconda/envs/r-4.2.3/bin/Rscript

version[['version.string']]
library(stringr)
library(Rsubread)


path = "~/"
setwd(path)


args = commandArgs(trailingOnly=TRUE)
rep = args[1]
folder= args[2]
rep=0
folder="/realist/"
### Get parameters


parallel::mclapply(c(0:49) , function(rep){
  

  transcript_path=paste0(path,"database_clean/transcript_database/",folder,rep,"/")
  
  transcripts_database=list.files(paste0(transcript_path,"/dedup_transcripts_ID/"), 
                                  full.names = TRUE, pattern = '.fasta')
  
for (i in transcripts_database){
  
  # phylum=str_remove(i,paste0(path,"database_clean/transcript_database/dedup_transcripts_ID/"))
  phylum=str_remove(i, transcript_path)

  phylum=str_remove(phylum,'_transcripts.fasta')
  print(phylum)

  
  fasta=scanFasta(i)
  fasta$boolarray=(! fasta$Duplicate & fasta$Length > 200 )
  

  write.table(fasta[,c(1,2,7)], paste0(transcript_path,phylum,"_boolaray.tsv"),
      quote = FALSE, row.names = FALSE, sep = "\t" )
}

},mc.cores = 10)





