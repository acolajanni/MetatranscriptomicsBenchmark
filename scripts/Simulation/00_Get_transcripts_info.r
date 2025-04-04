#!/shared/ifbstor1/software/miniconda/envs/r-4.2.3/bin/Rscript

version[['version.string']]
library(stringr)
library(Rsubread)

path = "/run/user/1001/gvfs/sftp:host=core.cluster.france-bioinformatique.fr,user=acolajanni/shared/projects/microbiome_translocation/"
path = "/shared/projects/microbiome_translocation/"
setwd(path)


args = commandArgs(trailingOnly=TRUE)
rep = args[1]
family= args[2]
#rep=0
#family=TRUE
### Get parameters




parallel::mclapply(c(0:9) , function(rep){
  
  if (family == "FALSE") {transcript_path=paste0(path,"database_clean/transcript_database/",rep,"/")
  } else {transcript_path=paste0(path,"database_clean/transcript_database/family/",rep,"/")}
  
  
  # transcripts_database=list.files(paste0(path,"database_clean/transcript_database/dedup_transcripts_ID/"), 
  #                                 full.names = TRUE, pattern = '.fasta')
  
  transcripts_database=list.files(paste0(transcript_path,"/dedup_transcripts_ID/"), 
                                  full.names = TRUE, pattern = '.fasta')
  
for (i in transcripts_database){
  
  # phylum=str_remove(i,paste0(path,"database_clean/transcript_database/dedup_transcripts_ID/"))
  phylum=str_remove(i, transcript_path)

  phylum=str_remove(phylum,'_transcripts.fasta')
  print(phylum)

  # Scan through the fasta file to get transcript names and lengths, and selectionable transcripts
  # dir.create(paste0(data_dir,sp,"/") , showWarnings = FALSE) 
  
  ### !! CAREFULL !! renaming sequences leads to duplicated ID needs to take care of that
  ## ==> seeking them in /dedup_transcripts ID
  
  fasta=scanFasta(i)
  fasta$boolarray=(! fasta$Duplicate & fasta$Length > 200 )
  
  # write.table(fasta[,c(1,2,7)], paste0(path,"database_clean/transcript_database/",phylum,"_boolaray.tsv"),
  #       quote = FALSE, row.names = FALSE, sep = "\t" )

  write.table(fasta[,c(1,2,7)], paste0(transcript_path,phylum,"_boolaray.tsv"),
      quote = FALSE, row.names = FALSE, sep = "\t" )
  
}

},mc.cores = 10)





