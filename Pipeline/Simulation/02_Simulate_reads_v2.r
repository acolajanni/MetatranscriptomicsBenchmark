#!/shared/ifbstor1/software/miniconda/envs/r-4.2.3/bin/Rscript

version[['version.string']]
library(stringr)
library(Rsubread)
library(doParallel)

args = commandArgs(trailingOnly=TRUE)

path = "/run/user/1001/gvfs/sftp:host=core.cluster.france-bioinformatique.fr,user=acolajanni/shared/projects/microbiome_translocation/"
path = "/shared/projects/microbiome_translocation/"
# path='/home/acolajanni/Documents/work/'

setwd(path)

data_dir              = args[1]
transcript_dir        = args[2]
n_transcript          = as.numeric(args[3])
reads_per_transcripts = as.numeric(args[4])


num_cores <- 10



## Load transcripts - create final dir
# transcript_dir=paste0(path, "data/Simulation/family/transcripts/0/")
# data_dir=paste0(path, "data/Simulation/family/")
# n_transcript=10
# reads_per_transcripts=10


transcripts=list.files(path = transcript_dir, pattern = paste0("_", n_transcript,".fasta") )

dir.create(paste0(data_dir,"raw_reads_",n_transcript,"/"),showWarnings = FALSE)
output_dir=paste0(data_dir,"raw_reads_",n_transcript,"/")

ID_run=str_split(transcript_dir,fixed("/"))[[1]]
ID_run=ID_run[ID_run != ""]
ID_run=ID_run[length(ID_run)]

output_dir=paste0(output_dir,ID_run,"/")
dir.create(output_dir, showWarnings = FALSE)

set.seed(as.numeric(ID_run))


for (tr in transcripts) {
  phylum <- str_remove(tr, paste0("_", n_transcript, ".fasta"))
  print(phylum)
  
  # Create directories
  dir.create(paste0(output_dir, phylum, "/"), showWarnings = FALSE)
  dir.create(paste0(output_dir, phylum, "/transcript_reads/"), showWarnings = FALSE)
  reads_path <- paste0(output_dir, phylum, "/")
  fasta_path <- paste0(transcript_dir, tr)
  dir.create(paste0(reads_path, "/tmp_merge/"), showWarnings = FALSE)
  
  # Scan fasta file
  true_n_transcript <- nrow(scanFasta(fasta_path))
  expr <- rep(0, true_n_transcript)
  
  sim_list <- list()
  
  # Process in chunks of 100 iterations
  for (chunk_start in seq(1, n_transcript, by = 100)) {
    chunk_end <- min(chunk_start + 99, n_transcript)
    cat("Processing chunk:", chunk_start, "to", chunk_end, "\n")
    
    # Start parallel backend for this chunk
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    
    # Parallelized chunk
    chunk_results <- foreach(i = chunk_start:chunk_end, .packages = c("stringr","Rsubread")) %dopar% {
      expr[i] <- 1
      
      tmp <- simReads(
        fasta_path,
        expr,
        paste0(reads_path, "/transcript_reads/", phylum, "_", reads_per_transcripts, "rt_", i),
        library.size = reads_per_transcripts,
        read.length = 100,
        truth.in.read.names = TRUE,
        simulate.sequencing.error = TRUE,
        quality.reference = NULL,
        paired.end = TRUE,
        fragment.length.min = 100,
        fragment.length.max = 500,
        fragment.length.mean = 180,
        fragment.length.sd = 40,
        simplify.transcript.names = FALSE
      )
      
      expr[i] <- 0
      
      tmp[tmp$NReads > 0, ]
    }
    
    # Stop the parallel backend
    stopCluster(cl)
    
    # Collect chunk results
    sim_list <- c(sim_list, chunk_results)
    
    # Concatenate reads for this chunk
    reads <- list.files(path = paste0(reads_path, "/transcript_reads/"), 
                        pattern = paste0(phylum, "_", reads_per_transcripts, "rt_"), full.names = TRUE)
    R1 <- reads[str_detect(reads, "R1.fastq.gz")]
    R2 <- reads[str_detect(reads, "R2.fastq.gz")]
    
    system(paste("cat", paste(R1, collapse = " "), ">", 
                 paste0(reads_path, "/tmp_merge/", phylum, "_", reads_per_transcripts, "rt_", round(chunk_start / 100), "_R1.fastq.gz")))
    system(paste("cat", paste(R2, collapse = " "), ">", 
                 paste0(reads_path, "/tmp_merge/", phylum, "_", reads_per_transcripts, "rt_", round(chunk_start / 100), "_R2.fastq.gz")))
    
    file.remove(c(R1, R2))  # Remove temporary files
  }
  
  # Combine all results into a single data frame
  sim_list <- do.call(rbind, sim_list)
  
  # Write results
  write.table(sim_list, 
              file = paste0(reads_path, phylum, "_", reads_per_transcripts, "rt_simu_df.tsv"),
              sep = '\t', quote = FALSE, row.names = FALSE)
}
























### old
# 
# for (tr in transcripts) {
#   phylum <- str_remove(tr, paste0("_", n_transcript, ".fasta"))
#   print(phylum)
# 
#   # Create directories
#   dir.create(paste0(output_dir, phylum, "/"), showWarnings = FALSE)
#   dir.create(paste0(output_dir, phylum, "/transcript_reads/"), showWarnings = FALSE)
#   reads_path <- paste0(output_dir, phylum, "/")
#   fasta_path <- paste0(transcript_dir, tr)
#   dir.create(paste0(reads_path, "/tmp_merge/"), showWarnings = FALSE)
#   
#   # Scan fasta file
#   true_n_transcript <- nrow(scanFasta(fasta_path))
#   expr <- rep(0, true_n_transcript)
# 
#   sim_list <- list()
#   
#   # Process in chunks of 100 iterations
#   for (chunk_start in seq(1, 200, by = 100)) {
#     chunk_end <- min(chunk_start + 99, 200)
#     cat("Processing chunk:", chunk_start, "to", chunk_end, "\n")
#     
#     # Start parallel backend for this chunk
#     cl <- makeCluster(num_cores)
#     registerDoParallel(cl)
#     
#     # Parallelized chunk
#     chunk_results <- foreach(i = chunk_start:chunk_end, .packages = c("stringr")) %dopar% {
#       expr[i] <- 1
#       
#       tmp <- simReads(
#         fasta_path,
#         expr,
#         paste0(reads_path, "/transcript_reads/", phylum, "_", reads_per_transcripts, "rt_", i),
#         library.size = reads_per_transcripts,
#         read.length = 100,
#         truth.in.read.names = TRUE,
#         simulate.sequencing.error = TRUE,
#         quality.reference = NULL,
#         paired.end = TRUE,
#         fragment.length.min = 100,
#         fragment.length.max = 500,
#         fragment.length.mean = 180,
#         fragment.length.sd = 40,
#         simplify.transcript.names = FALSE
#       )
#       
#       expr[i] <- 0
#       
#       tmp[tmp$NReads > 0, ]
#     }
#     
#     # Stop the parallel backend
#     stopCluster(cl)
#     
#     # Collect chunk results
#     sim_list <- c(sim_list, chunk_results)
#     
#     # Concatenate reads for this chunk
#     reads <- list.files(path = paste0(reads_path, "/transcript_reads/"), 
#                         pattern = paste0(phylum, "_", reads_per_transcripts, "rt_"), full.names = TRUE)
#     R1 <- reads[str_detect(reads, "R1.fastq.gz")]
#     R2 <- reads[str_detect(reads, "R2.fastq.gz")]
#     
#     system(paste("cat", paste(R1, collapse = " "), ">", 
#                  paste0(reads_path, "/tmp_merge/", phylum, "_", reads_per_transcripts, "rt_", chunk_start / 100, "_R1.fastq.gz")))
#     system(paste("cat", paste(R2, collapse = " "), ">", 
#                  paste0(reads_path, "/tmp_merge/", phylum, "_", reads_per_transcripts, "rt_", chunk_start / 100, "_R2.fastq.gz")))
#     
#     file.remove(c(R1, R2))  # Remove temporary files
#   }
#   
#   # Combine all results into a single data frame
#   sim_list <- do.call(rbind, sim_list)
#   
#   # Write results
#   write.table(sim_list, 
#               file = paste0(reads_path, phylum, "_", reads_per_transcripts, "rt_simu_df.tsv"),
#               sep = '\t', quote = FALSE, row.names = FALSE)
# }
# 
