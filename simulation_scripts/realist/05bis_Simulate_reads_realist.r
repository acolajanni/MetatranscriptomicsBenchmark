#!/shared/ifbstor1/software/miniconda/envs/r-4.2.3/bin/Rscript

version[['version.string']]
library(stringr)
library(Rsubread)
library(doParallel)

args = commandArgs(trailingOnly=TRUE)


path = "~//"
# path='/home/acolajanni/Documents/work/'

setwd(path)

ID_run                = as.numeric(args[1])
n_transcript_label = "all"
num_cores <- 10

path = "~//"
transcript_dir=paste0(path, "/data/Simulation/realist/transcripts/",ID_run,"/")
data_dir=paste0(path, "data/Simulation/realist/")

path_parameters = "~//database_clean/selected_genomes/realist/"
transcripts_to_draw=read.delim(paste0(path_parameters,"transcript_to_draw.tsv"), sep="\t")
row.names(transcripts_to_draw) = sapply(c(0:49), toString)
colnames(transcripts_to_draw) = str_remove(colnames(transcripts_to_draw), "g_")

reads_to_simulate=read.delim(paste0(path_parameters,"reads_to_simulate.tsv"), sep="\t")
row.names(reads_to_simulate) = sapply(c(0:49), toString)
colnames(reads_to_simulate) = str_remove(colnames(reads_to_simulate), "g_")



## Load transcripts - create final dir
transcript_dir=paste0(path, "/data/Simulation/realist/transcripts/",ID_run,"/")
data_dir=paste0(path, "data/Simulation/realist/")


###
n_transcript="all"
reads_per_transcripts=10
###


transcripts=list.files(path = transcript_dir, pattern = paste0(".fasta") )
transcripts=transcripts[!str_detect(transcripts,"temp")]

dir.create(paste0(data_dir,"raw_reads/"),showWarnings = FALSE)
output_dir=paste0(data_dir,"raw_reads/")

# ID_run=str_split(transcript_dir,fixed("/"))[[1]]
# ID_run=ID_run[ID_run != ""]
# ID_run=ID_run[length(ID_run)]

output_dir=paste0(output_dir,ID_run,"/")
dir.create(output_dir, showWarnings = FALSE)

set.seed(as.numeric(ID_run))

# Save the original input for later use
n_transcript_label <- n_transcript


for (tr in transcripts) {
  
  # Extract taxon name
  taxon = sub("_[^_]+$", "", tr)
  
  # retrieve fasta infos
  n_transcript = transcripts_to_draw[toString(ID_run), taxon ]
  fasta_path <- paste0(transcript_dir, tr)
  true_n_transcript <- nrow(scanFasta(fasta_path))

  # Always simulate 100 reads per transcripts
  reads_per_transcripts = 100
  
  if (n_transcript_label == "all") {
    n_transcript <- true_n_transcript
  } else {
    n_transcript <- as.numeric(n_transcript_label)
    if (n_transcript > true_n_transcript) {
      stop(paste("Requested", n_transcript,
                 "transcripts, but only", true_n_transcript,
                 "are available in", tr))
    }
  }
  
  
  # Create directories
  dir.create(paste0(output_dir, taxon, "/"), showWarnings = FALSE)
  dir.create(paste0(output_dir, taxon, "/transcript_reads/"), showWarnings = FALSE)
  reads_path <- paste0(output_dir, taxon, "/")
  dir.create(paste0(reads_path, "/tmp_merge/"), showWarnings = FALSE)
  
  expr <- rep(0, true_n_transcript)
  sim_list <- list()
  
  # Process in chunks of 100 iterations
  n_iterations=ifelse(n_transcript >= 100, 100, n_transcript)
  for (chunk_start in seq(1, n_transcript, by = n_iterations)) {
    
    cat("true_n_transcript:", true_n_transcript, "\n")
    cat("n_transcript:", n_transcript, "\n")
    cat("Length of expr:", length(expr), "\n")
    
    chunk_end <- min(chunk_start + 99, n_transcript)
    cat("Processing chunk:", chunk_start, "to", chunk_end, "\n")
    
    # Start parallel backend for this chunk
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    
    # Parallelized chunk
    chunk_results <- foreach(i = chunk_start:chunk_end, .packages = c("stringr","Rsubread")) %dopar% {
      if (i > length(expr)) return(NULL)  # avoid index errors
      
      expr[i] <- 1
      
      tmp <- simReads(
        fasta_path,
        expr,
        paste0(reads_path, "/transcript_reads/", taxon, "_", reads_per_transcripts, "rt_", i),
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
                        pattern = paste0(taxon, "_", reads_per_transcripts, "rt_"), full.names = TRUE)
    R1 <- reads[str_detect(reads, "R1.fastq.gz")]
    R2 <- reads[str_detect(reads, "R2.fastq.gz")]
    
    system(paste("cat", paste(R1, collapse = " "), ">", 
                 paste0(reads_path, "/tmp_merge/", taxon, "_", reads_per_transcripts, "rt_", round(chunk_start / 100), "_R1.fastq.gz")))
    system(paste("cat", paste(R2, collapse = " "), ">", 
                 paste0(reads_path, "/tmp_merge/", taxon, "_", reads_per_transcripts, "rt_", round(chunk_start / 100), "_R2.fastq.gz")))
    
    file.remove(c(R1, R2))  # Remove temporary files
  }
  
  # Combine all results into a single data frame
  sim_list <- do.call(rbind, sim_list)
  
  # Write results
  write.table(sim_list, 
              file = paste0(reads_path, taxon, "_", reads_per_transcripts, "rt_simu_df.tsv"),
              sep = '\t', quote = FALSE, row.names = FALSE)
  
}




