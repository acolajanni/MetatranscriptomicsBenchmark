library(dplyr)
library(MCMCpack)
library(ggplot2)


simulate_reads_per_cell_stochastic <- function(counts_matrix, reads_per_transcript = 100) {
  reads_matrix <- matrix(0,
                         nrow = nrow(counts_matrix),
                         ncol = ncol(counts_matrix),
                         dimnames = dimnames(counts_matrix))
  
  for (i in 1:nrow(counts_matrix)) {
    for (j in 1:ncol(counts_matrix)) {
      N <- counts_matrix[i, j]
      if (N > 0) {
        # each read has probability 1/reads_per_transcript to count as 1 transcript
        reads_matrix[i, j] <- max(rbinom(1, N, 1/reads_per_transcript), 1)
      } else {
        reads_matrix[i, j] <- 0
      }
    }
  }
  return(as.data.frame(reads_matrix) )
}


path="~//data/Simulation/realist/"

set.seed(123)

# Bladder cancer
df1=readxl::read_excel(paste0(path,"mbio.01607-23-s0004.xlsx"))

colnames(df1)[2:ncol(df1)] = paste0("g_",colnames(df1)[2:ncol(df1)])
colnames(df1)[1]="Sample"
# Head neck cancer
df2=readxl::read_excel(paste0(path,"mbio.01607-23-s0005.xlsx"))

# Breast cancer
df3=readxl::read_excel(paste0(path,"mbio.01607-23-s0006.xlsx"))


# Merge dataset:
all_cols <- Reduce(union, list(names(df1), names(df2), names(df3)))

df1[setdiff(all_cols, names(df1))] <- 0
df2[setdiff(all_cols, names(df2))] <- 0
df3[setdiff(all_cols, names(df3))] <- 0

df_all <- rbind(
  df1[all_cols],
  df2[all_cols],
  df3[all_cols]
)

### Average taxon representation for each dataset
# proportion rowwise
get_most_frequent_taxa=function(data){
  data[,1]=NULL
  homo_reads=data$g_Homo
  data_prop_full <- sweep(data, 1, rowSums(data), FUN = "/")
  homo_reads_prop=data_prop_full$g_Homo
 
  depth_with_human=rowSums(data)
  data[,c("g_Homo")] = NULL
  m_depth=rowSums(data)
  
  # proportion from now on
  data_prop <- sweep(data, 1, rowSums(data), FUN = "/")
  data_prop=data_prop[! is.na(data_prop[,1]) , ]
  mean_taxa=apply(data_prop, MARGIN = 2, mean)
  hist(mean_taxa)
  df=as.data.frame(mean_taxa)
  df$taxon=row.names(df)
  df=df[order(df$mean_taxa, decreasing = TRUE) , ]
  df$cum_freq <- cumsum(df$mean_taxa)
  l=nrow(df[df$mean_taxa > 0.001,])
  plot(
    df$cum_freq,
    type = "l", main = l,
    xlab = "Observation (ordered by frequency)",
    ylab = "Cumulative frequency"
  )
  abline(v = l, lty = 2, lwd = 2)
  
  taxa_prop=df[df$mean_taxa > 0.001,]

  return(list("taxa_prop"=taxa_prop,
              "otu_mat"=data[,taxa_prop$taxon],
              "depth"= depth_with_human,
              "homo_read_prop"=homo_reads_prop,
              "microbial_depth"=m_depth))
}
res=get_most_frequent_taxa(df_all)

lib_sizes=res$depth
lib_sizes2=res$microbial_depth
hist(lib_sizes, breaks = 20)
hist(lib_sizes2, breaks = 20)
summary(res$homo_read_prop)





ggplot(data.frame(x=res$homo_read_prop), aes(x = x)) +
  geom_histogram(
    bins = 20,
    color = "black",
    fill = "steelblue"
  ) +
  labs(
    title = "Distribution of Read classified as humain",
    x = "Human Read Proportion",
    y = ""
  ) +
  theme_minimal(base_size = 14)

### Set human frac to 50%
human_frac <- 0.5
bact_frac <- 1-human_frac

# Assume otu_top is 729 x 81 (samples x taxa)
otu_top=res$otu_mat
otu_mat <- as.matrix(otu_top)

# convert to relative abundances (proportions)
otu_prop <- otu_mat / rowSums(otu_mat)
otu_prop = otu_prop[!is.na(otu_prop[,1]) , ]


n_sim <- 50
simu_depth <- round(mean(lib_sizes)) * bact_frac

# estimate mean proportion per taxon
mu <- colMeans(otu_prop)  
alpha_hat = feralpack::dirichlet_params(mu, apply(otu_prop, 2, sd))


# Dirichlet draws
p_sim <- MCMCpack::rdirichlet(n_sim, alpha_hat)

# simulate counts
counts_sim_df <- t(sapply(1:n_sim, function(i) {
  rmultinom(1, simu_depth, prob = p_sim[i, ])
}))

# check totals
summary(rowSums(counts_sim))

# Add human reads to the mix
simu_depth_human <- round(mean(lib_sizes)) * human_frac
counts_sim_df$g_Homo = simu_depth_human

# Remove taxon for which no read has been simulated for them
summary(counts_sim_df)
counts_sim_df = counts_sim_df[,colSums(counts_sim_df) > 0]

# Binomial distribution with 1/100 probability of drawing 1 transcript for a read
reads_per_transcript <- 100
transcripts_sim <- simulate_reads_per_cell_stochastic(counts_sim_df, 100)

save_dir="~//database_clean/selected_genomes/realist/"

# write.table(transcripts_sim, paste0(save_dir,"transcript_to_draw.tsv"), quote=FALSE, sep="\t", row.names=FALSE)

# write.table(counts_sim_df, paste0(save_dir,"reads_to_simulate.tsv"), quote=FALSE, sep="\t", row.names=FALSE)

tmp=read.table(paste0(save_dir,"reads_to_simulate.tsv"), sep="\t", header=TRUE)


# writeLines(stringr::str_remove(colnames(counts_sim_df),"g_"), paste0(path,"genus_to_select.txt"))