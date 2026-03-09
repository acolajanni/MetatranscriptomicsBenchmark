
filepath="/shared/projects/microbiome_translocation/"
load(paste0(filepath,"database_clean/scripts/selected_genomes_list.RData"))


lapply(names(random_drawing), function(rep){
  x=random_drawing[[rep]]
  print(rep)
  current_path=paste0(filepath, "database_clean/")
  dir.create(paste0(current_path,rep),showWarnings = FALSE)
  write.table(x,paste0(current_path,rep,"/selected_genomes.tsv"), 
              sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  }