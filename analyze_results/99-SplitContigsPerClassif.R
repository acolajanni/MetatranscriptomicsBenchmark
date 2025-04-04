#!/shared/ifbstor1/software/miniconda/envs/r-4.2.3/bin/Rscript
library(stringr)
library(plyr)


args = commandArgs(trailingOnly=TRUE)

#path = "/run/user/1001/gvfs/sftp:host=core.cluster.france-bioinformatique.fr,user=acolajanni/shared/projects/microbiome_translocation/"

path = "/shared/projects/microbiome_translocation/"

setwd(path)

SRA_id = args[1]
result_dir = args[2]
WhatToDo = args[3]


#SRA_id="SRR14418854"
#result_dir=paste0(path,"/results/Douek_Cleveland/Contigs/",SRA_id,"/")
#WhatToDo="Split"



result_dir=paste0(result_dir,SRA_id,"/")
contigs_classif=read.csv(file = paste0(result_dir,"/Contig_classification.csv"))
dir.create(paste0(result_dir,"contigs_classification"),showWarnings=FALSE)

if (WhatToDo=="Split"){

  for (sk in unique(contigs_classif$superkingdom)){
    
    print(sk)
    
    current_classif=contigs_classif[contigs_classif$superkingdom == sk,]
    current_dir=paste0(result_dir,"contigs_classification/")
    dir.create(current_dir,showWarnings = FALSE)
    writeLines(unique(current_classif$contig_ID), paste0(current_dir,sk,"Contigs.txt"))
    
    
    for (phy in unique(current_classif$phylum)){
      
      
      phy_classif=current_classif[current_classif$phylum == phy,]
      current_dir=paste0(result_dir,"contigs_classification/",sk,"/")
      dir.create(current_dir,showWarnings = FALSE)
      
      phy=str_replace_all(phy," ","_")
      
      writeLines(unique(phy_classif$contig_ID), paste0(current_dir,phy,"Contigs.txt"))
      
      
      
    }
    
  }
} else if (WhatToDo=="Merge"){
  
  ### merge reads ==> contigs ==> classif
  
  classif_reads=read.table(paste0(result_dir,"contigs_classification/ReadsToContigs.txt"))
  merged_classif=unique(merge(classif_reads, contigs_classif[,colnames(contigs_classif) %in% c("contig_ID", "superkingdom","phylum","class")], by.x="V2", by.y = "contig_ID"  ))
  
  
  write.table(merged_classif, file=paste0(result_dir,"contigs_classification/ReadsClassif.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")

}

stop()
