#!/shared/ifbstor1/software/miniconda/envs/r-4.2.3/bin/Rscript

library(plyr)
library(dplyr)
library(reshape2)
library(stringr)

path="~/"



args = commandArgs(trailingOnly=TRUE)


result_dir          =   args[1]
sra_list            =   args[2]
lowest_rank         =   args[3]
threads             =   args[4]
use_cleanifier      =   as.logical(args[5])
dataset             =   args[6]
remove_unclasified  =   as.logical(args[7])
remove_human        =   as.logical(args[8])


# result_dir=paste0(path,"results/Simulation/realist/Contigs_rnaSpades/")
# sra_list=paste0(path,"data/Simulation/realist/patient_name.txt")
# lowest_rank = "genus"
# threads = 16
# use_cleanifier = FALSE
# dataset="Simulation/realist"
# remove_unclasified=FALSE
# remove_human=FALSE


# result_dir=paste0(path,"results/Contigs_rnaSpades/")
# sra_list=paste0(path,"data/sra_list_RNA.txt")
# lowest_rank = "genus"
# threads = 23
# use_cleanifier = FALSE
# dataset="Douek_Cleveland"
# remove_unclasified=TRUE
# remove_human=TRUE

### V A R I A B L E S ###
#dataset=str_split_fixed(result_dir, "/", 7)[1,6]
sra_list=readLines(sra_list)
# sra_metadata=read.table(paste0(path_data,"/sra_metadata.tsv"),header=TRUE)


if (basename(result_dir) %in% c("Contigs","Contigs_rnaSpades") ){ 
  type="Blast"
  
} else if ( (str_detect(basename(result_dir), "hybrid"))) {
  type="hybrid"
  
} else{type="kraken"}




### L O A D I N G ###
ranks=c("superkingdom","phylum", "class", "order", "family", "genus", "species", "strain")
ranks_to_keep = ranks[1:which(ranks == lowest_rank)]


### 1) load read per read classif
classif_ku=parallel::mclapply(sra_list, function(sra){
  ranks=c("superkingdom","phylum", "class", "order", "family", "genus", "species")
  
  if(type=="kraken"){
    x=read.table(paste0(result_dir,sra,"/ReadsClassif.txt"), header=FALSE, fill = TRUE, sep = "\t",quote = "")
    colnames(x)=c("read","taxid", "superkingdom","phylum", "class", "order", "family", "genus", "species")
    x=x[, colnames(x) %in% c("read","taxid",ranks_to_keep)]
    x[x$taxid %in% c(1,198431,155900,131567), ranks_to_keep] = "Unclassified"
    
  } else if( type == "Blast"){
    x=read.table(paste0(result_dir,sra,"/ContigsToReads/ReadsClassif.txt"), header=FALSE, fill = TRUE, sep = "\t",quote = "")
    colnames(x)=c("read","superkingdom","phylum", "class", "order", "family", "genus", "species", "strain")
    
    if(x[1,1] == "read_name"){
      x = x[2:nrow(x),] }
    
  } else if( type == "hybrid"){
    x=read.table(paste0(result_dir,sra,"/ReadsClassif.txt"), header=FALSE, fill = TRUE, sep = "\t",quote = "")
    colnames(x)=c("read","superkingdom","phylum", "class", "order", "family", "genus", "species")
    x=x[, colnames(x) %in% c("read",ranks_to_keep) ]


  }
  
  x[x==""]="Unclassified"
  
  if (remove_unclasified){
    x=x[x[[lowest_rank]] != "unclassified" ,]
    x=x[ ! str_detect(x[[lowest_rank]] , "unclassified") ,]
  } else{ # if we let unclassified in, need to make results cleaner
    cols <- rev(ranks_to_keep)
    x[cols] <- lapply(x[cols], function(x) {
      ifelse(grepl("^unclassified", x), "Unclassified", x) })
  }  
  
  return(x)
},mc.cores = threads)
names(classif_ku) = sra_list




# ### 1) a) counting reads before filtrating unclassified 
# classif_table=as.data.frame(do.call(rbind,lapply(classif_ku, function(x) nrow(x)) ))
# classif_table$Sample=row.names(classif_table)
# colnames(classif_table) = c("human_filtered_Alignment","Sample")
# 
# write.table(classif_table, file=paste0(result_dir,"human_filtered_Alignment_reads.csv"),
#             quote=FALSE, row.names=FALSE, sep=",")
# 
# ### 1) b) Filtrating unclassified 
# 
# classif_ku=parallel::mclapply(sra_list, function(x){
#   x=x[x[[lowest_rank]] != "unclassified" ,]
#   x=x[ ! str_detect(x[[lowest_rank]] , "unclassified") ,]
#   return(x)
# },mc.cores = threads)
# names(classif_ku) = sra_list

if (use_cleanifier) {
  cleanifier_path=paste0(path,"data/",dataset,"/cleaned/")
  
  classif_cleanifier=parallel::mclapply(sra_list, function(sra){
    ranks=c("superkingdom","phylum", "class", "order", "family", "genus", "species")
    x=read.table(paste0(cleanifier_path,sra,"/ReadsClassif_cleanifier.tsv"), header=FALSE, fill = TRUE, sep = "\t",quote = "")
    colnames(x)=c("read","taxid", "superkingdom","phylum", "class", "order", "family", "genus", "species")
    
    x=x[, colnames(x) %in% c("read","taxid",ranks_to_keep)]

    ### Load classif by kraken
    x2=classif_ku[[sra]]
    human_read_cleanifier=x$read
    
    ### Merge both classification
    x_filtered=x2[ ! x2$read %in% human_read_cleanifier,]
    x_filtered=rbind(x_filtered, x)
    
    table_classif=as.data.frame(table(x2[x2$read %in% human_read_cleanifier,]$genus))
    
    
    return(list(x_filtered,table_classif ))
  },mc.cores = threads)
  
  ### 1) load readstats per read classif
  
  stat_list=parallel::mclapply(sra_list, function(sra){
    x=read.table(paste0(cleanifier_path,sra,"/Counting_reads.tsv"), header=FALSE, fill = TRUE, sep = "\t",quote = "")
    df_out <- data.frame(
      .id     = sra,
      filter  = sub("\\.fastq$", "", sub(".*_", "", basename(x$V1))),
      reads   = x$V2,
      stringsAsFactors = FALSE
    )
    return(df_out)
  },mc.cores = threads)
  names(stat_list) = sra_list
  stat_list=do.call(rbind,stat_list)
  
  write.table(stat_list, file=paste0(result_dir,"cleanifier_stats.csv"),
              quote=FALSE, row.names=FALSE, sep=",")
  
  count_cleanifier_classif=lapply(classif_cleanifier, function(x) x[[2]])
  names(count_cleanifier_classif) = sra_list
  
  count_cleanifier_classif=ldply(count_cleanifier_classif)
  count_cleanifier_classif$.id = NULL
  count_cleanifier_classif=aggregate(Freq ~ Var1, count_cleanifier_classif, sum)
  
  
  
  classif_ku=lapply(classif_cleanifier, function(x) x[[1]])
  names(classif_ku) = sra_list
  
}

### 2) Register the diversity in the dataset
classification=parallel::mclapply(names(classif_ku), function(sra){
  x = classif_ku[[sra]]
  
  x$read=NULL
  x$taxid=NULL
  x=unique(x)
  return(x)
},mc.cores=threads)
names(classification) = sra_list

classification=unique(do.call(rbind, classification))

write.table(classification, file=paste0(result_dir,"classification_",lowest_rank,".csv"),
            quote=FALSE, row.names=FALSE, sep=",")


### 3) filter out some taxa
classif_ku2=parallel::mclapply(classif_ku, function(x){
  
  x$classification <- apply(x[, ranks_to_keep], 1, paste0, collapse = "|")
  t=aggregate(superkingdom ~ classification, x, length )
  colnames(t) = c("classification","count")

  t=t[!str_detect(t$classification, "unclassified cellular organisms|unclassified root superkingdom"),]
  if(remove_human){ t=t[!str_detect(t$classification, "Mammalia"),]  }
    
  t$freq = t$count / sum(t$count)
  return(t)
},mc.cores = threads)

names(classif_ku2) = sra_list

# Remove gorganvirus
if(remove_unclasified){
  Gorganvirus="Viruses|Uroviricota|Caudoviricetes|unclassified Caudoviricetes order|unclassified Caudoviricetes family|Gorganvirus"
  classif_list=lapply(classif_ku2, function(x){return(x[x$classification != Gorganvirus,])})
}else{
  classif_list=classif_ku2
}

save(classif_list, file=paste0(result_dir,"Quantification_list_",lowest_rank,".RData"))
load(paste0(result_dir,"Quantification_list_",lowest_rank,".RData"))

classif_df=ldply(classif_list)
classif_df=classif_df[,c(1:3)]
total=aggregate(count ~ .id, classif_df,sum)


write.table(classif_df, file=paste0(result_dir,"classif_df_",lowest_rank,".csv"),
            quote=FALSE, row.names=FALSE, sep=",")





if (use_cleanifier) {
  colnames(count_cleanifier_classif) = c("genus","Freq")
  #cl=unique( classification[,c("superkingdom","phylum","")]) 
  # Create the new row
  new_row <- data.frame(
    superkingdom = c("Archaea","Bacteria","Bacteria","Viruses","Eukaryota"),
    phylum = c("Promethearchaeota","Pseudomonadota","Bacillota","Kitrinoviricota","Microsporidia"),
    class = c("Promethearchaeia","Alphaproteobacteria","Bacilli","Flasuviricetes","unclassified" ),
    order=c("Promethearchaeales","Acetobacterales","Bacillales","Amarillovirales","unclassified" ),
    family=c("Promethearchaeaceae","Acetobacteraceae","Bacillaceae"," Flaviviridae","Nosematidae"),
    genus=c("Candidatus Prometheoarchaeum","Novacetimonas","Paenalkalicoccus","Pestivirus","Vittaforma" ),
    stringsAsFactors = FALSE )
  
  # Add it to the table
  cl <- rbind(cl, new_row)
  count_cleanifier_classif=unique(merge(count_cleanifier_classif, cl))
  #count_cleanifier_classif[!count_cleanifier_classif$genus %in% cl$genus,]
  count_cleanifier_classif$classif=ifelse(str_detect(count_cleanifier_classif$phylum,"myco"),"Fungi",count_cleanifier_classif$superkingdom)
  count_cleanifier_classif$classif=ifelse(str_detect(count_cleanifier_classif$phylum,"Chordata"),"Homo sapiens",count_cleanifier_classif$classif)
  
  # Plot classif at superkingdom level
  to_plot=aggregate(Freq ~ classif, count_cleanifier_classif, sum)
  to_plot$condition="cleanifier classif"
  to_plot$prop = to_plot$Freq / sum(to_plot$Freq)
  to_plot$label <- sprintf("%.2f%%", to_plot$prop * 100)
  to_plot$log_prop= log10(to_plot$prop)
  to_plot$label = paste0("n_read=",to_plot$Freq,"\n",to_plot$label)
  
  
  ggplot(to_plot, aes(x=classif,y=prop, fill=classif))+
    geom_bar(position="stack", stat="identity") + 
    geom_text(aes(label = label), size = 4 , position = position_stack(vjust = 1.05),color='black')  + 
    #geom_text(aes(label = Freq), size = 4 , position = position_stack(vjust = 1.15),color='black')  + 
    
    ggtitle("KrakenUniq classif of the human predicted reads by Cleanifier")
    
  # Zoom at Bacteria, Fungi and virus level:
  
  clean=count_cleanifier_classif[count_cleanifier_classif$genus != "Homo",]
  clean$label=ifelse(clean$Freq > 50, as.character(clean$genus), "" )
  clean = clean[,c("genus","classif","label","superkingdom","Freq")]
  clean$genus=ifelse(clean$label == "", "Other", as.character(clean$genus))
  clean=aggregate(Freq ~ ., clean, sum)
  clean$label=ifelse(clean$label == "", paste0(clean$superkingdom," - Other"), as.character(clean$label))
  clean$genus=ifelse(clean$genus == "Other", paste0(clean$superkingdom," - Other"), as.character(clean$genus))
  clean[clean$classif == "Fungi" & clean$genus == "Eukaryota - Other",]$genus = "Fungi - Other"
  
  
  clean$genus=factor(clean$genus, levels =   clean$genus[order(-clean$Freq)] )
  
  ggplot(clean, aes(x=classif,y=Freq, fill=genus))+
    geom_bar(position="stack", stat="identity") + 
    geom_text(aes(label = label), size = 4 , position = position_stack(vjust = 0.5),color='black')  + 
    #geom_text(aes(label = Freq), size = 4 , position = position_stack(vjust = 1.15),color='black')  + 
    
    ggtitle("Non-human predicted reads by KrakenUniq genus classification of human predicted reads by Cleanifier")
  
  

}
