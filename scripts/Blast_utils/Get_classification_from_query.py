################################################################################                                                                           
# > september 2023                                                                                                             
# > Script : Get_classification_from_query.py                                                                                                        
# > Function : After blast/kaiju query, filter the most likely taxon for each contigs/reads                               
# @ COLAJANNI Antonin                                                          
################################################################################

import numpy as np
import pandas as pd
from ete3 import NCBITaxa
import os, re, os.path
import glob
import argparse

ncbi = NCBITaxa()

def lineage_from_taxID(taxID):
    lineage = ncbi.get_lineage(taxID)
    taxonomy_dict=ncbi.get_taxid_translator(lineage)
    tax_levels = np.vectorize(ncbi.get_rank(taxonomy_dict).get)(lineage).tolist()
    human_lineage_name = np.vectorize(taxonomy_dict.get)(lineage).tolist()
    
    tax_df = pd.DataFrame([human_lineage_name,lineage,tax_levels])
    tax_lvl_to_keep = np.logical_and(tax_df.loc[2] != "no rank",  tax_df.loc[2] != "clade").to_list()
    tax_df = tax_df[tax_df.columns[tax_lvl_to_keep]]
    tax_df.columns = tax_df.loc[2]
    tax_df.index = [taxID, "taxID","rank"]    
    return(tax_df.loc[taxID])

def get_lineage_df_old(taxID_list):
    taxonomy_df = pd.Series()
    
    for index, taxID in enumerate(taxID_list):
        if taxID == 0: # in theory, not possible causes errors
            continue
        
        if len(taxonomy_df) == 0 : 
            taxonomy_df = lineage_from_taxID(taxID)
        else: 
            tmp = lineage_from_taxID(taxID)
            taxonomy_df=pd.concat([taxonomy_df,tmp],axis=1)
            
    #col_to_remove = ["species group", "subclass" , "superfamily","subfamily", "infraclass" ,
    #                 "cohort" , "tribe" , "subgenus" ,"subcohort" , "subkingdom", 
    #                 "infraorder","parvorder","subphylum","superclass","superorder",
    #                 "suborder", "subspecies"]
    
    # ATTENTION : Si un niveau taxonomique n'est pas retrouvé dans les données au moins une fois, peut causer une erreur
    
    col_to_keep = ["superkingdom","kingdom","phylum","class","order",
                   "family","genus","species",
                   "strain",
                   "staxids"]
    taxonomy_df = taxonomy_df.copy()
    
    #taxonomy_df = taxonomy_df.drop(col_to_remove, axis = 0)
    taxonomy_df = taxonomy_df.transpose()
    taxonomy_df["staxids"] = taxonomy_df.index.to_list()
    taxonomy_df = taxonomy_df[col_to_keep]
    return(taxonomy_df)

def get_lineage_df(taxID_list):
    taxonomy_df = pd.Series()
    
    for index, taxID in enumerate(taxID_list):
        if taxID == 0: # in theory, not possible causes errors
            continue
        
        if len(taxonomy_df) == 0 : 
            taxonomy_df = lineage_from_taxID(taxID)
        else: 
            # In case we have duplicate index (ex: taxID 1027253)
            tmp = lineage_from_taxID(taxID)
            tmp = tmp.to_frame()
            tmp["index_test"] = tmp.index
            tmp = tmp.drop_duplicates(subset="index_test",keep='first').set_index("index_test")
            
            taxonomy_df=pd.concat([taxonomy_df,tmp],axis=1)
    
    col_to_keep = ["superkingdom","kingdom","phylum","class","order",
                   "family","genus","species","staxids"]
    taxonomy_df = taxonomy_df.copy()
    
    #taxonomy_df = taxonomy_df.drop(col_to_remove, axis = 0)
    taxonomy_df = taxonomy_df.transpose()
    taxonomy_df["staxids"] = taxonomy_df.index.to_list()
    
    # If no strain is found, needs to handle this case otherwise, it causes errors
    if "strain" in taxonomy_df.columns:
        col_to_keep.append("strain")
        taxonomy_df = taxonomy_df[col_to_keep]

    else:
        taxonomy_df = taxonomy_df[col_to_keep]
        taxonomy_df["strain"] = np.nan
    
    taxonomy_df = taxonomy_df[
        ["superkingdom","kingdom","phylum","class","order",
         "family","genus","species","strain", "staxids"]
    ]
 
    return(taxonomy_df)

def get_lineage_df_kaiju(taxID_list):
    taxonomy_df = pd.Series()    
    taxID_unique = list(set(taxID_list))
    
    for index, taxID in enumerate(taxID_unique):
        
        if taxID ==0 :
            continue

        if len(taxonomy_df) == 0 : # If dataframe empty creating it
            taxonomy_df = lineage_from_taxID(taxID)
        else:  # Else, add rows
            tmp = lineage_from_taxID(taxID)
            taxonomy_df=pd.concat([taxonomy_df,tmp],axis=1)       
            
    
    if 0 in taxID_unique : 
        taxonomy_df['0'] = np.nan
        
    col_to_keep = ["superkingdom","kingdom","phylum","class","order","family","genus","species","strain","staxids"]
    taxonomy_df = taxonomy_df.copy()
    
    taxonomy_df = taxonomy_df.transpose()
    taxonomy_df["staxids"] = taxonomy_df.index.to_list()
    taxonomy_df = taxonomy_df[col_to_keep]
    return(taxonomy_df)

def get_query_result_from_Blast(path_to_query):
    hit = pd.read_csv(path_to_query, header=None, sep="\t",na_values='N/A', comment = "#")
    hit = hit.fillna(0)
    #taxID=hit[5].to_list()
    taxID=hit[9].to_list()
    taxID = [int(str(item).split(";")[0]) for item in taxID]
    #hit[5] = taxID
    taxid2name = ncbi.get_taxid_translator(taxID)
    hit[len(hit.columns)+1] = np.vectorize(taxid2name.get)(taxID)
    #hit = hit.drop([7,8,9,10],axis=1)
    #hit.columns = ["contig_ID","subject_ID","evalue","pident","subject_title","staxids","subject_species"]
    hit.columns = ["contig_ID","subject_ID","evalue","pident","qcov_hsp","subject_acc","subject_acc", "align_len","subject_title","staxids","subject_species"]

    # necesary because if one ids, we have a Pandas.Series and not a Pandas.DataFrame, causing an error
    if len( [*set(hit["staxids"])] ) == 1 : 
        ids = [*set(hit["staxids"])]*2
    else:
        ids = [*set(hit["staxids"])]
        
    lineage_df = get_lineage_df( ids )
    lineage_df = lineage_df.fillna(0)

    hit = hit.merge(lineage_df, left_on='staxids', right_on='staxids', how='left') 
    return(hit.drop_duplicates())

def get_query_result_from_Blast_v2(path_to_query):
    hit = pd.read_csv(path_to_query, header=None, sep="\t",na_values='N/A', comment = "#",low_memory=False)
    hit = hit.fillna(0)
    col_taxid=8
    taxID = list(set(hit[col_taxid]))    
    taxID = [str(i) for i in taxID]


    taxID_multiple = [i for i in taxID if ";" in str(i)]

    for elem in taxID_multiple:
        # slice dataframe
        tmp = hit[hit[col_taxid] == elem]
        hit = hit[hit[col_taxid] != elem]
        tmp2 = tmp.copy()
        # split the element with ';' in the middle
        splited = elem.split(";")
        #print(splited)
        tmp2.loc[:,col_taxid] = splited[0]
        tmp.loc[:,col_taxid] = splited[1]
        # duplicate the entry, one with each of the 2 elements separeted with the ';'
        hit = pd.concat([hit, tmp, tmp2])
        
    taxID = hit[col_taxid].to_list()
    if isinstance(taxID[0],float):
        taxID = [int(i) for i in taxID]
    taxID = [int(str(item)) for item in taxID ]

    
    taxid2name = ncbi.get_taxid_translator(taxID)
    hit['species_name'] = np.vectorize(taxid2name.get)(taxID)
    #hit = hit.drop([7,8,9,10],axis=1)
    #hit.columns = ["contig_ID","subject_ID","evalue","pident","subject_title","staxids","subject_species"]
    hit.columns = ["contig_ID","subject_ID","evalue","pident","qcov_hsp","subject_acc", "align_len","subject_title","staxids","titles","subject_species"]

    # necesary because if one ids, we have a Pandas.Series and not a Pandas.DataFrame, causing an error
    if len( [*set(hit["staxids"])] ) == 1 : 
        ids = [*set(hit["staxids"])]*2
        lineage_df = get_lineage_df( ids )
        lineage_df = lineage_df.fillna(0)
        hit = hit.merge(lineage_df, left_on='staxids', right_on='staxids', how='left') 
        return(hit.drop_duplicates())

    else:
        ids = [*set(hit["staxids"])]
        
    lineage_df = get_lineage_df( ids )
    lineage_df = lineage_df.fillna(0)
    hit = hit.merge(lineage_df, left_on='staxids', right_on='staxids', how='left') 
    return(hit)

def get_query_result_from_Blast_clean(path_to_query):
    hit = pd.read_csv(path_to_query, header=None, sep="\t",na_values='N/A',low_memory=False)
    hit = hit.fillna(0)
    col_taxid=8
    taxID = list(set(hit[col_taxid]))    
    taxID = [str(i) for i in taxID]

    taxID_multiple = [i for i in taxID if ";" in str(i)]

    for elem in taxID_multiple:
        # slice dataframe
        tmp = hit[hit[col_taxid] == elem]
        hit = hit[hit[col_taxid] != elem]
        tmp2 = tmp.copy()
        # split the element with ';' in the middle
        splited = elem.split(";")
        #print(splited)
        tmp2.loc[:,col_taxid] = splited[0]
        tmp.loc[:,col_taxid] = splited[1]
        # duplicate the entry, one with each of the 2 elements separeted with the ';'
        hit = pd.concat([hit, tmp, tmp2])
        
    hit.columns = ["contig_ID","subject_ID","evalue","pident","qcov_hsp","subject_acc", "align_len","subject_title","staxids","titles"]
    return(hit[["contig_ID","evalue","pident","qcov_hsp","align_len","staxids"]])


def filter_reads_from_query(query_result):
    contig_list = [*set(query_result["contig_ID"])]
    curated_df = pd.DataFrame(columns=query_result.columns.to_list())
    discard_df = pd.DataFrame(columns=query_result.columns.to_list())
    
    for count,ID in enumerate(contig_list):

        # Order by evalue to take the first (closest evalue to 0)
        tmp_df = query_result[query_result["contig_ID"] == ID].head(3) 
        tmp_df = tmp_df.sort_values(by=["contig_ID","evalue"])
        
        ## First filter when the 3 first prediction are the same / similar
        if len(set(tmp_df["class"] )) == 1 and tmp_df["superkingdom"].to_list()[0] == "Eukaryota" :
            # Discard plantae and mammals
            if tmp_df["class"].to_list()[0] == "Mammalia" or tmp_df["kingdom"].to_list()[0] == "Viridiplantae" :
                discard_df = pd.concat([discard_df, tmp_df.head(1)])
             
            # we capture the rest of eukaryotes
            else:
                curated_df = pd.concat([curated_df, tmp_df.head(1)])

        # if synthetic construct : discard
        elif len(set(tmp_df["subject_species"] )) == 1 and tmp_df["subject_species"].to_list()[0] == 'eukaryotic synthetic construct' :
            discard_df = pd.concat([discard_df, tmp_df.head(1)])
        
        elif  len(set(tmp_df["species"] )) == 1 :
            curated_df = pd.concat([curated_df, tmp_df.head(1)])
        
        elif  len(set(tmp_df["genus"] )) == 1 :
            curated_df = pd.concat([curated_df, tmp_df.head(1)])
            
        elif  len(set(tmp_df["superkingdom"] )) == 1 and tmp_df["superkingdom"].to_list()[0] == "Bacteria":
            tmp_df = tmp_df[tmp_df.genus != 0]
            curated_df = pd.concat([curated_df, tmp_df.head(1)])
        
        elif  len(set(tmp_df["superkingdom"] )) == 1 and tmp_df["superkingdom"].to_list()[0] == "Viruses":
            tmp_df = tmp_df[tmp_df.genus != 0]
            curated_df = pd.concat([curated_df, tmp_df.head(1)])
        

    tmp_df = tmp_df.head(1) 
    
    if tmp_df["superkingdom"].to_list()[0] == "Bacteria":
        curated_df = pd.concat([curated_df, tmp_df.head(1)])

    elif tmp_df["superkingdom"].to_list()[0] == "Viruses":
        curated_df = pd.concat([curated_df, tmp_df.head(1)])

    elif tmp_df["superkingdom"].to_list()[0] == "Eukaryota" and tmp_df["class"].to_list()[0] != "Mammalia" and tmp_df["kingdom"].to_list()[0] != "Viridiplantae" :
        curated_df = pd.concat([curated_df, tmp_df.head(1)])
    
    else : 
        discard_df = pd.concat([discard_df, tmp_df.head(1)])
      
    discard_df = pd.concat([discard_df, curated_df[curated_df['subject_species'] == "Human ORFeome Gateway entry vector"] ])          
    curated_df = curated_df[curated_df['subject_species'] != "Human ORFeome Gateway entry vector"]
    return( curated_df, discard_df)

def create_dir(directory):
    # checking if the directory demo_folder 
    # exist or not.
    if not os.path.exists(directory):
        # if the demo_folder directory is not present 
        # then create it.
        os.makedirs(directory)
    return(directory)

def reverse_list(arr):
    left = 0
    right = len(arr)-1
    while (left < right):
        # Swap
        temp = arr[left]
        arr[left] = arr[right]
        arr[right] = temp
        left += 1
        right -= 1
 
    return arr




################# Parser

parser = argparse.ArgumentParser()
parser.add_argument('--dir', type=str, 
                    help='Specify the directory')

parser.add_argument('--FileName', type=str, 
                    help='Specify the exact name of the query.txt file (blast output)')

parser.add_argument('--QueryType', type=str,
                    help="select one of 'kaiju', 'blast' or translate_taxID ")

parser.add_argument('--OutName', type=str,
                    help="Name of the output file") 






args=parser.parse_args()




os.chdir(args.dir)


#dir_path = r'/shared/projects/microbiome_translocation/data/Douek_cell2021/Contigs/SRR**/*_Blast_query.txt'
#path = "/shared/projects/microbiome_translocation/"

#dir_path = r'/run/user/1001/gvfs/sftp:host=core.cluster.france-bioinformatique.fr,user=acolajanni/shared/projects/microbiome_translocation/data/Douek_cell2021/Contigs/SRR**/*_Blast_query.txt'
#path = "/run/user/1001/gvfs/sftp:host=core.cluster.france-bioinformatique.fr,user=acolajanni/shared/projects/microbiome_translocation/"


path = args.dir
filename=args.FileName
outname=args.OutName


file_path = f'{path}{filename}'

os.chdir(path)


if args.QueryType == "blast":
    
    #dir_path = r'./**/**_Blast_query*.txt'
    #dir_path = r'./SRR**/*_Blast_query.txt'
    #dir_path = r'./SRR**/*_MegaBlast_query*.txt'
    #dir_path = r'./SRR14418929/test_format/SRR14418929_Blast_query_test_format.txt'
    
    #for file in reverse_list(glob.glob(dir_path, recursive=True)):
        
    filename = os.path.basename(file_path)
    SRA_id = filename.split('_Blast_query_VRC_no_comment.txt')[0] 
    print(" \n", SRA_id, " \n")
        
    if os.stat(file_path).st_size == 0 : 
        print("file is empty")
        exit()

        
    #dir_result = os.path.dirname(file)
    #print(dir_result)

    #query = get_query_result_from_Blast_v2(file_path)
    query = get_query_result_from_Blast_clean(file_path)
    query.staxids=query.staxids.astype(int)
    query.to_csv(f'{path}/{SRA_id}_{outname}.tsv', sep='\t', index=False) 


    ####### Last version: No fitlering at this step (no execution of bellow lines) #####

    #query = pd.read_csv(f'{path}/{SRA_id}_{outname}.csv', index_col=0)   
        
    #curated_query, discared_query = filter_reads_from_query(query)
        

    #discared_query.to_csv(path+"/"+SRA_id+"_discarded_query.csv")    
    #with open(path+"/"+SRA_id+"_contigs_filtered_id.txt", 'w') as file:    
    #        file.writelines( list( "%s\n" % item for item in curated_query["contig_ID"].to_list() ) )
        
        
        
#dir_path = r'/shared/projects/microbiome_translocation/results/Douek_cell2021/SRR**/*_kaiju.out'
#path = "/shared/projects/microbiome_translocation/"

elif args.QueryType == "kaiju":

    dir_path = r'./SRR**/*_kaiju.out'

    for file in reverse_list(glob.glob(dir_path, recursive=True)):
        if os.stat(file).st_size == 0 : 
            print("file is empty")
            continue
        
        filename = os.path.basename(file)
        SRA_id = filename.split('_')[0] 
        print(" \n", SRA_id, " \n")
        
        dir_result = os.path.dirname(file)
        
        query = pd.read_table(file,header=None)
        query.columns = ["state","contig_ID","staxids"]
        taxID_list = query['staxids'].to_list()
        lineage_unique = get_lineage_df_kaiju(taxID_list)
        lineage_full = pd.merge(query, lineage_unique )

        lineage_full.to_csv(dir_result+"/"+SRA_id+"_kaiju_lineage.csv") 

        
        
elif args.QueryType == "translate_taxID":
    
    filename = os.path.basename(file_path)
        
    if os.stat(file_path).st_size == 0 : 
        print("file is empty")
        exit()

