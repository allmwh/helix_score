import re
import json
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
from Bio.SeqRecord import SeqRecord

def seq_to_fasta(seq_id,seq,path):
    '''
    string的string做成fasta檔案

    seq: ACDERRKKRRK
    path: fasta檔案的路徑

    return: path
    '''
    
    fasta_list = []
    record = SeqRecord(Seq(seq),
                       description='',
                       id = str(seq_id))
    fasta_list.append(record)
        
    with open(path, "w") as output_handle:
        SeqIO.write(fasta_list, output_handle, "fasta")


def seq_aa_check(sequence):
    '''
    check 20 amino acid abbreviations are usual used. If not, replace with similar letter 
    see: https://www.ncbi.nlm.nih.gov/Class/MLACourse/Modules/MolBioReview/iupac_aa_abbreviations.html
    
    sequence: str, seqeunce
    '''
    amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V','X']

    sequence = sequence.replace('B', 'D').replace('Z', 'Q')
            
    return sequence

def get_uniprot_rawdata(path):
    """
    load protein data download from uniprot, the columns of uniprot must be
    "Entry, Gene names (primary), Protein names, Sequence, Organism ID" and save as "Tab-separated format (.tab)"
    uniprot: https://www.uniprot.org/uploadlists/

    path: str, the uniprot tab file path
    """
    df = pd.read_csv(path, sep="\t", names=["uniprot_id", "gene_name", "protein_name", "protein_sequence", "taxonomy"])
    df = df.drop(0).reset_index().drop(axis=1, labels="index")
    return df

def s2jgood(string):
    '''
    json清理器，json格式只能有 " ，不能有 '
    '''
    return json.loads(string.replace("\'","\""))

def get_fasta_uniprot_id(path):
    '''
    好像get_fasta_seq_info會讀比較久，所以如果不需homologous_info，就用這個
    '''
    path = Path(path)
    human_uniprot_id = path.parts[-1].split(".")[0]

    return {"human_uniprot_id":human_uniprot_id}

def get_fasta_seq_info(path):
    '''
    get fasta sequence info, include uniprot_id(fasta file name), all seqeunces tax info in fasta
    '''
    path = Path(path)
    path = path.absolute()
    fasta_list = list(SeqIO.parse(path, "fasta"))
    seq_info_list = []
    for i in fasta_list:
        sequence_info = i.description
        sequence_info = sequence_info.split('|')
        
        oma_protein_id = sequence_info[0].strip()
        species = sequence_info[1].strip()
        taxon_id = sequence_info[2].strip()
        oma_cross_reference = sequence_info[3].strip()
        
        seq_info_list.append({"oma_protein_id":oma_protein_id,
                              "species":species,
                              "taxon_id":taxon_id,
                              "oma_cross_reference":oma_cross_reference})
    
    human_uniprot_id = path.parts[-1].split(".")[0]
        
    return {"human_uniprot_id":human_uniprot_id,
            "homologous_info":seq_info_list}


def fasta_to_seqlist(path):
    '''
    輸入fasta讀成list

    path: fasta檔案

    return: list
    '''
    return list(SeqIO.parse(str(path),'fasta'))

def fasta_seq_iterator(fasta_path):
    '''
    不要一次讀出來，用generator一個個蹦出來

    fasta_path: fasta路徑

    return: generator 會一個個吐序列
    '''
    with open(fasta_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            yield str(record.seq)


def find_human_sequence(path):
    """
    find human sequence from fasta file
    """
    path = Path(path)
    all_sequence_list = fasta_to_seqlist(path)
    uniprot_id = get_fasta_seq_info(path)['human_uniprot_id']

    for i in all_sequence_list:
        tax_id = i.description.split("|")[2]
        if int(tax_id) == 9606:
            sequence = str(i.seq)
            return {"uniprot_id": uniprot_id,
                    "sequence_name": i.description,
                    "remove_gap_sequence":sequence.replace("-",""),
                    "sequence": sequence}

    # error handle for no human sequence
    raise Exception("{}, fasta path {} does not have human sequence".format(uniprot_id, str(path)))
            

def get_subset(human_df, subset):
    if subset == "rbp":
        subset_list_path = Path("./rawdata/rbp_uniprotid_list.txt")
    elif subset == "mrbp":
        subset_list_path = Path("./rawdata/mrbp_uniprotid_list.txt")
    elif subset == "human":
        return human_df
    else:
        raise ValueError("subset must be 'rbp', 'mrbp', 'human'")

    with open(subset_list_path, "r") as tf:
        subset_list = tf.read().split("\n")
    subset_df = human_df[human_df["uniprot_id"].isin(subset_list)].reset_index(drop=True)

    return subset_df


def seqlist_to_fasta(seqlist,path):
    '''
    //只用在輸入是單獨序列而已＝＝
    
    把序列的list放到fasta檔案裡面，暫時id 123...亂填

    seqlist: list[seq1,seq2...]
    path: fasta檔案的路徑

    return: path
    '''

    fasta_list = []
    for index,element in enumerate(seqlist):
        record = SeqRecord(Seq(element),
                           id = str(index))
        fasta_list.append(record)
        
    with open(path, "w") as output_handle:
        SeqIO.write(fasta_list, output_handle, "fasta")

def get_uniprot_id_from_fasta(path):
    '''
    從fastapath拿到uniprotid
    '''
    path = Path(path)
    uniprot_id = path.parts[-1].split('.')[0]
    return uniprot_id
