"""
fasta_parsing.py              
"""

from Bio import SeqIO
from Bio.SeqIO import parse 
from Bio.SeqRecord import SeqRecord
import Bio.SwissProt
import re
import glob
import os
from src import markers 
from src import sequences as seq
from src import taxonomy as ta
from src import limit as lim
from src import message

def update_field(regex, header, label,dict_sequence):
    m = regex.search(header)
    if m==None:
        raise ValueError() 
    res=m.group()
    res=res[4:]
    res=res.strip()
    if label in {"OX", "GN"}:
        res=res.upper()
    #if '=' in res:
    res=re.sub(r'\b\w+\s*=$', '', res)
    dict_sequence[label]=res

def parse_fasta_uniprot_header(header):
    """ parsing fasta uniprot headers  """
    message.debug(header)
    # Remove spaces before and after '='
    header = re.sub(r'\s+(?==)', '', header)
    header = re.sub(r'(?<==)\s+', '', header)
    dict_sequence={}
    re_taxid=re.compile(' OX=[^\s]*?(?=\s|$)',re.IGNORECASE)
    update_field(re_taxid, header, "OX", dict_sequence)
    re_protein=re.compile(' GN=[^\s]*?(?=\s|$)',re.IGNORECASE)
    update_field(re_protein, header, "GN", dict_sequence)
    re_taxon_name=re.compile('( OS=[a-zA-Z\s]*=)|( OS=[a-zA-Z\s]*?(?=$))',re.IGNORECASE)
    update_field(re_taxon_name, header, "OS", dict_sequence)
    fields=header.split(" ")
    seqid=fields[0]
    if '|' in seqid:
        re_seqid=re.compile('|[a-zA-Z0-9]*|')
        m=re_seqid.search(seqid)
        seqid=m.group()
        seqid=seqid.replace('|','')
        seqid=seqid.replace(' ','')
    dict_sequence["SeqID"]=seqid
    new_sequence=seq.Sequence(field=dict_sequence)
    return new_sequence
    

def build_set_of_sequences_from_fasta_file(fasta_file_name):
    set_of_sequences=set()
    fasta_file = open(fasta_file_name)
    is_empty=True
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        is_empty=False
        try:
            current_sequence=parse_fasta_uniprot_header(seq_record.description)
        except ValueError:
            message.warning("File "+fasta_file_name+", header "+seq_record.description+":\nBad header. Sequence ignored. See PAMPA documentation." )
            continue
        seq=str(seq_record.seq)
        if len(seq)==0 :
              message.warning("File "+fasta_file_name+", header "+seq_record.description+":\nEmpty sequence. Sequence ignored. See PAMPA documentation." )
              continue
        current_sequence.field["Sequence"]=seq
        set_of_sequences.add(current_sequence)
    if is_empty:
        message.warning("File "+fasta_file_name+" is empty.")
    fasta_file.close()
    return set_of_sequences


# prend en entrée un fichier texte qui contient une liste de fichiers au format Fasta
def build_set_of_sequences_from_fasta_files(list_of_files):
    set_of_sequences=set()
    summary_file=open(list_of_files).read().splitlines()
    for fasta_file in summary_file:
        set_of_new_sequences=build_set_of_sequences_from_fasta_file(fasta_file)
        set_of_sequences.update(set_of_new_sequences)
    return set_of_sequences

def build_set_of_sequences_from_fasta_dir(fasta_dir):
    set_of_sequences=set()
    fasta_files = [file for file in glob.glob(os.path.join(fasta_dir, "*.f*a"))]
    for fasta_file in fasta_files:
        set_of_new_sequences= build_set_of_sequences_from_fasta_file(fasta_file)
        set_of_sequences.update(set_of_new_sequences)
    return set_of_sequences

def build_set_of_sequences(fasta, directory, list_of_constraints, taxonomy):

    list_of_hard_constraints=[]
    list_of_soft_constraints=[]
        
    if len(list_of_constraints)>0:
        list_of_hard_constraints=[constraints[key] for constraints in list_of_constraints for key in constraints if key=="FileName"]
        list_of_soft_constraints=[constraints for constraints in list_of_constraints if "FileName" not in constraints]

    if fasta:
        if  len(list_of_constraints)==0 or fasta in list_of_hard_constraints:
            set_of_hard_sequences=build_set_of_sequences_from_fasta_file(fasta)
            set_of_soft_sequences=set()
        else:
            set_of_hard_sequences=set()
            set_of_soft_sequences=build_set_of_sequences_from_fasta_file(fasta)
    elif directory:
        set_of_hard_sequences=set()
        if len(list_of_constraints)==0: # TO DO: add taxonomy
            set_of_hard_sequences=build_set_of_sequences_from_fasta_dir(directory)
            set_of_soft_sequences=set()
        else:
            for file_name in list_of_hard_constraints:
                set_of_hard_sequences.update(build_set_of_sequences_from_fasta_file(os.path.join(directory, file_name)))          
            set_of_soft_sequences=build_set_of_sequences_from_fasta_dir(directory)
            set_of_soft_sequences=lim.apply_limits(list_of_soft_constraints, set_of_soft_sequences, taxonomy, False)
   
    set_of_sequences=set_of_hard_sequences | set_of_soft_sequences
    
    message.debug("constraints:" + str(len(list_of_constraints)))
    message.debug("hard sequences:" + str(len(set_of_hard_sequences)))
    message.debug("soft sequences:" + str(len(set_of_soft_sequences)))
    message.debug("soft constraints:" + str(len(list_of_soft_constraints)))
    
    return set_of_sequences
    
