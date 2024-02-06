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

def parse_fasta_NCBI_header(header, gene_name, taxonomy):
    """ parsing fasta NCBI headers """
    print(header)
    new_sequence=seq.Sequence()
    re_taxon_name=re.compile('\[.*\]')
    m = re_taxon_name.search(header)
    taxon_name=m.group()
    taxon_name=taxon_name.replace('[','')
    taxon_name=taxon_name.replace(']','')
    new_sequence.taxon_name=taxon_name
    new_sequence.taxid=ta.search_taxid_from_taxon_name(taxon_name, taxonomy)
    if new_sequence.taxid==None:
        print("Missing taxon: "+taxon_name)
    new_sequence.protein=gene_name
    fields=header.split(" ") # pour le SeqID - ne marche pas avec le "vrai" format uniprot
    new_sequence.seqid=fields[0]
    return new_sequence

 
def parse_fasta_uniprot_header(header):
    """ parsing fasta uniprot headers  """
    print(header)
    new_sequence=seq.Sequence()
    re_taxid=re.compile('OX=[0-9]*\s')
    re_protein=re.compile('GN=[^\s]*\s')
    re_taxon_name=re.compile('OS=[a-zA-Z\s]*[=\n]') #en cours
    m = re_taxid.search(header)
    taxid=m.group()
    taxid=taxid.lstrip('OX=') # taxid
    taxid=taxid.strip()
    new_sequence.taxid=taxid
    m = re_taxon_name.search(header)
    taxon_name=m.group()
    taxon_name=taxon_name.replace('OS=','')
    if '=' in taxon_name:
        taxon_name=taxon_name[:-3]
    taxon_name=taxon_name.strip()
    new_sequence.taxon_name=taxon_name
    m = re_protein.search(header)
    prot=m.group()
    prot=prot.lstrip('GN=')
    prot=prot.strip()
    prot=prot.upper()
    new_sequence.protein=prot
    fields=header.split(" ") # pour le SeqID - ne marche pas avec le "vrai" format uniprot
    new_sequence.seqid=fields[0]
    return new_sequence

def build_set_of_sequences_from_fasta_file(fasta_file_name, gene_name="", taxonomy=None):
    set_of_sequences=set()
    fasta_file = open(fasta_file_name)
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        if "OS=" in seq_record.description :
            current_sequence=parse_fasta_uniprot_header(seq_record.description)
        else:
            current_sequence=parse_fasta_NCBI_header(seq_record.description,gene_name,taxonomy)
        current_sequence.sequence=str(seq_record.seq) 
        set_of_sequences.add(current_sequence)
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

def build_set_of_sequences(fasta, directory, limit, taxonomy):
    if fasta:
        set_of_sequences=build_set_of_sequences_from_fasta_file(fasta)
    elif directory:
        set_of_sequences=build_set_of_sequences_from_fasta_dir(directory)
    if limit:
        set_of_sequences=lim.apply_limits(limit, set_of_sequences, taxonomy, False)
    return set_of_sequences
            
