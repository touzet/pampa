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
from os.path import join

import markers 
import sequences 
from taxonomy import *


def parse_fasta_NCBI_header(header, gene_name, taxonomy):
    """ parsing fasta NCBI headers """
    new_sequence=sequences.Sequence()
    re_taxon_name=re.compile('\[.*\]')
    m = re_taxon_name.search(header)
    taxon_name=m.group()
    taxon_name=taxon_name.replace('[','')
    taxon_name=taxon_name.replace(']','')
    new_sequence.taxon_name=taxon_name
    new_sequence.taxid=search_taxid_from_taxon_name(taxon_name, taxonomy)
    if new_sequence.taxid==None:
        print("Taxon manquant: "+taxon_name)
    new_sequence.protein=gene_name
    fields=header.split(" ") # pour le SeqID - ne marche pas avec le "vrai" format uniprot
    new_sequence.seqid=fields[0]
    return new_sequence

  
def parse_fasta_uniprot_header(header):
    """ parsing fasta uniprot headers  """
    new_sequence=sequences.Sequence()
    re_taxid=re.compile('OX=[0-9]*\s')
    re_protein=re.compile('GN=[a-zA-Z0-9]*\s')
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


def build_set_of_sequences_from_fasta_files(list_of_files):
    set_of_sequences=set()
    summary_file=open(list_of_files).read().splitlines()
    for fasta_file in summary_file:
        set_of_new_sequences=build_set_of_sequences_from_fasta_file(fasta_file)
        set_of_sequences.update(set_of_new_sequences)
    return set_of_sequences

def build_set_of_sequences_from_line_fasta_files(line, dir):
    set_of_sequences=set()
    entry= re.split('\t', line)
    entry = [item for item in entry if len(item)>0] 
    if len(entry)==1:
        files = [file for file in glob.glob(join(dir, entry[0]+"*.fasta"))]
        #print(files)
    elif len(entry)==2:
        files = [file for file in glob.glob(join(dir, entry[0]+"_"+entry[1]+"*.fasta"))]
    elif len(entry)==3:
        files= [join(dir, entry[0]+"_"+entry[1]+"_"+entry[2]+".fasta")]
    else:
        return set()
    for fasta_file in files:
        set_of_new_sequences= build_set_of_sequences_from_fasta_file(fasta_file)
        set_of_sequences.update(set_of_new_sequences)
    return  set_of_sequences
    

def build_set_of_sequences_from_fasta_dir(fasta_dir):
    set_of_sequences=set()
    fasta_files = [file for file in glob.glob(join(fasta_dir, "*.fasta"))]
    for fasta_file in fasta_files:
        file_name= os.path.join(fasta_dir, fasta_file)
        set_of_new_sequences= build_set_of_sequences_from_fasta_file(file_name)
        set_of_sequences.update(set_of_new_sequences)
    return set_of_sequences


def build_set_of_sequences_from_TSV_fasta_file(entry_file, fasta_file):
    summary_file=open(entry_file).read().splitlines()
    set_of_taxa=set()
    set_of_prot=set()
    set_of_pairs=set()
    for line in summary_file:
        e= re.split('\t', line)
        e = [item for item in e if len(item)>0] 
    if len(e)==1: # species
        set_of_taxta.add(entry[0])
    elif len(entry)==2:
        if len(entry[0]==0): #protein
            set_of_prot.add(entry[1])
        else:
            set_of_pairs=set_of_pairs.add((entry[0], entry[1]))
    else :
        print("\n problem \n")
    summary_file.close()
    set_of_sequences=set()
    f_file = open(fasta_file)
    for seq_record in SeqIO.parse(f_file, "fasta"):
        #on parse l'entête FASTA pour extraire le seqID, le taxID, la proteine et l'espèce
        if "OS=" in seq_record.description :
            current_sequence=parse_fasta_uniprot_header(seq_record.description)
        else:
            current_sequence=parse_fasta_NCBI_header(seq_record.description,GN,taxonomy)
        if (current_sequence.taxon_name).replace(" ","") in set_of_taxa or  (current_sequence.protein) in set_of_protein or ((current_sequence.taxon_name).replace(" ",""), current_sequence.protein) in set_of_pairs:
            current_sequence.sequence=str(seq_record.seq) # seq est un peptide simple
            set_of_sequences.add(current_sequence)
    f_file.close()
    return set_of_sequences
    
def build_set_of_sequences(entry, fasta):
    if not fasta:
        fasta="."
        from_dir=True
    else:
         from_dir = os.path.isdir(fasta)
    if from_dir:
        if entry:
            return build_set_of_sequences_from_TSV_fasta_dir(entry, fasta)
        else:
            return build_set_of_sequences_from_fasta_dir(fasta)
    else:
        if entry:
            return build_set_of_sequences_from_TSV_fasta_file(entry, fasta)
        else:
            return build_set_of_sequences_from_fasta_file(fasta, gene_name="",taxonomy=None)
    
