"""" 
   peptide_table.py                 

   All manipulations related to peptide tables: 
  parse_peptide_table, parse_peptide_tables, build_peptide_table_from_set_of_markers

"""              


import csv
import shutil
import sys
from src import markers as ma
from src import limit as lim
#from src import utils as ut

    
def parse_one_row_from_peptide_table(row):
    (rank, taxid, taxid_name, peptide_seq, ptm, code, mass, protein_name, seqid, begin_pos, end_pos,comment) = row
    if len(peptide_seq)==0 and len(mass)==0:
        raise Exception("Both petide sequence and mass are missing. You should provide at least one of those two elements.")
    ptm=ptm.replace(" ","")
    peptide_seq=peptide_seq.replace(" ","") # + passer en majuscule et controler l'alphabet
    code=code.replace(" ","")
    taxid=taxid.replace(" ","")
    mass=mass.replace(" ","")
    new_marker=ma.Marker(rank, taxid, taxid_name, peptide_seq, ptm, code, mass, protein_name, seqid, begin_pos, end_pos,comment)
    
    return new_marker

def parse_peptide_table(peptide_table_file_name):
    set_of_markers=set()
    with open(peptide_table_file_name, "r") as peptide_file:
        tsv_reader = csv.reader(peptide_file, delimiter="\t")
        # Skip the first row, which is the header
        next(tsv_reader)
        for row in tsv_reader:
            new_marker=parse_one_row_from_peptide_table(row)
            set_of_markers.add(new_marker)  
    return set_of_markers


def parse_peptide_tables(list_of_peptide_tables, limit, taxonomy):
    set_of_markers=set()
    for peptide_table  in list_of_peptide_tables:
        set_of_new_markers=parse_peptide_table(peptide_table)
        set_of_markers.update(set_of_new_markers)
    if limit:
        set_of_markers=lim.apply_limits(limit, set_of_markers, taxonomy,True)
    return set_of_markers

def build_peptide_table_from_set_of_markers(set_of_markers,tsv_outfile_name, append_file=""):
    if len(append_file)==0:
        tsv_file = open(tsv_outfile_name, "w")
        tsv_file.write("Rank \t Taxid \t Taxon name \t Sequence \t PTM \t Code \t Mass \t Gene \t SeqId \t Begin \t End \t Comment\n")
    else:
        shutil.copyfile(append_file, tsv_outfile_name)
        tsv_file = open(tsv_outfile_name, "a") 
    for marker in set_of_markers:
        tsv_file.write(str(marker)+"\n")
    tsv_file.close()



