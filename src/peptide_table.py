"""" 
   peptide_table.py                 

   All manipulations related to peptide tables: 
  parse_peptide_table, parse_peptide_tables, build_peptide_table_from_set_of_markers

"""              



import csv
import shutil
import sys
import re
from src import markers as ma
from src import limit as lim
from src import message
from src import utils 


def standard(s):
    if s==None:
        return None
    s.replace(" ","")
    if len(s)==0:
        return None
    else:
        return s

def floating(s):
    if s==None:
        return None
    s=s.replace(" ","")
    if len(s)==0:
        return None
    fl=float(s)
    if fl<0:
        raise ValueError()
    return fl

def integer(s):
    if s==None:
        return None
    s.replace(" ","")
    if len(s)==0:
        return None
    pos=int(s)
    if pos<0:
        raise ValueError()
    return pos
        
def process_one_row_from_peptide_table(row, i, file):

    if "mass" not in row and "sequence" not in row:          
        message.escape("File "+file+ "(peptide table): Both peptide sequence and mass columns are missing in the peptide table. You should provide at least one of those two elements.")

    rank=standard(row.get("rank"))
    taxid=standard(row.get("taxid"))
    if taxid==None:
         message.warning("File "+file+", line "+str(i)+": missing taxid. Ignored.")
         return set()
    taxid_name=row.get("taxonname")
    peptide_seq=standard(row.get("sequence")) # + passer en majuscule et controler l'alphabet
    code=standard(row.get("code"))
    ptm=row.get("ptm")
    if not utils.is_PTM(ptm,{'O', 'D', 'P'}):
        message.warning("File "+file+", line "+str(i)+": wrong PTM, "+row["ptm"]+ ". Ignored.")
        ptm=None 
    try:
        mass=floating(row.get("mass"))
    except ValueError:
        message.warning("File "+file+", line "+str(i)+": wrong mass, "+row["mass"]+ ". Ignored.")
        mass=None    
    protein_name=standard(row.get("gene"))
    helical=standard(row.get("hel"))
    try:
        begin_pos=integer(row.get("begin"))
    except ValueError:
        message.warning("File "+file+", line "+str(i)+": wrong begin position, "+row["begin"]+ ". Ignored.")
        begin_pos=None
    try:    
        end_pos=integer(row.get("end"))
    except ValueError:
        message.warning("File "+file+", line "+str(i)+": wrong end position, "+row["end"]+ ". Ignored.")
        end_pos=None
    comment=row.get("comment")

    if mass==None and peptide_seq==None:
         message.warning("File "+file+", line "+str(i)+": missing mass and sequence. Ignored.")
         return set()
    
    if "seqid" not in row:
        new_marker=ma.Marker(rank, taxid, taxid_name, peptide_seq, ptm, code, mass, protein_name, helical, None, begin_pos, end_pos,comment)
        return {new_marker}
    seqid=standard(row["seqid"])
    if len(seqid)==0:
         new_marker=ma.Marker(rank, taxid, taxid_name, peptide_seq, ptm, code, mass, protein_name, helical, None, begin_pos, end_pos,comment)
         return {new_marker}
    seqids=seqid.split()
    new_markers=set()
    for id in seqids:
        id=standard(id)
        new_marker=ma.Marker(rank, taxid, taxid_name, peptide_seq, ptm, code, mass, protein_name, helical, id, begin_pos, end_pos,comment)
        new_markers.add(new_marker)
    return new_markers

def parse_peptide_table(peptide_table_file_name):
    set_of_markers=set()
    peptide_table = csv.DictReader(open(peptide_table_file_name), delimiter="\t")
    for i,row in enumerate(peptide_table):
        cleaned_row = {key.replace(" ","").lower():value for key, value in row.items()}
        new_markers=process_one_row_from_peptide_table(cleaned_row,i+2, peptide_table_file_name)
        set_of_markers.update(new_markers)
    if len(set_of_markers)==0:
        message.warning("File "+peptide_table_file_name+": no valid data found.") 
    return set_of_markers

def parse_peptide_tables(list_of_peptide_tables, limit, taxonomy):
    set_of_markers=set()
    for peptide_table  in list_of_peptide_tables:
        set_of_new_markers=parse_peptide_table(peptide_table)
        set_of_markers.update(set_of_new_markers)
    if limit:
        set_of_markers=lim.apply_limits(lim.parse_limits(limit), set_of_markers, taxonomy,True)
    return set_of_markers

def build_peptide_table_from_set_of_markers(set_of_markers,tsv_outfile_name, append_file=""):
    if len(append_file)==0:
        tsv_file = open(tsv_outfile_name, "w")
        tsv_file.write("Rank \t Taxid \t Taxon name \t Sequence \t PTM \t Code \t Mass \t Gene \t Hel \t SeqId \t Begin \t End \t Comment\n")
    else:
        shutil.copyfile(append_file, tsv_outfile_name)
        tsv_file = open(tsv_outfile_name, "a") 
    for marker in set_of_markers:
        tsv_file.write(str(marker)+"\n")
    tsv_file.close()



