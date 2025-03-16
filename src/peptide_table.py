"""" 
   peptide_table.py                 

   All basic manipulations related to peptide tables:
  parse_peptide_table, parse_peptide_tables, build_peptide_table_from_set_of_markers

"""              

import csv
import shutil
import sys
import re
from functools import cmp_to_key, partial

from src import markers as ma
from src import limit as lim
from src import message
from src import utils
from src import config

def rename_field(field):
        clean={
            "GENE":"GN",
            "PROTEIN":"GN",
            "TAXID":"OX",
            "TAXON":"OS",
            "TAXON NAME":"OS",
            "SEQID":"SeqID"
                }
        if field.upper() in clean:
            return clean[field.upper()]
        else:
            return field
            
def restitute_field(field):
    restitute={
        "OX": "TaxID",
        "OS": "Taxon name",
        "GN": "Gene"
        }
    if field in restitute:
        return restitute[field]
    else:
        return field
    

def process_fields_of_a_row(row):
    clean_row={}
    for key, value in row.items():
        clean_key=utils.clean(key)
        clean_value=utils.clean(value)
        if clean_key is None or clean_value is None:
            continue
        clean_key=rename_field(clean_key)
        clean_row[clean_key]=clean_value
    return clean_row
        
def process_one_row_from_peptide_table(row, i, file=None, warning_on=False):
    clean_row =process_fields_of_a_row(row)
    if warning_on :
        # this one should be elsewhere
        #if "Mass" not in clean_row and "Sequence" not in clean_row:
        #message.escape("File "+file+ "(peptide table): Both peptide sequence and mass columns are missing in the peptide table. You should provide at least one of those two elements.")
        if "Mass" not in clean_row and "Sequence" not in clean_row:
            message.warning("File "+file+", line "+str(i)+": missing mass and peptide sequence. Ignored.")
            return set()
        if "OS" not in clean_row and "OX" not in clean_row:
            message.warning("File "+file+", line "+str(i)+": missing taxid and taxon name. Ignored.")
            return set()
            
        if "PTM" in clean_row and not utils.is_PTM(row.get("PTM"),{'H', 'D', 'P'}): # config
            message.warning("File "+file+", line "+str(i)+": wrong PTM, "+clean_row["PTM"]+ ". Ignored.")
            del clean_row["PTM"]
    
    try:
        if "Mass" in clean_row:
            clean_row["Mass"]=utils.floating(clean_row["Mass"])
    except ValueError:
        message.warning("File "+file+", line "+str(i)+": wrong mass, "+clean_row["Mass"]+ ". Ignored.")
        del clean_row["Mass"]
    try:
        if "Hel" in clean_row:
            clean_row["Hel"]=utils.integer(clean_row["Hel"])
    except ValueError:
        message.warning("File "+file+", line "+str(i)+": wrong helical position, "+clean_row["Hel"]+ ". Ignored.")
        del clean_row["Hel"]
    try:
        if "Length" in clean_row:
            clean_row["Length"]=utils.integer(clean_row["Length"])
    except ValueError:
        message.warning("File "+file+", line "+str(i)+": wrong peptide length, "+clean_row["Length"]+ ". Ignored.")
        del clean_row["Length"]
    try:
        if "Begin" in clean_row:
            clean_row["Begin"]=utils.integer(clean_row["Begin"])
    except ValueError:
        message.warning("File "+file+", line "+str(i)+": wrong begin position, "+clean_row["Begin"]+ ". Ignored.")
        del clean_row["Begin"]
    try:
        if "End" in clean_row:
            clean_row["End"]=utils.integer(clean_row["End"])
    except ValueError:
        message.warning("File "+file+", line "+str(i)+": wrong end position, "+clean_row["End"]+ ". Ignored.")
        del clean_row["End"]
                  
    if "SeqID" not in clean_row :
        new_marker=ma.Marker(field=clean_row)
        return {new_marker}
        
    seqids=clean_row["SeqID"].split()
    new_markers=set()
    for id in seqids:
        clean_row["SeqID"]=utils.standard(id)
        new_marker=ma.Marker(field=clean_row)
        new_markers.add(new_marker)
    return new_markers


def parse_peptide_table(peptide_table_file_name, warning_on):
    set_of_markers=set()
    peptide_table = csv.DictReader(open(peptide_table_file_name), delimiter="\t")
    list_of_headers=list(map(utils.clean, peptide_table.fieldnames))
    for i,row in enumerate(peptide_table):
       # cleaned_row = {key.replace(" ","").lower():value for key, value in row.items()}
        new_markers=process_one_row_from_peptide_table(row,i+2, peptide_table_file_name, warning_on)
        set_of_markers.update(new_markers)
    if len(set_of_markers)==0:
        message.warning("File "+peptide_table_file_name+": no valid data found.") 
    return set_of_markers, list(map(rename_field, list_of_headers))


def parse_peptide_tables(list_of_peptide_tables, list_of_constraints, taxonomy, warning_on=True):
    set_of_markers=set()
    for peptide_table  in list_of_peptide_tables:
        set_of_new_markers, list_of_headers=parse_peptide_table(peptide_table, warning_on)
        set_of_markers.update(set_of_new_markers)
    if list_of_constraints is not None and len(list_of_constraints)>0 :
        set_of_markers=lim.apply_limits(list_of_constraints, set_of_markers, taxonomy,True)
    if len(list_of_peptide_tables)>1:
        list_of_headers=[]
    return set_of_markers, list_of_headers


def marker_order(m1, m2, list_of_codes):
    if m1.taxon_name()<m2.taxon_name():
        return -1
    if m1.taxon_name()>m2.taxon_name():
        return 1
    if list_of_codes.index(m1.code())< list_of_codes.index(m2.code()):
        return -1
    if list_of_codes.index(m1.code())> list_of_codes.index(m2.code()):
        return 1
    if utils.none_str(m1.PTM())<utils.none_str(m2.PTM()):
        return -1
    if utils.none_str(m1.PTM())>utils.none_str(m2.PTM()):
        return 1
    return 0
    
    
def build_peptide_table_from_set_of_markers(set_of_markers, outfile_name, list_of_headers=None, append_file=""):
    TSV_file = open(outfile_name, "w")
    if not list_of_headers:
        headers={restitute_field(key) for m in set_of_markers for key in m.field}
        list_of_headers=(config.sort_headers(headers))
    else:
        list_of_headers=list(map(restitute_field, list_of_headers))
    writer = csv.DictWriter(TSV_file, fieldnames=list_of_headers, delimiter="\t")
    writer.writeheader()
    # ordering markers
    set_of_codes=[m.code() for m in set_of_markers]
    list_of_codes=config.sort_headers(set_of_codes)
    list_of_markers=list(set_of_markers)
    list_of_markers.sort(key=cmp_to_key(partial(marker_order, list_of_codes=list_of_codes)))
    for m in list_of_markers:
        dict={h:m.field[key] for key in m.field for  h in list_of_headers if utils.equiv(h,restitute_field(key))}
        writer.writerow(dict)
    TSV_file.close()
    
    """
    if len(append_file)==0:
        
        tsv_file = open(tsv_outfile_name, "w")
        tsv_file.write("Rank \t Taxid \t Taxon name \t Sequence \t PTM \t Marker \t Mass \t Gene \t Hel \t SeqId \t Begin \t End \t Status\t Comment\n")
    else:
        shutil.copyfile(append_file, tsv_outfile_name)
        tsv_file = open(tsv_outfile_name, "a") 
    for marker in set_of_markers:
        tsv_file.write(str(marker)+"\n")
    tsv_file.close()
    """


