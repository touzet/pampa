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
import pandas as pd

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
            

def integer(s):
    if pd.isna(s) or s is None or s=="":
        return None
    n=int(float(s))
    if n<0:
        raise ValueError()
    return n
    
def process_fields_of_a_row(row):
    clean_row={}
    for key, value in row.items():
        clean_key=utils.clean(key)
        if value=="nan":
            clean_value=None
        else:
            clean_value=utils.clean(value)
        if clean_key is None or clean_value is None  :
            continue
        clean_key=rename_field(clean_key)
        clean_row[clean_key]=clean_value
    return clean_row
   
def check_marker(row, index, file=None, warning_on=False):
    clean_row = process_fields_of_a_row(row)
    if warning_on :
        # this one should be elsewhere
        # if "Mass" not in clean_row and "Sequence" not in clean_row:
        # message.escape("File "+file+ "(peptide table): Both peptide sequence and mass columns are missing in the peptide table. You should provide at least one of those two elements.")
        if "Mass" not in clean_row and "Sequence" not in clean_row:
            message.warning("File "+file+", line "+str(index)+": missing mass and peptide sequence. Ignored.")
            return set()
        if "OS" not in clean_row and "OX" not in clean_row:
            message.warning("File "+file+", line "+str(index)+": missing taxid and taxon name. Ignored.")
            return set()
            
        if "PTM" in clean_row and not utils.is_PTM(row.get("PTM"),{'H', 'D', 'P'}): # config
            message.warning("File "+file+", line "+str(index)+": wrong PTM, "+clean_row["PTM"]+ ". Ignored.")
            del clean_row["PTM"]
    try:
        if "Mass" in clean_row:
            clean_row["Mass"]=utils.floating(clean_row["Mass"])
    except ValueError:
        message.warning("File "+file+", line "+str(index)+": wrong mass, "+clean_row["Mass"]+ ". Ignored.")
        del clean_row["Mass"]
    try:
        if "Hel" in clean_row:
            clean_row["Hel"]=integer(clean_row["Hel"])
    except ValueError:
        message.warning("File "+file+", line "+str(index)+": wrong helical position, "+clean_row["Hel"]+ ". Ignored.")
        del clean_row["Hel"]
    try:
        if "Length" in clean_row:
            clean_row["Length"]=integer(clean_row["Length"])
    except ValueError:
        message.warning("File "+file+", line "+str(index)+": wrong peptide length, "+str(clean_row["Length"])+ ". Ignored.")
        del clean_row["Length"]
    try:
        if "Begin" in clean_row:
            clean_row["Begin"]=integer(clean_row["Begin"])
    except ValueError:
        message.warning("File "+file+", line "+str(index)+": wrong begin position, "+clean_row["Begin"]+ ". Ignored.")
        del clean_row["Begin"]
    try:
        if "End" in clean_row:
            clean_row["End"]=integer(clean_row["End"])
    except ValueError:
        message.warning("File "+file+", line "+str(index)+": wrong end position, "+clean_row["End"]+ ". Ignored.")
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
    # Read the TSV file
    df = pd.read_csv(peptide_table_file_name, sep="\t")
    df = df.astype(str)
    # Drop empty columns
    df = df.drop(columns=[col for col in df.columns
    if 'Unnamed' in col and df[col].isna().all()
    ]) # isna or isnull ?
    list_of_headers= list(map(rename_field,list(map(utils.clean,df.columns.tolist()))))
    list_of_markers = df.to_dict(orient='records')
    set_of_markers=set()
    for index,m in enumerate(list_of_markers):
        new_markers=check_marker(m,index+1,peptide_table_file_name)
        set_of_markers.update(new_markers)
    if len(set_of_markers)==0:
        message.warning("File "+peptide_table_file_name+": no valid data found.")
    return set_of_markers, list_of_headers
     

def parse_peptide_tables(list_of_peptide_tables, list_of_constraints, taxonomy, warning_on=True):
    set_of_markers=set()
    for peptide_table  in list_of_peptide_tables:
        set_of_new_markers, list_of_headers=parse_peptide_table(peptide_table, warning_on)
        set_of_markers.update(set_of_new_markers)
    if list_of_constraints is not None and len(list_of_constraints)>0 :
        set_of_markers=lim.apply_limits(list_of_constraints, set_of_markers, taxonomy,True)
    if len(list_of_peptide_tables)>1: # TO DO
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
    
    
def build_peptide_table_from_set_of_markers(set_of_markers, outfile_name, sorted_headers, sorted_markers=None):
    set_of_headers={utils.restitute_field(key) for m in set_of_markers for key in m.field}
    sorted_headers=list(map(utils.restitute_field, sorted_headers))
    sorted_headers=utils.sort_headers(sorted_headers, set_of_headers)
    TSV_file = open(outfile_name, "w")
    writer = csv.DictWriter(TSV_file, fieldnames=sorted_headers, delimiter="\t")
    writer.writeheader()
    # ordering markers
    set_of_codes=[m.code() for m in set_of_markers]
    list_of_codes=utils.sort_headers(sorted_markers,set_of_codes)
    list_of_markers=list(set_of_markers)
    list_of_markers.sort(key=cmp_to_key(partial(marker_order, list_of_codes=list_of_codes)))
    for m in list_of_markers:
        dict={h:m.field[key] for key in m.field for  h in sorted_headers if utils.equiv(h,utils.restitute_field(key))}
        writer.writerow(dict)
    TSV_file.close()
    
    
