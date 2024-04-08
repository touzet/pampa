import re
from src import markers 
from src import sequences
from src import taxonomy as ta
from src import message
from src import utils 



def reduce(seq):
    #seq=seq.lower()
    return  seq.strip(' ')

def parse_sequence_info(line, limit):
    dict={}
    if len(line)==0:
        return dict
    # Specific case: SeqID with no "SeqID=" flag
    if "=" not in line:
        dict['FileName']=line
        return dict
    # Define a regular expression pattern with all possible fields
    SUFFIXES = [" OS", " OX", " Length", " SeqID", " GN", " PTM"]
    pattern = re.compile(f'(?=(OS|OX|GN|PTM|SeqID|Length)=([^=]+))')
    # Find all matches in the line 
    matches = pattern.findall(line)
    if len(matches)==0:
        message.warning("File "+limit+", line "+line+" is incorrect. Ignored.")
    for field, value in matches:
        new_value=next((value[: -len(suffix)] for suffix in SUFFIXES if value.endswith(suffix)), value)     
        dict[field]={reduce(v) for v in new_value.split(',')}
    return dict


def parse_limits(limit):
    limit_file=open(limit).read().splitlines()
    list_of_constraints=[constraints for constraints in [parse_sequence_info(line,limit) for line in limit_file]  if len(constraints)>0]
    return list_of_constraints

def apply_limits(list_of_constraints,  set_of_items, taxonomy, is_marker):
    
    set_of_new_items=set()

    for i in range(len(list_of_constraints)): # why i ?
        dict=list_of_constraints[i]
        current_set=set_of_items
        if "OS" in dict:
            if taxonomy:
                current_set={s for s in current_set for t in dict["OS"] if s.taxid in taxonomy.descendants[ta.search_taxid_from_taxon_name(t,taxonomy)]}
            else:
                current_set={s for s in current_set for t in dict["OS"] if reduce(t) in reduce(s.taxon_name)}
        if "OX" in dict:
            if taxonomy:
                current_set={s for s in current_set for t in dict["OX"] if s.taxid in taxonomy.descendants[t]}
            else:
                current_set={s for s in current_set for t in dict["OX"] if s.taxid==reduce(t)}
        if "GN" in dict:
            current_set={s for s in current_set if s.protein in dict["GN"]}
        if "SeqID" in dict:
            current_set={s for s in current_set if s.seqid in dict["SeqID"]}
        if is_marker and "PTM" in dict:                     
            current_set={s for s in current_set if utils.is_PTM(s.ptm,dict["PTM"])}
        
        set_of_new_items.update(current_set)

    return set_of_new_items
