import re
from src import markers 
from src import sequences
from src import taxonomy as ta

def authorized_PTM(PTM_string, set_of_PTM):
    found_number=""
    for char in PTM_string:
        if char.isdigit():
            found_number += char
        else:
          if int(found_number)>0 and char not in set_of_PTM:
            return False
          found_number=""
    
    return True

def reduce(seq):
    #seq=seq.lower()
    return  seq.strip(' ')

def parse_sequence_info(line):
    dict={}
    SUFFIXES = [" OS", " OX", " Length", " SeqID", " GN", "PTM"]
    # Specific case: SeqID with no <SeqID=> flag
    re_seqid = re.compile('\S+\s')
    seqid=re_seqid.match(line)
    if seqid:
        if "=" not in seqid.group(0): 
            dict['SeqID']={seqid.group(0)}
    # Define a regular expression pattern with all possible fields
    pattern = re.compile(f'(?=(OS|OX|GN|PTM|SeqID|Length)=([^=]+))')
    # Find all matches in the line 
    matches = pattern.findall(line)
    for field, value in matches:
        new_value=next((value[: -len(suffix)] for suffix in SUFFIXES if value.endswith(suffix)), value)     
        dict[field]={reduce(v) for v in new_value.split(',')}
    return dict



    #for m in matches:
    #    print(m.group(1)+" : "+ m.group(2))
    

def apply_limits(limit_file,  set_of_items, taxonomy, is_marker):
    # is_marker: True marker, False sequence
    file=open(limit_file).read().splitlines()
    list_of_constraints=[parse_sequence_info(line) for line in file  if len(parse_sequence_info(line))>0]
    for i in range(len(list_of_constraints)):
        dict=list_of_constraints[i]
        for key in dict:
            print(str(key)+"->"+str(dict[key]), end=" ")
        print("")

    set_of_new_items=set()

    for i in range(len(list_of_constraints)):
        dict=list_of_constraints[i]
        current_set=set_of_items
        if "OS" in dict:
            if taxonomy:
                current_set={s for s in current_set for t in dict["OS"] if s.taxid in taxonomy.descendants[ta.search_taxid_from_taxon_name(t,taxonomy)]}
            else:
                current_set={s for s in current_set for t in dict["OS"] if reduce(s.taxon_name)==reduce(t)}
        if "OX" in dict:
            if taxonomy:
                print("dict")
                print(dict["OX"])
                current_set={s for s in current_set for t in dict["OX"] if s.taxid in taxonomy.descendants[t]}
            else:
                current_set={s for s in current_set for t in dict["OX"] if s.taxid==reduce(t)}
        if "GN" in dict:
            current_set={s for s in current_set if s.protein in dict["GN"]}
        if "SeqID" in dict:
            current_set={s for s in current_set if s.seqid in dict["SeqID"]}
        if is_marker and "PTM" in dict:                     
            current_set={s for s in current_set if authorized_PTM(s.ptm,dict["PTM"])}
        print("----")
        for m in current_set:
            print(str(m))
        set_of_new_items.update(current_set)

    return set_of_new_items
