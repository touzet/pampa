import re
from src import taxonomy as ta
from src import message
from src import utils 


def parse_sequence_info(line, limit):
    constraint_dict={}
    if len(line)==0:
        return constraint_dict
    # Specific case: SeqID with no "SeqID=" flag
    if "=" not in line:
        constraint_dict['FileName']=line
        return constraint_dict
    # Split the text into parts by spaces, while respecting '=' as the delimiter
    pattern = r'(\w+)\s*=\s*(.+?)(?=\s+\w+\s*=|$)'
    matches = re.findall(pattern, line)    # Create a dictionary by splitting each pair at '='
    if len(matches)==0:
        message.warning("File "+limit+", line "+line+" is incorrect. Ignored.")
    constraint_dict = {field.strip(): {utils.standard(v) for v in sentence.split(',')}  for field, sentence in matches}
    return constraint_dict

def parse_limits(limit):
    if limit is None:
        return []
    limit_file=open(limit).read().splitlines()
    list_of_constraints=[constraints for constraints in [parse_sequence_info(line,limit) for line in limit_file]  if len(constraints)>0]
    return list_of_constraints

#items can be sequences or peptide markers
def apply_limits(list_of_constraints,  set_of_items, taxonomy, is_marker):
    set_of_item_fields={key for s in set_of_items for key in s.field}
    global_constraint=False
    set_of_new_items=set()
    for constraint_dict in list_of_constraints:
        constraint=False
        current_set=set_of_items
        if "OS" in constraint_dict:
            constraint=True
            if taxonomy:
                current_set={s for s in current_set for t in constraint_dict["OS"] if s.taxid() in taxonomy.descendants[ta.search_taxid_from_taxon_name(t,taxonomy)]}
            else:
                current_set={s for s in current_set for t in constraint_dict["OS"] if utils.equiv(t, s.taxon_name())}
        if "OX" in constraint_dict:
            constraint=True
            if taxonomy:
                current_set={s for s in current_set for t in constraint_dict["OX"] if s.taxid() in taxonomy.descendants[t]}
            else:
                current_set={s for s in current_set for t in constraint_dict["OX"] if utils.equiv(t, s.taxid())}
        if is_marker and "PTM" in constraint_dict:
            constraint=True
            current_set={s for s in current_set if utils.is_PTM(s.PTM(),constraint_dict["PTM"])}
        for field in constraint_dict.keys()-{"OS","OX","PTM"}:
            if field not in set_of_item_fields:
                continue
            constraint=True
            current_set={s for s in current_set if field in s.field and s.field[field] in constraint_dict[field]}
            
        if constraint:
            set_of_new_items.update(current_set)
            global_constraint=True
            
    if global_constraint:
        return set_of_new_items
    else:
        return set_of_items

def extract_taxonomical_constraints(list_of_constraints):
    constraint=set()
    for d in list_of_constraints:
        if "Taxonomy" in d :
            constraint.update(d["Taxonomy"])
    if len(constraint)==0:
        return None
    else:
        return constraint

def deamidated_codes(list_of_constraints, deamidation):
    set_of_codes=set()
    if not deamidation:
        return set_of_codes
    for constraint in list_of_constraints:
        if 'Deamidation' in constraint:
         set_of_codes.update(constraint['Deamidation'])
    if len(set_of_codes)>0:
        return set_of_codes
    else:
        return None