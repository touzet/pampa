import csv
import shutil
from src import utils
from src import sequences




def search_for_incomplete_markers(set_of_markers, required_fields):
    set_of_incomplete_markers=set()
    set_of_complete_markers=set()
    set_of_incomplete_fields=set()
    for m in set_of_markers:
        if required_fields.issubset(m.field):
            set_of_complete_markers.add(m)
        else:
            set_of_incomplete_markers.add(m)
            for f in required_fields:
                if f not in m.field:
                    set_of_incomplete_fields.add(f)
            
    return set_of_incomplete_markers, set_of_complete_markers, set_of_incomplete_fields

def add_digestion_status(set_of_markers, set_of_sequences):
    set_of_seqid={m.seqid() for m in set_of_markers if "Digestion" not in m.field}
    for seqid in set_of_seqid:
        compatible_sequences={s.sequence() for s in set_of_sequences if s.seqid()==seqid}
        if len(compatible_sequences)==0:
            continue
        seq=next(iter(compatible_sequences))
        set_of_peptides=sequences.raw_in_silico_digestion(seq)
        for m in set_of_markers:
            if m.seqid()==seqid:
                if m.sequence() in set_of_peptides:
                    m.field["Digestion"]="Yes"
                else:
                    m.field["Digestion"]="No"
    return set_of_markers
        
def add_length(set_of_markers):
    for m in set_of_markers:
        if m.length() is None and m.sequence() is not None:
            m.field["Length"]=len(m.sequence())
    return set_of_markers
