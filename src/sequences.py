"""

sequences.py                         


"""
from pyteomics import parser

from src import markers, collagen, supplement, utils

class Sequence(object):
    def __init__(self, field=None):
        if field is None:
            self.field={}
        else:
            self.field = field

    def __len__(self):
        return len(self.field["Sequence"])
        
    def __str__(self):
        return str(self.field)

    def seqid(self):
        return self.field.get("SeqID")
        
    def taxid(self):
        return self.field.get("OX")
        
    def taxon_name(self):
        return self.field.get("OS")
        
    def protein(self):
        return self.field.get("GN")
        
    def sequence(self):
        return self.field.get("Sequence")


def raw_in_silico_digestion(seq, config_digestion, min_length=None, max_length=None):
    """ build a set of peptides from a sequence by in silico digestion"""
    if min_length is None:
        min_length=config_digestion["min_peptide_length"]
    if max_length is None:
        max_length=config_digestion["max_peptide_length"]
    enzyme=config_digestion["enzyme"]
    number_of_misscleavages=config_digestion["number_of_missed_cleavages"]
    set_of_peptides=parser.icleave(seq, parser.expasy_rules[enzyme], number_of_misscleavages, min_length, max_length)
    return {y for (x,y) in set_of_peptides}


def in_silico_digestion(set_of_sequences, config_digestion, min_length=None, max_length=None, mature=True):
    """ build a set of markers from a set of sequences by in silico digestion"""
    if min_length is None:
        min_length=config_digestion["min_peptide_length"]
    if max_length is None:
        max_length=config_digestion["max_peptide_length"]
    enzyme=config_digestion["enzyme"]
    number_of_misscleavages=config_digestion["number_of_missed_cleavages"]
    set_of_markers=set()
    for s in set_of_sequences:
        min_pos, max_pos=None, None
        if mature:
            (min_pos,max_pos)=collagen.helical_region(s)
        if min_pos is None:
            helical=False
            (min_pos,max_pos)=(0,len(s))
        else:
            helical = True
            (min_pos,max_pos)=(min_pos-1, max_pos-1)
        mature_seq=s.sequence()[min_pos:max_pos]
        not_collagen=False
        for j in range(0,len(mature_seq),3):
            if mature_seq[j] not in {'G','X'}: ## to change
                not_collagen=True
        set_of_peptides=  parser.icleave(mature_seq, parser.expasy_rules[enzyme], number_of_misscleavages, min_length, max_length)
        for (pos, peptide) in set_of_peptides:
            if not utils.is_aa_sequence(peptide) or (helical and not collagen.is_collagen_peptide(peptide)):
                continue
            dict_marker = {
                "Sequence": peptide,
                "OX": s.taxid(),
                "OS": s.taxon_name(),
                "GN": s.protein(),
                "Rank": "species",
                "SeqID": s.seqid(),
                "Length": len(peptide),
                "Begin": min_pos + pos + 1,
                "End": min_pos + pos + len(peptide),
                "Status": "Genetics",
                "Digestion": "Yes",
            }
            if helical:
                dict_marker["Hel"] = pos + 1
            new_marker=markers.Marker(field=dict_marker)
            set_of_markers.add(new_marker)
    supplement.add_marker_names(set_of_markers)
    return set_of_markers

def is_digested_peptide(peptide, config_digestion):
    #enzyme = config_digestion["enzyme"]
    #print(parser.expasy_rules[enzyme])
    return peptide[0]!='P'






