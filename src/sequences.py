"""

sequences.py                         


"""

from pyteomics import parser
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from utils import *
from markers import *
from peptide_table import *

class Sequence(object):
    def __init__(self, sequence="", taxid="", taxon_name="", seqid="", protein="", rank="", comment=""):
        self.sequence = sequence
        self.taxid = taxid
        self.taxon_name = taxon_name
        self.seqid=seqid
        self.protein=protein
        self.rank = rank
        self.comment = comment

    def __len__(self):
        return (len(self.sequence))
        
    def __str__(self):
        return f"Sequence{self.seqid} \t{self.taxid}\t{self.taxon_name.lstrip()}\t{self.sequence}\t{self.protein}\t{self.rank}\t{self.comment}"


def in_silico_digestion(set_of_sequences, number_of_misscleavages, min_length, mature):
    """ build a set of markers from a set of sequences by in silico digestion"""
    set_of_markers=set()
    for s in set_of_sequences:
        if (mature):
            (min,max)=mature_sequence(s)
        else:
            (min,max)=(0, len(s))
        mature_seq=s.sequence[min:max]
        set_of_peptides=parser.icleave(mature_seq, parser.expasy_rules['trypsin'], number_of_misscleavages, min_length)
        for (pos, peptide) in set_of_peptides:
            if ('Z' in peptide or 'B' in peptide or 'X' in peptide):
                continue
            new_marker=Marker()
            new_marker.sequence=peptide
            new_marker.taxid=s.taxid
            new_marker.taxon_name=s.taxon_name
            new_marker.protein=s.protein
            new_marker.seqid=s.seqid
            new_marker.begin=min+pos
            new_marker.end=new_marker.begin+len(peptide)-1
            new_marker.rank="species"
            new_marker.comment="in silico digestion"
            if mature:
                new_marker.comment=new_marker.comment + " - mature "
            set_of_markers.add(new_marker)
    return set_of_markers
      
def mature_sequence(seq):
    human_seq_col1A1_start="QLSYGYDEKSTGGISVPGPMGPSGPRGLPG"
    human_seq_col1A1_end="PGPPSAGFDFSFLPQPPQEKAHDGGRYYRA"
    human_seq_col1A2_start="QYDGKGVGLGPGPMGLMGPRGPPGAAGAPG"
    human_seq_col1A2_end="GPPGPPGPPGPPGVSGGGYDFGYDGDFYRA"

    min=0
    max=len(seq.sequence)
    
    if seq.protein=="COL1A1":
        al1=pairwise2.align.localms(human_seq_col1A1_start, seq.sequence,1,-1,-1,-1, one_alignment_only=True)
        al2=pairwise2.align.localms(human_seq_col1A1_end, seq.sequence,1,-1,-1,-1, one_alignment_only=True)
        if al1[0].score>20.0:
            min=al1[0].start
        if al2[0].score>20.0:
            max=al2[0].end
    if seq.protein=="COL1A2":
        al1=pairwise2.align.localms(human_seq_col1A2_start, seq.sequence,1,-1,-1,-1, one_alignment_only=True)
        al2=pairwise2.align.localms(human_seq_col1A2_end, seq.sequence,1,-1,-1,-1, one_alignment_only=True)
        if al1[0].score>17.0:
            min=al1[0].start
        if al2[0].score>17.0:
            max=al2[0].end

    return (min,max)
        
