"""

sequences.py                         


"""
import re
from pyteomics import parser
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

from src import markers 


class Sequence(object):
    def __init__(self, sequence="", taxid="", taxon_name="", seqid="", helical=None, protein="", rank="", comment=""):
        self.sequence = sequence
        self.taxid = taxid
        self.taxon_name = taxon_name
        self.seqid=seqid
        self.helical=helical
        self.protein=protein
        self.rank = rank
        self.comment = comment

    def __len__(self):
        return (len(self.sequence))
        
    def __str__(self):
        return f"Sequence{self.seqid} \t{self.taxid}\t{self.taxon_name.lstrip()}\t{self.sequence}\t{self.protein}\t{self.rank}\t{self.comment}"


def helical_region(seq):
    # positions are 1-based
    pattern=re.compile('(G\w{2}){5,}')
    matches=re.finditer(pattern, seq)
    positions=[match.span() for match in matches]
    if (len(positions)==0):
          return None, None
    start_match=positions[0][0]
    end_match=positions[0][1]
    for segment in positions[1:]:
        if segment[0]-end_match>6:# allowing for errors in the aa sequence 
            start_match= segment[0]
        end_match=segment[1]
    pos_GPM=seq[start_match:].find("GPM")
    if pos_GPM==-1 or pos_GPM>6:
        return start_match+1, end_match+1  
    else:
        return pos_GPM+start_match +1, end_match+1 
    

def in_silico_digestion(set_of_sequences, number_of_misscleavages=1, min_length=12, max_length=33, mature=True):
    """ build a set of markers from a set of sequences by in silico digestion"""
    set_of_markers=set()
    for s in set_of_sequences:
        min=None
        if mature:
            (min,max)=helical_region(s.sequence)
        if min==None:
            helical=False
            (min,max)=(0,len(s.sequence))
        else:
            (min,max)=(min-1, max-1)
            helical=True               
        mature_seq=s.sequence[min:max]
       
        set_of_peptides=parser.icleave(mature_seq, parser.expasy_rules['trypsin'], number_of_misscleavages, min_length)
        for (pos, peptide) in set_of_peptides:
            if ('Z' in peptide or 'B' in peptide or 'X' in peptide): ## amino acids
                continue
            if len(peptide)>max_length:
                continue
            new_marker=markers.Marker()
            new_marker.sequence=peptide
            new_marker.taxid=s.taxid
            new_marker.taxon_name=s.taxon_name
            new_marker.protein=s.protein
            if not helical :
                new_marker.helical=None
            else:
                new_marker.helical=str(pos+1)
            new_marker.seqid=s.seqid
            new_marker.begin=min+pos+1
            new_marker.end=new_marker.begin+len(peptide)-1
            new_marker.rank="species"
            new_marker.comment="in silico digestion"
            if mature:
                new_marker.comment=new_marker.comment + " - mature "
            set_of_markers.add(new_marker)
    return set_of_markers


def alignment(target, query):
    return pairwise2.align.localms(target,query,1,-1,-1,-1, one_alignment_only=True)

def mature_sequence(seq):
    HUMAN_COL1A1_START="QLSYGYDEKSTGGISVPGPMGPSGPRGLPG"
    CHICKEN_COL1A1_START="QMSYGYDEKSAGVAVPGPMGPAGPRGLPG"
    HUMAN_COL1A1_END="PGPPSAGFDFSFLPQPPQEKAHDGGRYYRA"
    CHICKEN_COL1A1_END="PGPPSGGFDLSFLPQPPQEKAHDGGRYYRA"
    HUMAN_COL1A2_START="QYDGKGVGLGPGPMGLMGPRGPPGAAGAPG"
    CHICKEN_COL1A2_START="QYDPSKAADFGPGPMGLMGPRGPPGASGPPG"
    HUMAN_COL1A2_END="GPPGPPGPPGPPGVSGGGYDFGYDGDFYRA"
    CHICKEN_COL1A2_END="GPPGPPGPPGPPGPNGGGYEVGFDAEYYR"
    
    HUMAN_COL1A1_MATURE="QLSYGYDEKSTGGISVPGPMGPSGPRGLPGPPGAPGPQGFQGPPGEPGEPGASGPMGPRGPPGPPGKNGDDGEAGKPGRPGERGPPGPQGARGLPGTAGLPGMKGHRGFSGLDGAKGDAGPAGPKGEPGSPGENGAPGQMGPRGLPGERGRPGAPGPAGARGNDGATGAAGPPGPTGPAGPPGFPGAVGAKGEAGPQGPRGSEGPQGVRGEPGPPGPAGAAGPAGNPGADGQPGAKGANGAPGIAGAPGFPGARGPSGPQGPGGPPGPKGNSGEPGAPGSKGDTGAKGEPGPVGVQGPPGPAGEEGKRGARGEPGPTGLPGPPGERGGPGSRGFPGADGVAGPKGPAGERGSPGPAGPKGSPGEAGRPGEAGLPGAKGLTGSPGSPGPDGKTGPPGPAGQDGRPGPPGPPGARGQAGVMGFPGPKGAAGEPGKAGERGVPGPPGAVGPAGKDGEAGAQGPPGPAGPAGERGEQGPAGSPGFQGLPGPAGPPGEAGKPGEQGVPGDLGAPGPSGARGERGFPGERGVQGPPGPAGPRGANGAPGNDGAKGDAGAPGAPGSQGAPGLQGMPGERGAAGLPGPKGDRGDAGPKGADGSPGKDGVRGLTGPIGPPGPAGAPGDKGESGPSGPAGPTGARGAPGDRGEPGPPGPAGFAGPPGADGQPGAKGEPGDAGAKGDAGPPGPAGPAGPPGPIGNVGAPGAKGARGSAGPPGATGFPGAAGRVGPPGPSGNAGPPGPPGPAGKEGGKGPRGETGPAGRPGEVGPPGPPGPAGEKGSPGADGPAGAPGTPGPQGIAGQRGVVGLPGQRGERGFPGLPGPSGEPGKQGPSGASGERGPPGPMGPPGLAGPPGESGREGAPGAEGSPGRDGSPGAKGDRGETGPAGPPGAPGAPGAPGPVGPAGKSGDRGETGPAGPAGPVGPVGARGPAGPQGPRGDKGETGEQGDRGIKGHRGFSGLQGPPGPPGSPGEQGPSGASGPAGPRGPPGSAGAPGKDGLNGLPGPIGPPGPRGRTGDAGPVGPPGPPGPPGPPGPPSAGFDFSFLPQPPQEKAHDGGRYYRA"

    HUMAN_COL1A2_MATURE="QYDGKGVGLGPGPMGLMGPRGPPGAAGAPGPQGFQGPAGEPGEPGQTGPAGARGPAGPPGKAGEDGHPGKPGRPGERGVVGPQGARGFPGTPGLPGFKGIRGHNGLDGLKGQPGAPGVKGEPGAPGENGTPGQTGARGLPGERGRVGAPGPAGARGSDGSVGPVGPAGPIGSAGPPGFPGAPGPKGEIGAVGNAGPAGPAGPRGEVGLPGLSGPVGPPGNPGANGLTGAKGAAGLPGVAGAPGLPGPRGIPGPVGAAGATGARGLVGEPGPAGSKGESGNKGEPGSAGPQGPPGPSGEEGKRGPNGEAGSAGPPGPPGLRGSPGSRGLPGADGRAGVMGPPGSRGASGPAGVRGPNGDAGRPGEPGLMGPRGLPGSPGNIGPAGKEGPVGLPGIDGRPGPIGPAGARGEPGNIGFPGPKGPTGDPGKNGDKGHAGLAGARGAPGPDGNNGAQGPPGPQGVQGGKGEQGPPGPPGFQGLPGPSGPAGEVGKPGERGLHGEFGLPGPAGPRGERGPPGESGAAGPTGPIGSRGPSGPPGPDGNKGEPGVVGAVGTAGPSGPSGLPGERGAAGIPGGKGEKGEPGLRGEIGNPGRDGARGAPGAVGAPGPAGATGDRGEAGAAGPAGPAGPRGSPGERGEVGPAGPNGFAGPAGAAGQPGAKGERGAKGPKGENGVVGPTGPVGAAGPAGPNGPPGPAGSRGDGGPPGMTGFPGAAGRTGPPGPSGISGPPGPPGPAGKEGLRGPRGDQGPVGRTGEVGAVGPPGFAGEKGPSGEAGTAGPPGTPGPQGLLGAPGILGLPGSRGERGLPGVAGAVGEPGPLGIAGPPGARGPPGAVGSPGVNGAPGEAGRDGNPGNDGPPGRDGQPGHKGERGYPGNIGPVGAAGAPGPHGPVGPAGKHGNRGETGPSGPVGPAGAVGPRGPSGPQGIRGDKGEPGEKGPRGLPGLKGHNGLQGLPGIAGHHGDQGAPGSVGPAGPRGPAGPSGPAGKDGRTGHPGTVGPAGIRGPQGHQGPAGPPGPPGPPGPPGVSGGGYDFGYDGDFYRA"
    
    CHICKEN_COL1A1_MATURE="QMSYGYDEKSAGVAVPGPMGPAGPRGLPGPPGAPGPQGFQGPPGEPGEPGASGPMGPRGPAGPPGKNGDDGEAGKPGRPGQRGPPGPQGARGLPGTAGLPGMKGHRGFSGLDGAKGQPGPAGPKGEPGSPGENGAPGQMGPRGLPGERGRPGPSGPAGARGNDGAPGAAGPPGPTGPAGPPGFPGAAGAKGETGPQGARGSEGPQGSRGEPGPPGPAGAAGPAGNPGADGQPGAKGATGAPGIAGAPGFPGARGPSGPQGPSGAPGPKGNSGEPGAPGNKGDTGAKGEPGPAGVQGPPGPAGEEGKRGARGEPGPAGLPGPAGERGAPGSRGFPGADGIAGPKGPPGERGSPGAVGPKGSPGEAGRPGEAGLPGAKGLTGSPGSPGPDGKTGPPGPAGQDGRPGPAGPPGARGQAGVMGFPGPKGAAGEPGKPGERGAPGPPGAVGAAGKDGEAGAQGPPGPTGPAGERGEQGPAGAPGFQGLPGPAGPPGEAGKPGEQGVPGNAGAPGPAGARGERGFPGERGVQGPPGPQGPRGANGAPGNDGAKGDAGAPGAPGNEGPPGLEGMPGERGAAGLPGAKGDRGDPGPKGADGAPGKDGLRGLTGPIGPPGPAGAPGDKGEAGPPGPAGPTGARGAPGDRGEPGPPGPAGFAGPPGADGQPGAKGETGDAGAKGDAGPPGPAGPTGAPGPAGZVGAPGPKGARGSAGPPGATGFPGAAGRVGPPGPSGNIGLPGPPGPAGKZGSKGPRGETGPAGRPGEPGPAGPPGPPGEKGSPGADGPIGAPGTPGPQGIAGQRGVVGLPGQRGERGFPGLPGPSGEPGKQGPSGASGERGPPGPMGPPGLAGPPGEAGREGAPGAEGAPGRDGAAGPKGDRGETGPAGPPGAPGAPGAPGPVGPAGKNGDRGETGPAGPAGPPGPAGARGPAGPQGPRGDKGETGEQGDRGMKGHRGFSGLQGPPGPPGAPGEQGPSGASGPAGPRGPPGSAGAAGKDGLNGLPGPIGPPGPRGRTGEVGPVGPPGPPGPPGPPGPPSGGFDLSFLPQPPQEKAHDGGRYYRA"

    CHICKEN_COL1A2_MATURE="QYDPSKAADFGPGPMGLMGPRGPPGASGPPGPPGFQGVPGEPGEPGQTGPQGPRGPPGPPGKAGEDGHPGKPGRPGERGVAGPQGARGFPGTPGLPGFKGIRGHNGLDGQKGQPGTPGTKGEPGAPGENGTPGQPGARGLPGERGRIGAPGPAGARGSDGSAGPTGPAGPIGAAGPPGFPGAPGAKGEIGPAGNVGPTGPAGPRGEIGLPGSSGPVGPPGNPGANGLPGAKGAAGLPGVAGAPGLPGPRGIPGPPGPAGPSGARGLVGEPGPAGAKGESGNKGEPGAAGPPGPPGPSGEEGKRGSNGEPGSAGPPGPAGLRGVPGSRGLPGADGRAGVMGPAGNRGASGPVGAKGPNGDAGRPGEPGLMGPRGLPGQPGSPGPAGKEGPVGFPGADGRVGPIGPAGNRGEPGNIGFPGPKGPTGEPGKPGEKGNVGLAGPRGAPGPEGNNGAQGPPGVTGNQGAKGETGPAGPPGFQGLPGPSGPAGEAGKPGERGLHGEFGVPGPAGPRGERGLPGESGAVGPAGPIGSRGPSGPPGPDGNKGEPGNVGPAGAPGPAGPGGIPGERGVAGVPGGKGEKGAPGLRGDTGATGRDGARGLPGAIGAPGPAGGAGDRGEGGPAGPAGPAGARGIPGERGEPGPVGPSGFAGPPGAAGQPGAKGERGPKGPKGETGPTGAIGPIGASGPPGPVGAAGPAGPRGDAGPPGMTGFPGAAGRVGPPGPAGITGPPGPPGPAGKDGPRGLRGDVGPVGRTGEQGIAGPPGFAGEKGPSGEAGAAGPPGTPGPQGILGAPGILGLPGSRGERGLPGIAGATGEPGPLGVSGPPGARGPSGPVGSPGPNGAPGEAGRDGNPGNDGPPGRDGAPGFKGERGAPGNPGPSGALGAPGPHGQVGPSGKPGNRGDPGPVGPVGPAGAFGPRGLAGPQGPRGEKGEPGDKGHRGLPGLKGHNGLQGLPGLAGQHGDQGPPGNNGPAGPRGPPGPSGPPGKDGRNGLPGPIGPAGVRGSHGSQGPAGPPGPPGPPGPPGPNGGGYEVGFDAEYYR"
    
    #COL1A1=[(HUMAN_COL1A1_START, HUMAN_COL1A1_END), (CHICKEN_COL1A1_START, CHICKEN_COL1A1_END)]
    COL1A2=[HUMAN_COL1A2_MATURE, CHICKEN_COL1A2_MATURE]
    COL1A1=[HUMAN_COL1A1_MATURE, CHICKEN_COL1A1_MATURE]
    
    start_query=0
    end_query=len(seq.sequence)

    if seq.protein=="COL1A1":
        list_of_target_sequences=COL1A1
    elif seq.protein=="COL1A2":
        list_of_target_sequences=COL1A2
    else:
        return (start_query,end_query)
    
    max_length=50
    for target in list_of_target_sequences:
        al=alignment(target, seq.sequence)
        if al[0].end-al[0].start>max_length:
            max_length= al[0].end-al[0].start
            end_query=al[0].end
            start_query=al[0].start 
   
    return (start_query,end_query)
        


    
def reduce(seq):
    seq=seq.replace(" ", "")
    seq=seq.lower()
    return  seq


# still useful ? 
def limit_sequences(set_of_sequences, limit_file):
    """
    build a new set of sequences that is compliant with the specification of the limit_file (TSV file, -l option) 
    """
    new_set=set_of_markers
    summary_file=open(limit_file).read().splitlines()
    for line in summary_file:      
        if line[:3]=="GN=":
            pattern = re.compile(r'GN=([\w\s,]+)')
            list_of_matches = pattern.findall(line)
            matches= list_of_matches[0].split(',')
            new_set={s for s in new_set if s.protein in matches}

        if line[:3]=="OS=":
            pattern = re.compile(r'OS=([\w\s,]+)')
            list_of_matches = pattern.findall(line)
            matches= list_of_matches[0].split(',')
            matches=[reduce(match) for match in matches]
            new_set={s for s in new_set if reduce(s.taxon_name) in matches}   

        if line[:6]=="SeqID=":
            pattern = re.compile(r'SeqID=([\w\s,]+)')
            list_of_matches = pattern.findall(line)
            matches= list_of_matches[0].split(',')
            matches=[reduce(match) for match in matches]
            new_set={s for s in new_set if reduce(s.seqid) in matches} 

    summary_file.close()
    return new_set

  
    
