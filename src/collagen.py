import re
from collections import Counter
import itertools

from src import message, compute_masses

def helical_region(seq):
    """
    input: COLLAGENE sequence
    output: (start position, end position) of helical region
    """
    # positions are 1-based
    pattern=re.compile('(G\w{2}){5,}')
    matches=re.finditer(pattern, seq.sequence())
    positions=[match.span() for match in matches]
    if (len(positions)==0):
          return None, None
    start_match=positions[0][0]
    end_match=positions[0][1]
    for segment in positions[1:]:
        if segment[0]-end_match>4:# allowing for errors in the aa sequence
            start_match= segment[0]
        end_match=segment[1]
    pos_GPM=seq.sequence()[start_match:].find("GPM")
    if pos_GPM==-1 :
        return start_match+1, end_match+1
    elif pos_GPM>6:
        return start_match+1, end_match+1
    else:
        return pos_GPM+start_match +1, end_match+1


def mature(sequence):
    min, max = helical_region(sequence)
    if min >1 or max < len(sequence.sequence()):
        sequence.field["Sequence"] = sequence.sequence()[min - 1:max - 1]
        sequence.field["Mature"] = min, max
    else:
        sequence.field["Mature"] = None
    return sequence

def mature_sequences(set_of_sequences):
    return {mature(seq) for seq in set_of_sequences}

def prohibited_trimers():
    return {'GAC', 'GAW', 'GAY', 'GCA', 'GCC', 'GCD', 'GCE', 'GCF', 'GCG', 'GCH', 'GCI', 'GCK', 'GCL', 'GCM', 'GCN', 'GCP', 'GCQ', 'GCR', 'GCS', 'GCT', 'GCV', 'GCW', 'GCY', 'GDC', 'GDE', 'GDH', 'GDM', 'GDN', 'GDW', 'GDY', 'GEC', 'GEW', 'GEY', 'GFC', 'GFD', 'GFE', 'GFF', 'GFG', 'GFI', 'GFM', 'GFR', 'GFW', 'GFY', 'GGC', 'GGF', 'GGH', 'GGI', 'GGW', 'GGY', 'GHC', 'GHD', 'GHE', 'GHF', 'GHI', 'GHL', 'GHM', 'GHT', 'GHW', 'GHY', 'GIC', 'GIF', 'GIH', 'GII', 'GIW', 'GIY', 'GKC', 'GKF', 'GKK', 'GKL', 'GKM', 'GKW', 'GLC', 'GLE', 'GLF', 'GLW', 'GLY', 'GMC', 'GME', 'GMF', 'GMG', 'GMH', 'GMI', 'GML', 'GMM', 'GMQ', 'GMV', 'GMW', 'GMY', 'GNC', 'GNE', 'GNF', 'GNG', 'GNH', 'GNM', 'GNW', 'GNY', 'GPC', 'GPW', 'GPY', 'GQC', 'GQE', 'GQF', 'GQG', 'GQW', 'GQY', 'GRC', 'GRF', 'GRH', 'GRK', 'GRL', 'GRM', 'GRQ', 'GRR', 'GRW', 'GRY', 'GSC', 'GSF', 'GSW', 'GSY', 'GTC', 'GTE', 'GTF', 'GTG', 'GTI', 'GTM', 'GTW', 'GTY', 'GVC', 'GVE', 'GVG', 'GVH', 'GVW', 'GWA', 'GWC', 'GWD', 'GWE', 'GWF', 'GWG', 'GWH', 'GWI', 'GWK', 'GWL', 'GWM', 'GWN', 'GWP', 'GWQ', 'GWR', 'GWS', 'GWT', 'GWV', 'GWW', 'GWY', 'GYC', 'GYD', 'GYE', 'GYF', 'GYG', 'GYH', 'GYI', 'GYK', 'GYL', 'GYM', 'GYQ', 'GYR', 'GYT', 'GYV', 'GYW', 'GYY'}

def check_GXY_pattern(sequence):
    seq=sequence.sequence()
    match = re.search(r'G.{2}G.{2}G', seq)
    phase = match.start() if match else 3
    if phase > 2:
        message.warning(sequence.seqid()+" : Not a collagen.")
        return False
    position=phase
    while position <len(seq)  and seq[position] in {'G','X'}:
        position += 3
    if position < len(seq):
        message.warning(
                    sequence.seqid() + " : error in the GXY pattern at position " + str(position) + ". Sequence removed." + seq[position - 5:position + 6])
        return False
    else:
        return True


def check_and_correct_GXY_pattern(sequence):
    seq=sequence.sequence()
    match = re.search(r'G.{2}G.{2}G', seq)
    phase = match.start() if match else 3
    if phase > 2:
        message.warning(sequence.seqid()+" : Not a collagen.")
        return sequence, -1
    position=phase
    while position<len(seq):
        while position <len(seq)  and seq[position] in {'G','X'}:
            position += 3
        if position>=len(seq):
            continue
        if seq[position] not in {'G','X'}:
            match = re.search(r'G.{2}G.{2}G', seq[position-3:])
            new_phase = match.start() if match else 7
            if new_phase < 7:
                message.warning(
                    sequence.seqid() + " : edition of the GXY pattern at position " + str(position) + ", " + seq[position - 5:position + 6])
            else:
                message.warning(sequence.seqid() + " : end of GXY pattern at position " + str(position-1)+", "+seq[position-3:position+4])
                seq=seq[:position]
            if new_phase == 1:
                seq=seq[:position-3]+"XX"+ seq[position-2:]
            if new_phase == 2:  # deletion
                seq = seq[:position-2] + 'X' + seq[position-1:]
            if new_phase == 4:
                seq = seq[:position - 2] + 'XX' + seq[position+1:]
            if new_phase == 5:
                seq = seq[:position-2] + 'XX' + seq[position+2 :]
            if new_phase == 6 :
                seq=seq[:position]+'X'+seq[position+1:]
            position+=3
    sequence.field["Sequence"]=seq
    return sequence, phase

def check_collagen_sequences(set_of_sequences):
    set_of_mature_sequences={mature(seq) for seq in set_of_sequences}
    set_of_mature_sequences={seq for seq in set_of_mature_sequences if check_GXY_pattern(seq)}
    return set_of_mature_sequences

def phasing_GXY_pattern(sequence):
    match = re.search(r'G.{2}G.{2}G', sequence)
    phase = match.start() if match else 3
    if phase > 2:
        message.warning("Not a collagen sequence: " + sequence)
        return -1
    for position in range(phase, len(sequence), 3):
        if sequence[position] != 'G':
            message.warning("Not a collagen sequence, disruption of the phase at position " +str(position) + "\n" +sequence)
            continue
    return phase

def count_trimers(set_of_sequences):
    # TO DO: phaser les séquences et vérifier la longueur
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    absent_trimers = {"G" + "".join(p) for p in itertools.product(amino_acids, repeat=2)}
    set_of_sequences= {check_GXY_pattern(sequence)[0] for sequence in set_of_sequences}
       # p=phasing_GXY_pattern(seq.sequence())
    for prot in ["COL1A1", "COL1A2"]:
        set_of_strings={s.sequence() for s in set_of_sequences if s.protein()==prot}
        seq_length=len(next(iter(set_of_strings)))
        trimers=Counter()
        X=Counter()
        Y=Counter()
        #p=phasing_GXY_pattern(s)
        for i in range(0, seq_length-2, 3):
            set_of_trimers={seq[i:i+3] for seq in set_of_strings if len(seq)>i+2}
            for trimer in set_of_trimers:
                trimers[trimer] += 1
                X[trimer[1]] += 1
                Y[trimer[2]] += 1
        total_trimers = sum(trimers.values())
        total_X = sum(X.values())
        total_Y = sum(Y.values())
        # Generate all 3-mers starting with G
        absent_trimers =  absent_trimers - trimers.keys()
        for w, count in sorted(trimers.items()):
            a=count/total_trimers
            b=X[w[1]]*Y[w[2]]/(total_X*total_Y)
            print(f"{w}\t{count}\t{a}\t{b}\t{a/b}")
    absent_trimers = list(absent_trimers)
    absent_trimers.sort()
    return absent_trimers


def P_pattern(seq):
    p=phasing_GXY_pattern(seq)
    if p<0:
        return (-1, -1)
    weak_P=0
    strong_P=0
    for i in range((p+2)%3, len(seq), 3):
        if seq[i] == 'P':
            strong_P+= 1
    for i in range((p+1)%3, len(seq), 3):
        if seq[i] == 'P':
            weak_P+= 1
    return (strong_P, weak_P)
    
def Pmask_distance(seq1, seq2):
    if len(seq1)!=len(seq2):
        return -1
    r=0
    for c1, c2 in zip(seq1, seq2):
        if c1=='P' and c2=='P':
            None
        elif c1=='P'or c2=='P':
            r=r+1
    return r

# a changer: prendre en paramètre un string, plutôt qu'une list of strings
def compute_X_Y_positions(seq):
    # mettre les 3 lignes ci-dessous ailleurs
    phase=phasing_GXY_pattern(seq)
    return range(phase+1, len(seq), 3), range(phase + 2, len(seq), 3)

# m is a marker
def PTM_adequacy (m):
    strong_P, _ = P_pattern(m.sequence())
    nb_of_H = compute_masses.number_of_H(m.PTM())
    return strong_P == nb_of_H

