"""
known_markers.py            

"""

from src import utils as ut
from src import compute_masses as mass
from src import sequences 
from src import markers

def hamming_distance(seq1, seq2):
    if len(seq11)!=len(seq2):
        return -1
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))
        
def find_markers_single_sequence(seq,  set_of_digested_peptides, dict_of_model_markers, set_of_markers):
    # seq: protein (object defined in sequences.py)
    # dict_of_model_markers: dictionary, key: peptide sequence [string], value: set of 3-uplets (PTM, code, gene_name)

    # for each marker of dict_of_model_markers (characterized by a *code*) find the best location in the sequence seq.
    # output: set of markers

    set_of_found_codes=set() # contains the codes of  model markers found in the current sequence
    set_of_found_sequences=set() # contains the sequences of the model markers present in the current sequence
    set_of_exact_markers=set() # contains the name of markers that have been found exactly, with no error 
    found_markers={} # key: name of the marker, value: triplet (position, hamming distance, marker sequence)

    for marker_seq in dict_of_model_markers:
        pos=(seq.sequence).find(marker_seq)
        if (pos>=0):
            set_of_codes={s[1] for s in dict_of_model_markers[marker_seq] if seq.protein==s[2] }
            for code in set_of_codes:
                found_markers[code]=(pos,0,marker_seq)

    set_of_found_codes=found_markers.keys()
                
    for marker_seq in dict_of_model_markers.keys():      
        set_of_protein_names={s[2] for s in dict_of_model_markers[marker_seq]}
        if  seq.protein not in set_of_protein_names:
            continue
        set_of_codes={s[1] for s in dict_of_model_markers[marker_seq]  }
        if set_of_codes <   set_of_found_codes:
            continue
        l=len(marker_seq)
        for pos in range(0, len(seq.sequence)-l):
            d=hamming_distance((seq.sequence)[pos:pos+l], marker_seq)
            if d<len(marker_seq)/10+1:
                set_of_codes={s[1] for s in dict_of_model_markers[marker_seq]}
                for  code in set_of_codes:
                    if (code not in found_markers) or (d<found_markers[code][1]):
                        found_markers[code]=(pos,d,marker_seq)
                    
    set_of_new_markers=set()
    set_of_raw_digested_peptides={pep.sequence for pep in set_of_digested_peptides}
    for code in found_markers:
        (pos,d,marker_seq)=found_markers[code]
        l=len(marker_seq)
        new_sequence= seq.sequence[pos:pos+l]
        if new_sequence in set_of_raw_digested_peptides:
            found=True
            new_cleavage=False
        else:
            found=False
            for peptide in set_of_raw_digested_peptides:
                seq_pos=peptide.find(new_sequence)
                if (seq_pos>=0):
                    new_sequence=peptide
                    pos=pos+seq_pos
                    l=len(peptide)
                    found=True
                    new_cleavage=True
        if not found:
            continue
        # this loop seems weird, but is actually useful.
        # this is due to the fact that a given marker_name can have multiple PTM. In this case, each possibility gives rise to a distinct marker.  
        for model_marker in dict_of_model_markers[marker_seq]:
            if model_marker[1]!=code:
                continue
            new_marker=markers.Marker()
            new_marker.taxid=seq.taxid
            new_marker.taxon_name=seq.taxon_name
            new_marker.seqid=seq.seqid
            new_marker.sequence=new_sequence
            new_marker.begin=str(pos)
            new_marker.end=str(pos+l-1)
            new_marker.rank="species"
            new_marker.ptm=model_marker[0]
            new_marker.masses={mass.peptide_mass_with_PTM(seq.sequence[pos:pos+l],model_marker[0])} # à modifier pour conserver la masse
            new_marker.code=model_marker[1]
            new_marker.protein=model_marker[2]
            new_marker.comment=seq.sequence[pos-1]+" - peptide - "+ seq.sequence[pos+l]+", #mismatch "+str(d)
            if d>0:
                new_marker.comment=new_marker.comment +  " - coming from " + marker_seq
            if new_cleavage:
                new_marker.comment=new_marker.comment +" - new cleavage site"
            set_of_new_markers.add(new_marker)
            
    return set_of_new_markers

def find_markers_all_sequences(set_of_sequences, set_of_markers):
    # set_of_sequences: set of sequences[object defined in sequences.py]
    # dict_of_model_markers: dictionary, key: peptide sequence [string], value:set of 3-uplets (PTM, code, gene_name)
    dict_of_model_markers={}
    marker_count=0
    for m in set_of_markers: 
            if len(m.sequence)==0 :
                 continue
            if len(m.code)==0:
                marker_count+=1
                peptide_name="M"+str(marker_count)
            else:
                peptide_name=m.code
            model_marker=(m.ptm.strip(), peptide_name.strip(), m.protein.strip().upper())
            ut.update_dictoset(dict_of_model_markers, m.sequence.strip(), {model_marker})

    set_of_markers=set()
    for seq in set_of_sequences:
        set_of_digested_peptides=sequences.in_silico_digestion({seq},1, 12, 33, False)
        s=find_markers_single_sequence(seq, set_of_digested_peptides, dict_of_model_markers, set_of_markers)
        set_of_markers.update(s)

    return set_of_markers


