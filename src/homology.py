"""
homology.py
"""

from collections import Counter

from src import utils as ut
from src import compute_masses as mass
from src import taxonomy
from src import sequences 
from src import markers
from src import collagen
from src import config

def hamming_distance(seq1, seq2):
    if len(seq1)!=len(seq2):
        return -1
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

    
def find_markers_single_sequence(seq, set_of_digested_peptides, dict_of_model_markers, set_of_markers, taxo):
    # seq: protein (object defined in sequences.py)
    # dict_of_model_markers, dictionary, key: peptide sequence [string], value: set of 3-uplets (PTM, code, gene_name)
    # dict_of_taxid, key: (sequence,PTM), value: set of taxid
    #set_of_markers: ?

    # for each marker of dict_of_model_markers (characterized by a *code*) find the best location in the sequence seq.
    # output: set of markers

    helical_start=sequences.helical_region(seq)[0]
    set_of_found_codes=set() # contains the set of codes for  model markers found in the current sequence
    set_of_raw_digested_peptides={pep.sequence() for pep in set_of_digested_peptides}
    
    # markers already present in the set of markers (same taxid, PTM and code)
    set_of_new_markers={m for m in set_of_markers if seq.taxid()==m.taxid() and seq.protein()==m.protein()}
    # TO DO: add a test for GN
    for m in set_of_new_markers: # markers for the same taxid
        pos=(seq.sequence()).find(m.sequence())
        if pos>=0:
            m.field["Begin"]=pos+1
            m.field["SeqID"]=seq.seqid()
            m.field["Length"]=len(m.sequence())
            m.field["End"]=pos+len(m.sequence())
            if helical_start is None or helical_start>pos:
                m.field["Hel"]=None
            else:
                m.field["Hel"]=pos-helical_start+2
        if "Mass" not in m.field:
            m.field["Mass"]=mass.peptide_mass_with_PTM(m.sequence(),m.PTM())
        set_of_found_codes.add(m.code())
        if m.sequence() not in set_of_raw_digested_peptides:
            m.field["Digestion"]="No"
            
    found_markers={} # key: name of the marker, value: triplet (position, hamming distance, Pmask, marker sequence)

    #markers found with exact match in some other organism 
    for marker_seq in dict_of_model_markers:
        pos=(seq.sequence()).find(marker_seq)
        if (pos>=0):
            set_of_codes={s[1] for s in dict_of_model_markers[marker_seq] if seq.protein()==s[2] }
            for code in set_of_codes:
                if code not in set_of_found_codes:
                    found_markers[code]=(pos,0,0,marker_seq)

    set_of_found_codes.update(found_markers.keys())
        
    # other markers
    for marker_seq in dict_of_model_markers.keys():      
        set_of_protein_names={s[2] for s in dict_of_model_markers[marker_seq]}
        if  seq.protein() not in set_of_protein_names:
            continue
        set_of_codes={s[1]  for s in dict_of_model_markers[marker_seq]}
        if set_of_codes <=  set_of_found_codes:
            continue
        l=len(marker_seq)
        for pos in range(0, len(seq.sequence())-l):
            d=hamming_distance((seq.sequence())[pos:pos+l], marker_seq)
            if d<len(marker_seq)/10+1:
                d2=collagen.Pmask_distance((seq.sequence())[pos:pos+l], marker_seq)
                for  code in set_of_codes:
                    if (code not in found_markers) :
                        found_markers[code]=(pos,d,d2, marker_seq)
                    elif d<found_markers[code][1] or (d==found_markers[code][1] and d2<found_markers[code][2]):
                        found_markers[code]=(pos,d,d2, marker_seq)
    
    for code in found_markers:
        (pos,d,d2, marker_seq)=found_markers[code]
        l=len(marker_seq)
        new_sequence= seq.sequence()[pos:pos+l]
        # this loop seems weird, but is actually useful.
        # this is due to the fact that a given marker_name can have multiple PTMs. In this case, each possibility gives rise to a distinct marker.  
        for model_marker in dict_of_model_markers[marker_seq]:
            if model_marker[1]!=code:
                continue
            dict={}
            dict["OX"]=seq.taxid()
            dict["OS"]=seq.taxon_name()
            dict["SeqID"]=seq.seqid()
            dict["Sequence"]=new_sequence
            dict["Begin"]= pos+1
            dict["Length"]=l
            if helical_start is not None and pos>helical_start:
                dict["Hel"]=pos-helical_start+2
            dict["End"]=pos+l
            dict["Rank"]="species"
            dict["Marker"]=model_marker[1]
            dict["GN"]=model_marker[2]
            dict["Status"]="Genetic"
            if new_sequence not in set_of_raw_digested_peptides:
                dict["Digestion"]="No"
            else:
                dict["Digestion"]="Yes"
            if d==0:
                dict["Comment"]="Homology : "+ seq.sequence()[pos-1]+" - peptide - "+ seq.sequence()[pos+l]+", exact match "
            elif d==1:
                dict["Comment"]="Homology : "+ seq.sequence()[pos-1]+" - peptide - "+ seq.sequence()[pos+l]+ ", 1 mismatch with " + marker_seq+" "
            else:
                dict["Comment"]="Homology : "+ seq.sequence()[pos-1]+" - peptide - "+ seq.sequence()[pos+l]+ ", "+str(d)+ " mismatches with " + marker_seq+" "
            if taxo:
                set_of_potential_taxids={m.taxid() for m in set_of_markers if m.sequence()==marker_seq}
                dict["Comment"]=dict["Comment"]+taxonomy.find_closest_ID(seq.taxid(),set_of_potential_taxids,taxo)
            else:
                dict["Comment"]=dict["Comment"]+". "
            if collagen.P_pattern(new_sequence)[0]==collagen.P_pattern(marker_seq)[0]:
                dict["PTM"]=model_marker[0]
            else:
                dict["Comment"]=dict["Comment"]+ "PTM updated. "
                
            new_marker=markers.Marker(field=dict)
            set_of_new_markers.add(new_marker)
    return mass.add_PTM_or_masses_to_markers(set_of_new_markers)
    
def find_helical_position(set_of_markers):
    dict_codes={m.code():[] for m in set_of_markers}
    for m in set_of_markers:
        if m.helical() is None or str(m.helical())=="0":
            continue
        dict_codes[m.code()].append(m.helical())
    for code in dict_codes:
        if len(dict_codes[code])==0:
            dict_codes[code]=0
            continue
        counts = Counter(dict_codes[code])
        most_common, freq = counts.most_common(1)[0]  # Get the most frequent element and its count
        if freq > len(dict_codes[code]) / 2:
            dict_codes[code]=most_common
        else:
            dict_codes[code]=0
    return dict_codes
    
def check_quality(set_of_markers, set_of_new_markers):
    dict_codes=find_helical_position(set_of_markers)
    set_of_sequences={m.sequence() for m in set_of_markers}
    set_of_new_sequences={m.sequence() for m in set_of_new_markers}
    dict_of_new_sequences={seq:len({m.taxon_name() for m in set_of_new_markers if m.sequence()==seq}) for seq in set_of_new_sequences}
    for m in set_of_new_markers:
        exist_sequence=(m.sequence() in set_of_sequences) or (dict_of_new_sequences[m.sequence()]>1)
        exist_code= dict_codes[m.code()]>0 and m.helical() is not None
        if not exist_code:
            m.field["Quality"]=0
        elif m.helical()==dict_codes[m.code()]:
            m.field["Quality"]=2
        elif abs(m.helical()-dict_codes[m.code()])<4:
            m.field["Quality"]=1
        else:
            m.field["Quality"]=0
        if exist_sequence:
            m.field["Quality"]+=1
    return set_of_new_markers

def find_markers_all_sequences(set_of_sequences, set_of_markers, taxo):
    # set_of_sequences: set of target fasta sequences[object defined in sequences.py]
    
    # construction of dict_of_model_markers and dict_of_taxid
    #dict_of_model_markers, key: peptide sequence [string], value:set of 3-uplets (PTM, code, gene_name)
    #dict_of_taxid, key: (sequence,PTM), value: set of taxid 
    dict_of_model_markers={}
    dict_of_taxid={}
    marker_count=0
    for m in set_of_markers: 
        if m.sequence() is None or len(m.sequence())==0 :
            continue
        else:
            sequence=m.sequence().strip()
        if m.code() is None or len(m.code())==0 :
            marker_count+=1
            peptide_name="M"+str(marker_count)
        else:
            peptide_name=m.code()
        model_marker=(ut.standard(m.PTM()), ut.standard(peptide_name), ut.standard_upper(m.protein()))
        ut.update_dictoset(dict_of_model_markers, sequence, {model_marker})
        
    set_of_new_markers=set()
    number_of_missed_cleavages=int(config.parse_config_file()["number_of_missed_cleavages"])
    for seq in set_of_sequences:
        set_of_digested_peptides=sequences.in_silico_digestion({seq}, number_of_missed_cleavages, 12, 33, False)
        s=find_markers_single_sequence(seq, set_of_digested_peptides, dict_of_model_markers, set_of_markers, taxo)
        set_of_new_markers.update(s)
        
    list_of_new_markers=markers.sort_and_merge(set_of_new_markers)
    list_of_new_markers=check_quality(set_of_markers, list_of_new_markers)

    return list_of_new_markers


