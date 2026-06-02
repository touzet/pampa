import json
from src import  sequences, compute_masses, markers

def assign_sequences_to_internal_clades(taxid, set_of_sequences, taxo, sequence_dict):
    if taxo.is_leaf(taxid):
        for seq in set_of_sequences:
            if seq.taxid() == taxid:
                sequence_dict[taxid]=seq.sequence(), 1
                return
        sequence_dict[taxid] = None
        return
    elif taxo.number_of_children(taxid)==1:
        child_taxid = next(iter(taxo.children[taxid]))
        assign_sequences_to_internal_clades(child_taxid, set_of_sequences,  taxo, sequence_dict)
        sequence_dict[taxid]=sequence_dict[child_taxid]
        return
    else:
        set_of_child_sequences=set()
        for c in taxo.children[taxid]:
            assign_sequences_to_internal_clades(c, set_of_sequences, taxo, sequence_dict)
            if sequence_dict[c] is not None:
                set_of_child_sequences.add(sequence_dict[c][0])
        max_length=max({len(seq) for seq in set_of_child_sequences if seq is not None})
        consensus_sequence=""
        for position in range(max_length):
            set_of_aa={seq[position] for seq in set_of_child_sequences if position<len(seq)}
            if  len(set_of_aa)==1:
                consensus_sequence += set_of_aa.pop()
            else:
                consensus_sequence += 'X'
    sequence_dict[taxid]=consensus_sequence, sum(sequence_dict[c][1] for c in taxo.children[taxid])
    return

def identify_conserved_markers(set_of_sequences, config_digestion):
    CP_dict = {}
    for seq in set_of_sequences :
        set_of_peptides = sequences.raw_in_silico_digestion(seq.sequence(), config_digestion)
        set_of_markers = {markers.Marker(field={"Sequence": peptide}) for peptide in set_of_peptides if
                          'X' not in peptide}
        # TO DO : retirer les peptides qui correspondent à des peptides marqueurs
        set_of_markers = compute_masses.add_PTM_or_masses_to_markers(set_of_markers)
        CP_dict[seq.taxid()] = sorted({m.mass() for m in set_of_markers})
    return CP_dict



def identify_conserved_markers_old(taxid,  set_of_sequences, taxo, config_digestion):
    sequence_dict = {}
    CP_dict = {}
    set_of_taxid = {seq.taxid() for seq in set_of_sequences}
    taxo2, _ = taxo.intersection(set_of_taxid)
    assign_sequences_to_internal_clades(taxid, set_of_sequences, taxo2, sequence_dict)
    for t in sequence_dict:
        set_of_peptides = sequences.raw_in_silico_digestion(sequence_dict[t][0], config_digestion)
        set_of_markers = {markers.Marker(field={"Sequence": peptide}) for peptide in set_of_peptides if
                          'X' not in peptide}
        # TO DO : retirer les peptides qui correspondent à des peptides marqueurs
        set_of_markers = compute_masses.add_PTM_or_masses_to_markers(set_of_markers)
        CP_dict[t] = sorted({m.mass() for m in set_of_markers})
    return CP_dict

def conserved_markers_initialisation(set_of_sequences, config_digestion, seed_name):
    file_name="conserved_"+seed_name+".json"
    CP_dict=identify_conserved_markers(set_of_sequences,  config_digestion)
    with open(file_name, 'w') as f:
        json.dump(CP_dict, f, ensure_ascii=False, indent=4)
    return [file_name]

"""
# find conserved markers for a set of FASTA sequences
def conserved_markers(taxid,  set_of_sequences, taxo, config):
    config_digestion = conf.config_digestion(config)
    sequence_dict={}
    set_of_taxid = {seq.taxid() for seq in set_of_sequences}
    taxo2, _ = taxo.intersection(set_of_taxid)
        #for taxid in taxo2.root:
        #    assign_sequences_to_internal_clades(taxid, set_of_prot_sequences, taxo2, sequence_dict)
    assign_sequences_to_internal_clades(taxid, set_of_sequences, taxo2, sequence_dict)
    set_of_peptides = sequences.raw_in_silico_digestion(sequence_dict[taxid][0], config_digestion)
    set_of_markers = {markers.Marker(field={"Sequence": peptide}) for peptide in set_of_peptides if
                                  'X' not in peptide}
            # TO DO : retirer les peptides qui correspondent à des peptides marqueurs
    set_of_markers = compute_masses.add_PTM_or_masses_to_markers(set_of_markers)
    return {m.mass() for m in set_of_markers}
"""

def find_CP_nodes(taxid, taxo, CP_dict):
    parent=taxid
    case=0 # 0: not found , -1 root, 1 or more found
    while case==0 :
        if parent in taxo.root:
            case = -1
        else:
            parent=taxo.parent[parent]
            nodes= taxo.descendants[parent] & CP_dict.keys()
            case = len(nodes)
    if case==-1 :
        return set()
    else:
        return set().union(*(CP_dict[node] for node in nodes))


""""
def conserved_markers_target_taxid(target_taxid, close_species, set_of_sequences, taxo, config):
    set_of_proteins={seq.protein() for seq in set_of_sequences}
    set_of_masses=set()
    for prot in set_of_proteins:
        set_of_prot_sequences={seq for seq in set_of_sequences if seq.protein() == prot}
        set_of_taxid = {seq.taxid() for seq in set_of_prot_sequences}
        taxo2, _ = taxo.intersection(set_of_taxid | close_species)
        print(close_species | {target_taxid})
        parent=find_parent(close_species | {target_taxid}, taxo2)
        set_of_masses.update(conserved_markers(parent, set_of_prot_sequences, taxo2, config))
    return set_of_masses
"""

def load_conserved_markers(file_name):
    with open(file_name, "r") as f:
        dict_conserved = json.load(f)
    return dict_conserved

def compute_CPs(taxid, taxo, set_of_sequences, config,  file_name="custom"):
    CP_dict=identify_conserved_markers(set_of_sequences, config)
    with open("conserved_"+file_name+".json", 'w') as f:
        json.dump(CP_dict, f, ensure_ascii=False, indent=4)
    return find_CP_nodes(taxid, taxo, CP_dict)

#file_name: json CP file
def extract_CPs_from_file(taxid, file_name, taxo):
    CP_dict=load_conserved_markers(file_name)
    return find_CP_nodes(taxid, taxo, CP_dict)

