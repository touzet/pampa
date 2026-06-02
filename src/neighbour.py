from collections import defaultdict

from src import compute_masses
from src import markers
from src import message
from src import utils
from src import assignment
from src import collagen
from src import sequences
from src import reconstruct


def add_trimer_scores(set_of_markers):
    proh_trimers=collagen.prohibited_trimers()
    for m in set_of_markers:
        p=collagen.phasing_GXY_pattern(m.sequence())
        if p==-1:
            m.field["Proh"]=False
            continue
        set_of_trimers={m.sequence()[i:i+3] for i in range(p, len(m.sequence()) - 2, 3)}
        m.field["Proh"]=not(any(trimer in proh_trimers  for trimer in set_of_trimers))


def add_conserved_peak_scores(set_of_markers, set_of_masses):
    for m in set_of_markers:
        if len(set_of_masses)==0:
            m.field["CP"] = 1
        else :
            m.field["CP"] = 2
        for mass in set_of_masses:
            if utils.matching_masses(mass, m.mass(), 0.1):
                m.field["CP"] = 0

def add_existence_scores(set_of_markers, set_of_all_markers):
    set_of_peptides={m.sequence() for m in set_of_all_markers}
    for m in set_of_markers:
        m.field["Exists"]=m.sequence() in set_of_peptides

# tous les marqueurs sont pour le même code
def find_variable_positions(set_of_markers):
    set_of_sequences = {m.sequence() for m in set_of_markers}
    if len(set_of_sequences)<2 or len({len(seq) for seq in set_of_sequences})>1:
        return set()
    marker_length= len(next(iter(set_of_sequences)))
    return {i for i in range(marker_length) if len({seq[i] for seq in set_of_sequences}) > 1}

def compute_variable_positions(set_of_markers, gamma_count):
    dict_of_variable_positions={}
    set_of_codes={m.code() for m in set_of_markers}
    for code in set_of_codes:
        set_of_selected_markers={m for m in set_of_markers if m.code()==code}
        dict_of_variable_positions[code]=find_variable_positions(set_of_selected_markers) | gamma_count[code].keys()
    return dict_of_variable_positions

# set_of_close_markers est supposé inclus dans set_of_all_markers (ce qui garantit la longueur)
def neighbour_for_one_sequence(seq, dict_of_orphan_markers, dict_of_variable_positions, target_taxid, matrices, gamma_count, config_digestion, set_of_codes_for_deamidation):
    code=dict_of_orphan_markers[seq][0]
    prot = dict_of_orphan_markers[seq][1]
    deamidation= True if (set_of_codes_for_deamidation is None) else  code in set_of_codes_for_deamidation
    if prot not in matrices.keys():
        message.warning("No available substitution scoring matrix for " + str(prot))
        return set()
    set_of_variable_positions= dict_of_variable_positions[code]
    set_of_neighbours=set()
    phase=collagen.phasing_GXY_pattern(seq)
    # To do: add error when phase <0
    for position in set_of_variable_positions:
        pos_label= 'Y' if (position-phase)%3==2 else 'X'
        if pos_label in matrices[prot].keys():
            matrix = matrices[prot][pos_label]
        else:
            message.warning("No available substitution scoring matrix for "+str(prot)+" with label "+pos_label)
            continue
        for aa in list(matrix.columns):
            if aa != seq[position] and {aa, seq[position]} - {'I', 'L'} and matrix.loc[aa][seq[position]]>0.0 and matrix.loc[seq[position]][seq[position]]+matrix.loc[aa][aa] > 4:
                new_peptide = seq[:position] + aa + seq[position + 1:]
                if not sequences.is_digested_peptide(new_peptide, config_digestion):
                    continue
                set_of_ptm = compute_masses.update_PTM(new_peptide, seq, dict_of_orphan_markers[seq][2], deamidation)
                for ptm in set_of_ptm:
                    new_m=markers.Marker(field={})
                    new_m.field["OX"]=target_taxid
                    new_m.field["Status"]="MS"
                    new_m.field["Sequence"] = new_peptide
                    new_m.field["Marker"] = code
                    new_m.field["PTM"] = ptm
                    new_m.field["Mass"] = compute_masses.peptide_mass_with_PTM(new_peptide, ptm )
                    new_m.field["GN"]= prot
                    new_m.field["Subst"]= matrix.loc[aa][seq[position]]
                    new_m.field["Mutability"]= matrix.loc[seq[position]][seq[position]]+matrix.loc[aa][aa]
                    new_m.field["Gamma"]=abs(gamma_count[code][position]) if position in gamma_count[code] else 0
                    new_m.field["Comment"]=" One mutation with "+seq+ " from "+ ", ".join(dict_of_orphan_markers[seq][3])
                    set_of_neighbours.add(new_m)
        #set_of_neighbours=markers.sort_and_merge(set_of_neighbours)
    return set_of_neighbours

def peak_filter_neighbours(set_of_markers, list_of_spectra, resolution, min_nb_spectra):
    mass_list=markers.sort_markers_by_mass(set_of_markers)
    set_of_new_markers=set()
    dict_markers={}
    for spectrum in list_of_spectra:
        peak_to_markers = assignment.find_matching_peaks_and_markers(spectrum, mass_list, resolution)
        for p in peak_to_markers:
            for m in peak_to_markers[p]:
                utils.update_dictoset(dict_markers, m, {(spectrum, p)})
    for m in dict_markers:
        set_of_found_spectra = dict_markers[m]
        if  len(set_of_found_spectra)<len(list_of_spectra) * 0.12 :
            continue
        matching_spectra = ", ".join(["(" + str(round(float(sp[1].mass), 2)) + ", " + str(int(sp[1].intensity))+ ") in " + sp[0].name for sp in
                                    set_of_found_spectra])
        markers.update_comment(m, "Peak " + matching_spectra+". ")
        m.field["Spectra"] = len(set_of_found_spectra)
        m.field["Intensity"] = round(sum([p[1].intensity/p[0].median for p in set_of_found_spectra]) / len(set_of_found_spectra),1)
        m.field["PTM_adequacy"] = collagen.PTM_adequacy(m)
        m.field["Neighbour"] = 0
        m.field["Delta_mz"]= reconstruct.delta_mz(m.mass(), set_of_found_spectra, resolution)
        set_of_new_markers.add(m)
    return set_of_new_markers

def find_candidate_sequences_for_orphan_markers(dict_of_orphan_markers, target_taxid, list_of_orphan_spectra, resolution, min_nb_of_spectra, matrices, gamma_count, set_of_all_markers, dict_of_variable_positions, set_of_codes_for_deamidation, config_digestion, dict_of_conserved_peaks):
    set_of_new_markers=set()
    for seq in dict_of_orphan_markers:
        set_of_seq_markers=neighbour_for_one_sequence(seq, dict_of_orphan_markers, dict_of_variable_positions, target_taxid, matrices, gamma_count, config_digestion, set_of_codes_for_deamidation)
        set_of_seq_markers=peak_filter_neighbours(set_of_seq_markers, list_of_orphan_spectra, resolution, min_nb_of_spectra)
        set_of_new_markers.update(set_of_seq_markers)
    set_of_new_markers = reconstruct.find_paired_markers(set_of_new_markers)
    add_trimer_scores(set_of_new_markers)
    add_conserved_peak_scores(set_of_new_markers, dict_of_conserved_peaks)
    add_existence_scores(set_of_new_markers, set_of_all_markers)#set_of_new_markers = pareto.add_pareto_score(set_of_new_markers)
    return set_of_new_markers

def merge_markers(list_of_markers):
    score_max=max({m.field["score"] for m in list_of_markers})
    list_of_markers[0].field["score"]=score_max
    return list_of_markers[0]


def remove_duplicates(set_of_markers):
    best = {}
    for m in set_of_markers:
        seq, ptm = m.sequence(), m.PTM()
        if (seq,ptm) not in best or float(m.field["score"]) > float(best[(seq, ptm)].field["score"]):
            best[(seq, ptm)] = m
    return set(best.values())

def remove_duplicate(set_of_markers):
    groups = defaultdict(list)
    for m in set_of_markers:
        groups[(m.sequence(), m.PTM())].append(m)
    set_of_new_markers=set()
    for l in groups.values():
        if len(l)==1:
            set_of_new_markers.add(l[0])
        else:
            set_of_new_markers.add(merge_markers(l))
    return set_of_new_markers