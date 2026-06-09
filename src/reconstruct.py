import copy
import csv
from src import utils
from src import assignment
from src import markers
from src import taxonomy as ta
from src import message

def parse_target_taxid(targetfile, list_of_all_spectra):
    list_of_all_spectra_names = {s.name for s in list_of_all_spectra}
    with open(targetfile, newline='', encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        result = {}
        for row in reader:
            if len(row) < 2:
                continue  # skip incomplete lines
            value = row[0].strip()
            if value not in list_of_all_spectra_names:
                message.warning("Spectral file *"+value+"* not found in directory. Ignored.")
            else :
                key = row[1].strip()
                for s in list_of_all_spectra:
                    if s.name==value:
                        result.setdefault(key, []).append(s)
    return result

def initialize_target_taxids(target, target_taxid_file, list_of_spectra, full_taxonomy):
    if target:
        return {ta.find_taxid(target, full_taxonomy):list_of_spectra}
    else:
        return parse_target_taxid(target_taxid_file, list_of_spectra)

def create_matching_marker(target_taxid, m):
    new_m = markers.Marker(field={})
    new_m.field["Sequence"] = m.sequence()
    new_m.field["Marker"] = m.code()
    new_m.field["PTM"] = m.PTM()
    new_m.field["Mass"] = m.mass()
    new_m.field["GN"] = m.protein()
    new_m.field["Status"] = "MS"
    new_m.field["OX"] = target_taxid
    new_m.field["OS"] = None
    new_m.field["Comment"] = ", match with "+m.taxon_name()+", sequence " + m.sequence()
    return new_m

# construit l'ensemble des marqueurs des espèces proches. Les marqueurs de même code, séquence et PTM sont fusionnés
def build_candidate_markers(target_taxid, set_of_taxid, set_of_markers):
    set_of_selected_markers=set()
    for m in set_of_markers:
        if m.taxid() not in set_of_taxid:
            continue
        found=False
        for m2 in set_of_selected_markers:
            if (m.mass(), m.sequence(), m.code()) == (m2.mass(), m2.sequence(), m2.code()):
                markers.post_comment(m2, ", " + m.taxon_name())
                found=True
        if not found:
            new_m=create_matching_marker(target_taxid, m)
            set_of_selected_markers.add(new_m)
    return set_of_selected_markers

def delta_mz(mass, set_of_spectra, resolution):
    return round(utils.margin_tolerance(mass,resolution)- sum([abs(mass -p[1].mass) for p in set_of_spectra]) / len(set_of_spectra), 4)

def spectra_threshold(list_of_spectra, found_spectra, resolution):
    return len(found_spectra)>0.12 * len(list_of_spectra)

def select_one_marker_per_sequence(set_of_markers):
    selected = {}
    for m in set_of_markers:
        if m.sequence() not in selected:
            selected[m.sequence()] = m
    return set(selected.values())

def find_known_markers_in_spectra(target_taxid, set_of_markers, list_of_spectra, resolution):
    mass_list=markers.sort_markers_by_mass(set_of_markers)
    set_of_matching_markers=set()
    list_of_orphan_spectra=copy.deepcopy(list_of_spectra) # argument à retirer
    dict_markers={}
    for spectrum in list_of_spectra:
        peak_to_markers=assignment.find_matching_peaks_and_markers(spectrum, mass_list, resolution)
        for p in peak_to_markers:
            for m in peak_to_markers[p]:
                utils.update_dictoset(dict_markers, (m.mass(), m.code()), {(spectrum,p,m)})
    for (mass, code) in dict_markers:
        set_of_found_spectra={(spectrum,p) for (spectrum,p,m) in dict_markers[mass, code]}
        if not spectra_threshold(list_of_spectra, set_of_found_spectra, resolution):
            continue
        set_of_taxons={m.taxon_name() for (s,p,m) in dict_markers[mass, code]}
        set_of_selected_markers=select_one_marker_per_sequence({m for (s,p,m) in dict_markers[mass, code]})
        #set_of_ptm={m.PTM() for (s,p,m) in dict_markers[mass, code]}
        #set_of_sequences={m.sequence() for (s,p,m) in dict_markers[mass, code]}
        #set_of_proteins={m.protein() for (s,p,m) in dict_markers[mass, code]}
        matching_spectra="".join(["("+str(round(float(sp[1].mass),2))+", "+str(int(sp[1].intensity))+") in "+sp[0].name+" " for sp in set_of_found_spectra])
        #GN = set_of_proteins.pop()
        #PTM = set_of_ptm.pop()
        for m in set_of_selected_markers :
            new_m=markers.Marker(field={})
            new_m.field["OX"]=target_taxid
            new_m.field["Status"]="MS"
            new_m.field["Sequence"] = m.sequence()
            new_m.field["Marker"] = code
            new_m.field["Mass"] = mass
            new_m.field["GN"] = m.protein()
            new_m.field["PTM"] = m.PTM()
            new_m.field["Spectra"]=len(set_of_found_spectra)
            new_m.field["Intensity"]= round(sum([p[1].intensity/p[0].median for p in set_of_found_spectra]) / len(set_of_found_spectra),1)
            new_m.field["OS"]=None
            new_m.field["Comment"] = "Peak "+matching_spectra + "matching with "+", ".join(set_of_taxons)
            new_m.field["Neighbour"] = 10
            new_m.field["Delta_mz"] = delta_mz(mass, set_of_found_spectra, resolution)
            new_m.field["Exists"] = True
            set_of_matching_markers.add(new_m)
    return set_of_matching_markers, list_of_orphan_spectra

# set_of_markers : contains the set of markers coming from close_species
def find_orphan_peaks_and_peptides(target_taxid, list_of_spectra, set_of_markers, error):
    set_of_codes={m.code() for m in set_of_markers}
    set_of_matching_markers, list_of_orphan_spectra= find_known_markers_in_spectra(target_taxid, set_of_markers, list_of_spectra, error)
    #set_of_orphan_codes=set_of_codes.difference({m.code() for m in set_of_matching_markers})
    set_of_orphan_codes=set_of_codes
    set_of_orphan_sequences={m.sequence() for m in set_of_markers if m.code() in set_of_orphan_codes}
    dict_of_orphan_markers={}
    for sequence in set_of_orphan_sequences:
        PTMs=set()
        taxons=set()
        for m in set_of_markers:
            if m.sequence()==sequence:
                code=m.code()
                prot=m.protein()
                PTMs.add(m.PTM())
                taxons.add(m.taxon_name())
        dict_of_orphan_markers[sequence]=(code, prot, PTMs, taxons)
        set_of_matching_markers = find_paired_markers (set_of_matching_markers)
    return set_of_matching_markers, dict_of_orphan_markers, list_of_orphan_spectra

def find_paired_markers(set_of_markers):
    for m in set_of_markers:
        m.field["Paired"] = 0
        for m2 in set_of_markers:
            if  m.taxid()==m2.taxid() and m.mass()!= m2.mass() and m.sequence()==m2.sequence():
                m.field["Paired"] = 1
    return set_of_markers

#def process_overlap(marker1,marker2):

"""
def post_process(set_of_markers, overlap_list):
    if len(overlap_list)==0:
        return
    for (code1, code2) in overlap_list:
        for m1, m2 in combination(set_of_markers,2):
            if (m1.code(),m2.code())==(code1,code2):
                process_overlap(m1,m2)
"""