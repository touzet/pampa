from collections import defaultdict, Counter

from src import compute_masses, sequences, markers, message, utils, taxonomy

def is_complete_marker(m, required_fields):
    return set(required_fields).issubset(m.field)

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



def add_marker_comment(list_of_markers, comment):
    for m in list_of_markers:
        m.field["Comment"] = comment

# use mass to find sequence
def add_sequences_from_mass(set_of_markers, set_of_sequences, resolution, config_digestion, taxo, taxonomy_ranks):
    set_of_target_markers = {m for m in set_of_markers if (
                m.taxid() or m.taxon_name()) and m.mass() and m.sequence() is None}  # missing sequences for the set of markers
    target_mass_list = [(m.mass(), m) for m in set_of_target_markers]
    target_mass_list.sort(key=lambda x: x[0])
    set_of_target_sequences = set()
    for m in set_of_target_markers:
        set_of_target_sequences.update(markers.find_matching_sequences(m, set_of_sequences, taxo,taxonomy_ranks))
    set_of_new_markers = {m for m in set_of_markers if m not in set_of_target_markers}
    set_of_denovo_markers = compute_masses.add_PTM_or_masses_to_markers(
        sequences.in_silico_digestion(set_of_target_sequences, config_digestion), True, True)
    denovo_mass_list = [(m.mass(), m) for m in set_of_denovo_markers]  # masses of all tryptic peptides
    denovo_mass_list.sort(key=lambda x: x[0])
    founded_masses = {y: 0 for (x, y) in target_mass_list}
    # for target_mass_list
    for mass in target_mass_list:
        for j in range(0, len(denovo_mass_list)):
            if utils.matching_masses(mass[0], denovo_mass_list[j][0], resolution) and (
                    mass[1].taxid() == denovo_mass_list[j][1].taxid() or utils.equiv(mass[1].taxon_name(),
                                                                                  denovo_mass_list[j][
                                                                                      1].taxon_name())):  # on a un match
                m2 = mass[1]
                m = denovo_mass_list[j][1]
                founded_masses[m2] += 1
                dict = {x: m.field[x] for x in m.field}
                if m2.code() is None:
                    dict["Marker"] = mass[1].code()
                else:
                    dict["Marker"] = m2.code()
                dict["Comment"] = "Supplement from "+m.field["magic_number"]+". "+m.comment() + "Sequence deduced from target mass. "
                new_marker = markers.Marker(field=dict)
                set_of_new_markers.add(new_marker)
                continue
            if denovo_mass_list[j][0] - mass[0] > 1:
                continue
    for m in founded_masses:
        if founded_masses[m] == 0:
            markers.update_comment(m, "No matching peptide found. ")
            set_of_new_markers.add(m)
            message.warning(m.taxon_name() + ": no matching peptide found for mass " + str(m.mass()) + ".")
        elif founded_masses[m] > 1:
            message.warning(m.taxon_name() + ": multiple peptides found for mass " + str(m.mass()) + ".")
    return set_of_new_markers

def add_digestion_status(set_of_markers, set_of_sequences,config_digestion):
    set_of_seqid={m.seqid() for m in set_of_markers if m.sequence() is not None and "Digestion" not in m.field}
    set_of_markers_with_sequences={m for m in set_of_markers if m.sequence() is not None and m.seqid() in set_of_seqid}
    if len(set_of_markers_with_sequences)==0:
        return set_of_markers
    min_length = min(len(m.sequence()) for m in set_of_markers_with_sequences)
    max_length = max(len(m.sequence()) for m in set_of_markers_with_sequences)
    for seqid in set_of_seqid:
        compatible_sequences={s.sequence() for s in set_of_sequences if s.seqid()==seqid}
        if len(compatible_sequences)==0:
            continue
        seq=next(iter(compatible_sequences))
        set_of_peptides=sequences.raw_in_silico_digestion(seq, config_digestion, min_length, max_length)
        for m in set_of_markers:
            if m.seqid()==seqid:
                if m.sequence() in set_of_peptides:
                    m.field["Digestion"]="Yes"
                else:
                    m.field["Digestion"]="No"

def add_length(set_of_markers):
    for m in set_of_markers:
        if m.length() is None :
            if m.sequence() is not None:
                m.field["Length"]=len(m.sequence())
            elif m.begin() is not None and m.end() is not None:
                m.field["Length"]=int(m.end())-int(m.begin())+1

def nomenclature(nom):
    return nom[0] + "-" + str(nom[1]) + "-" + str(nom[2] + nom[1] - 1)

def code_name(m, sequence_dict):
    if m.helical() is not None and m.length() is not None:
        if m.protein() is not None :
            return nomenclature((m.protein(), m.helical(), m.length()))
        else:
            return nomenclature(("COL", m.helical(), m.length()))
    else:
        return None

def add_marker_names(set_of_markers):
    indice=0
    sequence_dict=defaultdict(int)
    for m in set_of_markers:
        if m.sequence() is None or m.sequence() in sequence_dict:
            continue
        indice+=1
        sequence_dict[m.sequence()]=indice
    peptide_to_position=defaultdict(Counter)
    position_to_code= defaultdict(Counter)
    peptide_to_code = defaultdict(Counter)
    for m in set_of_markers:
        if m.sequence() is not None and m.code() is not None:
            peptide_to_code[m.sequence()][m.code()] += 1
        if m.protein() is not None and m.helical() is not None and m.length() is not None :
            if m.sequence() is not None and m.code() is  None:
                peptide_to_position[m.sequence()][(m.protein(), m.helical(),m.length())]+=1
            if m.code() is not None:
                position_to_code[(m.protein(), m.helical(), m.length())][m.code()]+=1
    most_frequent_position_code = {
        a: counter.most_common(1)[0][0]
        for a, counter in position_to_code.items()
    }
    most_frequent_peptide_code = {
        a: counter.most_common(1)[0][0]
        for a, counter in peptide_to_code.items()
    }
    most_frequent_peptide_position = {
        a: counter.most_common(1)[0][0]
        for a, counter in peptide_to_position.items()
    }
    for m in set_of_markers:
        if m.code() is not None :
            continue
        if (m.protein(), m.helical(), m.length()) in position_to_code.keys():
            m.field["Marker"]= most_frequent_position_code[(m.protein(), m.helical(), m.length())]
        elif m.sequence() in peptide_to_code.keys():
            m.field["Marker"] = most_frequent_peptide_code[m.sequence()]
        elif m.sequence() in peptide_to_position.keys():
            m.field["Marker"] = nomenclature(most_frequent_peptide_position[m.sequence()])
        else:
            m.field["Marker"] = code_name(m, sequence_dict)

def add_taxid(set_of_markers, set_of_sequences=None, taxo=None):
    if set_of_sequences is None:
        set_of_sequences=set()
    dict_taxid={m.taxon_name():m.taxid() for m in set_of_markers | set_of_sequences if m.taxid() is not None and m.taxon_name is not None}
    dict_name={value:key for key, value in dict_taxid.items()}
    index=1
    for m in set_of_markers:
        if m.taxid() is not None and m.taxon_name() is not None:
            continue
        if m.taxid() is None and m.taxon_name() is not None:
            if taxo:
                m.field["OX"] = taxonomy.search_taxid_from_taxon_name(m.taxon_name(), taxo)
            else:
                for taxon_name in dict_taxid:
                    if utils.equiv(taxon_name, m.taxon_name()):
                        m.field["OX"]=dict_taxid[taxon_name]
                        continue
        if m.taxid() is not None and m.taxon_name() is None:
            if taxo:
                m.field["OS"]=taxo.name[m.taxid()]
            else:
                m.field["OS"]=dict_name.get(m.taxid())
        if m.taxid() is None and m.taxon_name() is None:
            if m.seqid() is not None:
                seqs=[seq for seq in set_of_sequences if m.seqid()==seq.seqid()]
                if len(seqs)>0:
                    m.field["OS"]=seqs[0].taxon_name()
                    m.field["OX"]=seqs[0].taxid()
        #if m.taxid() is None:
        #    m.field["OX"] = "TxID" + str(index)
        #    index += 1