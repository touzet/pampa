import itertools
import math
import os
import pandas as pd
import re
from collections import defaultdict
from collections import Counter

from src import  collagen, taxonomy, message, utils

def fetch_aa(taxid,  set_of_sequences, position):
    for seq in set_of_sequences:
        if seq.taxid() == taxid :
            if position < len(seq.sequence()) :
                return seq.sequence()[position]
            else:
                return 'X'


# leaf-to-root step
def fitch_rec_leaf_to_root(taxid, set_of_sequences, position, taxo, list_of_states_dict):
    if taxo.is_leaf(taxid):
        list_of_states_dict[taxid] = [fetch_aa(taxid, set_of_sequences, position)]
    else:
        for c in taxo.children[taxid]:
            fitch_rec_leaf_to_root(c, set_of_sequences, position, taxo, list_of_states_dict)
        list_of_aa_labels=[list_of_states_dict[t] for t in taxo.children[taxid]]
        list_of_aa_label_sets=[set(l) for l in list_of_aa_labels]
        intersect=set.intersection(*list_of_aa_label_sets)
        if len(intersect)>0:
            list_of_states_dict[taxid] = list(intersect)
        else:
            list_of_states_dict[taxid] = sum(list_of_aa_labels,[])
    return

def fitch_root(taxo, fitch_dict, list_of_states_dict):
    for taxid in taxo.root:
        aa_counts = Counter(list_of_states_dict[taxid])
        max_count = aa_counts.most_common(1)[0][1]
        fitch_dict[taxid] = [aa for aa, count in aa_counts.items() if count == max_count]

# root_to_leaf step
def fitch_rec_root_to_leaves(taxid, list_of_states_dict, taxo, fitch_dict):
    if taxo.is_leaf(taxid):
        return
    else:
        for c in taxo.children[taxid]:
            inters=set(list_of_states_dict[c]).intersection(set(fitch_dict[taxid]))
            if len(inters)==0:
                fitch_dict[c]=list_of_states_dict[c]
            else:
                fitch_dict[c]=list(inters)
            fitch_rec_root_to_leaves(c, list_of_states_dict, taxo, fitch_dict)

def print_substitution_tree(tree_list):
    tree_set=set(tree_list)
    dict_substitution={taxid:(aa, aa2) for (taxid,aa,aa2) in tree_set if aa!=aa2}

def generate_all_substitutions(parent, list_of_child_sets):
    list_of_substitutions = []
    for s in list_of_child_sets:
        inter=set(parent) & set(s)
        if len(inter)==0:
            list_of_substitutions.extend([(min(a, b), max(a,b)) for a in parent for b in s])
        else:
            list_of_substitutions.extend([(a, a) for a in inter])
    return list_of_substitutions

def root_substitutions(taxid, fitch_dict):
    L=fitch_dict[taxid]
    if len(L)==1:
        return []
    else :
        return [(min(x,y), max(x,y)) for i, x in enumerate(L) for y in L[i+1:]]

# compute the list of all substitutions
def substitutions(taxid, taxo, fitch_dict):
    if taxo.is_leaf(taxid):
        return []
    list_of_sets=[fitch_dict[c] for c in taxo.children[taxid]]
    substitution_list= generate_all_substitutions(fitch_dict[taxid], list_of_sets)
    for c in taxo.children[taxid]:
        substitution_list.extend(substitutions(c, taxo, fitch_dict))
    return substitution_list


# compute the list of all substitutions
def old_substitutions(taxid, taxo, fitch_dict):
    if taxo.is_leaf(taxid):
        parent=taxo.parent[taxid]
        if fitch_dict[taxid] <= fitch_dict[parent]:
            aa=next(iter(fitch_dict[taxid]))
            return[(aa,aa)]
        else:
            return []
    substitution_list=[]
    if len(taxo.children[taxid])>1:
        list_of_sets=[fitch_dict[c] for c in taxo.children[taxid]]
        substitution_list= generate_all_substitutions(fitch_dict[taxid], list_of_sets)
    for c in taxo.children[taxid]:
        substitution_list.extend(substitutions(c, taxo, fitch_dict))
    return substitution_list

def fitch(set_of_sequences, set_of_positions, taxoB):
    taxo_R=taxonomy.reduce_tree(taxoB)
    # A METTRE AVANT
    #set_of_lengths={len(seq.sequence()) for seq in set_of_sequences}
    #if len(set_of_lengths) != 1:
    #    message.escape("Error with input FASTA sequences: all sequences should have the same length.")
    #set_of_taxid={seq.taxid() for seq in set_of_sequences}
    #taxoB, _ = taxo.intersection(set_of_taxid)
    all_substitutions=[]
    for position in set_of_positions:
        set_of_aa={seq.sequence()[position] for seq in set_of_sequences if position<len(seq.sequence()) and seq.sequence()[position]!='X'}
        if len(set_of_aa) > 1:
            list_of_states_dict = {}
            for taxid in taxoB.root:
                fitch_rec_leaf_to_root(taxid, set_of_sequences, position, taxo_R, list_of_states_dict)
                fitch_dict={taxid:list_of_states_dict[taxid]}
                fitch_root(taxo_R, fitch_dict, list_of_states_dict)
                #root_states=fitch_root(taxo_R, fitch_dict, list_of_states_dict)
                #for state in root_states:
                #    fitch_dict[taxid]=state
                fitch_rec_root_to_leaves(taxid, list_of_states_dict, taxo_R, fitch_dict)
                s=substitutions(taxid, taxo_R, fitch_dict)
                s.extend(root_substitutions(taxid, fitch_dict))
                all_substitutions.extend(s)
        else:
            aa=next(iter(set_of_aa))
            s=[(aa,aa)]*len(set_of_sequences) # to check !
            all_substitutions.extend(s)
    all_substitutions=[aa_pair for aa_pair in all_substitutions if aa_pair[0]!='X' and aa_pair[1]!='X']
    aa_range = ['C', 'S', 'T', 'A', 'G', 'P', 'D', 'E', 'Q', 'N', 'H', 'R', 'K', 'M', 'I', 'L', 'V', 'Y', 'F', 'W']
    df = pd.DataFrame(0.0, index=aa_range, columns=aa_range)
    for (k,l) in all_substitutions:
        df.at[k, l] += 1
        if k!= l:
            df.at[l, k] += 1
    return df

def compute_substitution_matrix(list_of_sequences, set_of_positions, taxo):
    set_of_taxid = {seq.taxid() for seq in list_of_sequences}
    taxoB, _ = taxo.intersection(set_of_taxid)
    matrix = fitch(list_of_sequences, set_of_positions, taxoB)
    number_of_occurrences= {}
    number_of_substitutions = {}
    for aa in matrix.columns:
        number_of_occurrences[aa]=matrix[aa].sum()
        number_of_substitutions[aa]= number_of_occurrences[aa] - matrix.at[aa,aa]
        if number_of_occurrences[aa]>0:
            # mutability
            matrix.at[aa,aa]= round(number_of_substitutions[aa]*100/ number_of_occurrences[aa], 1)
        else:
            matrix.at[aa,aa] = -1
    sum_of_aa=sum(number_of_occurrences.values())
    sum_of_substitutions = sum(matrix.at[aa,aa2]/2 for aa in matrix.columns for aa2 in matrix.columns if aa!=aa2)
    for aa in matrix.columns:
        for aa2 in matrix.columns:
            if matrix.at[aa,aa2]>0 and aa!=aa2:
                matrix.at[aa,aa2]=round(matrix.at[aa,aa2]*100 / (number_of_substitutions[aa]+number_of_occurrences[aa2] - matrix.at[aa,aa2]),1)
    return matrix

def print_matrix(matrix):
    print(matrix)

def same_length(set_of_strings):
    return len({len(s) for s in set_of_strings}) == 1

def same_phase(set_of_strings):
    set_of_phases={collagen.phasing_GXY_pattern(s) for s in set_of_strings}
    return len(set_of_phases)==1 and -1 not in set_of_phases

def compute_collagen_matrices_from_set_of_sequences(set_of_sequences, taxo, file_name="custom"):
    matrices={}
    for prot in ["COL1A1", "COL1A2", "COL1A3"]:
        matrices[prot]={}
        list_of_sequences=[s for s in set_of_sequences if s.protein()==prot]
        if len(list_of_sequences)==0:
            message.warning("No Fasta sequence found for "+prot)
            continue
        # ajouter un test sur la longueur des séquences ?
        X_positions, Y_positions=collagen.compute_X_Y_positions(list_of_sequences[0].sequence())
        matrices[prot]['X']=compute_substitution_matrix(list_of_sequences, X_positions, taxo)
        matrices[prot]['Y']=compute_substitution_matrix(list_of_sequences, Y_positions, taxo)
    return matrices, create_matrix_files(matrices, ".", file_name)

def read_scoring_matrix(path):
    df = pd.read_csv(path, index_col=0)
    return df

def find_genes(list_of_files):
    """
    Find GN such that files
      prefix_GN_X.csv and prefix_GN_Y.csv both exist in dir_path
    """
    GN_X=set()
    GN_Y=set()
    for path in list_of_files:
        filename=os.path.basename(path)
        matrix_dir=os.path.dirname(path)
        parts = filename.split("_")
        prefix="_".join(parts[0:-2])
        if filename.endswith("_X.csv"):
            GN_X.add(parts[-2])  # extract GN
        if filename.endswith("_Y.csv"):
            GN_Y.add(parts[-2]) # extract GN
    return GN_X & GN_Y, matrix_dir, prefix

def load_collagen_matrices(list_of_files):
    matrices = {}
    set_of_genes, matrix_dir, matrix_file= find_genes(list_of_files)
    for gene in set_of_genes:
        matrices[gene] = {
            'Y': read_scoring_matrix(os.path.join(matrix_dir, f'{matrix_file}_{gene}_Y.csv')),
            'X': read_scoring_matrix(os.path.join(matrix_dir, f'{matrix_file}_{gene}_X.csv'))
        }
    return matrices

# tous les marqueurs sont sur la même protéine
def compute_collagen_matrix_from_set_of_markers(set_of_markers):
    set_of_taxid = {m.taxid() for m in set_of_markers}
    code_to_taxid, code_to_sequence = {}, {}
    for m in set_of_markers:
        utils.update_dictoset(code_to_sequence, m.code(), {m.sequence()})
        utils.update_dictoset(code_to_taxid, m.code(), {m.taxid()})
    set_of_codes={m.code() for m in set_of_markers if same_length(code_to_sequence[m.code()]) and same_phase(code_to_sequence[m.code()]) and code_to_taxid[m.code()] == set_of_taxid}
    mydict = {}
    for m in set_of_markers:
        if m.code() in set_of_codes:
            mydict[(m.taxid(), m.code())] = m.sequence()
    list_of_codes = list(set_of_codes)
    list_of_sequences = []
    for taxid in set_of_taxid:
        seq = ''.join([mydict[(taxid, code)] for code in list_of_codes])
        list_of_sequences.append(seq)
    X_all_positions, Y_all_positions = set(), set()
    taxid = set_of_taxid.pop()
    code_length=0
    for code in list_of_codes:
        X_positions, Y_positions = collagen.compute_X_Y_positions(mydict[(taxid, code)])
        X_all_positions.update({p + code_length for p in X_positions})
        Y_all_positions.update({p + code_length for p in Y_positions})
        code_length+=len(mydict[taxid, code])
    return create_scoring_matrix(list_of_sequences, X_all_positions), create_scoring_matrix(list_of_sequences, Y_all_positions)

def compute_collagen_matrices_from_set_of_markers(set_of_markers):
    matrices={}
    for prot in ["COL1A1","COL1A2"]:
        matrices[prot] = {}
        set_of_m= {m for m in set_of_markers if m.protein()==prot}
        matrices[prot]["X"], matrices[prot]["Y"]=compute_collagen_matrix_from_set_of_markers(set_of_m)
    create_matrix_files(matrices, ".")
    return matrices

def create_matrix_files(matrices, path, seed_name):
    list_of_files = []
    for prot in matrices:
        for c in matrices[prot]:
            file_name="matrix_" + seed_name + "_" + prot + "_" + c + ".csv"
            matrices[prot][c].to_csv(os.path.join(path, file_name))
            list_of_files.append(file_name)
    return list_of_files
