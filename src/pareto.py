from collections import defaultdict
import numpy as np
from itertools import combinations

from src import neighbour

# a and b are dictionaries
def domination_rank(a, b):
    if a.keys() != b.keys():
        return None
    smaller_than, equal = 0, 0
    for key in a:
        smaller_than += a[key] <= b[key]
        equal += a[key] == b[key]
    return smaller_than if equal<len(a) else -1

# attribut: dictionnaire
# statut: 0=dominé, 1=2-dominé, 2=mieux
def apply_pareto(list_of_attributes):
    n=len(list_of_attributes)
    nb_of_attributes=len(list_of_attributes[0])
    statut = [-2]*n
    dr_matrix=[]
    for i in range(n):
        mylist=[domination_rank(list_of_attributes[i], list_of_attributes[j]) for j in range(n)]
      # est-ce que i est dominé ?
        if nb_of_attributes in mylist:
            statut[i]=-1
            mylist=[-2]*n
        dr_matrix.append(mylist)
    for j in range(n):
        if statut[j]==-1:
            continue
        statut[j]=[dr_matrix[i][j] for i in range(n)].count(2)
    return statut

def pareto_comparison(marker1, marker2, list_of_attributes):
    greater1, greater2 = 0,0
    for attribute in list_of_attributes:
        if attribute in marker1.field and attribute in marker2.field:
            if marker1.field[attribute] > marker2.field[attribute]:
                greater1 += 1
            elif marker1.field[attribute] < marker2.field[attribute]:
                greater2 += 1
    return greater1, greater2

# all markers are assumed to have the same code
def pareto_score(m, set_of_markers, list_of_attributes):
    score = len(list_of_attributes)
    for m2 in set_of_markers:
        if m2.code() == m.code() :
            greater, equal = pareto_comparison(m, m2, list_of_attributes)
            if greater == 0 and equal < len(list_of_attributes):
                score = 0
            elif greater+equal < score :
                score = greater+equal
    return score

# all markers are for the same code and taxon
def aggregative_score(set_of_markers,  list_of_attributes):
    for (attribute, weight) in list_of_attributes:
        values={float(m.field[attribute]) for m in set_of_markers if attribute in m.field}
        if len(values)>0:
            max_value=max(values)
        if len(values)==0 or max_value==0 :
            continue
        for m in set_of_markers:
            if attribute in m.field :
                to_add=m.field[attribute]/max_value*weight
                m.field["score"]+= to_add


# all markers are assumed to have the same code
def add_pareto_score(set_of_markers):
    if len(set_of_markers)==0:
        return set()
    list_of_attributes=[ "Subst", "Mutability", "Paired", "Intensity", "Gamma", "Spectra", "Proh", "CP", "Exists"]
    for m in set_of_markers:
        m.field["Pareto"] =  pareto_score(m, set_of_markers, list_of_attributes)
        if m.field["Neighbour"] == 1:
            m.field["Pareto"] += 10
        elif m.field["Subst"] == 0.0:
            m.field["Pareto"] = 0
    #return set_of_markers

def pareto(set_of_markers,  list_of_attributes):
    attributes=[a[0] for a in list_of_attributes]
    for m in set_of_markers:
        m.field["Pareto"] = 1
    for m,m2 in combinations(set_of_markers, 2) :
        if m2.code()!=m.code() or m2.sequence()==m.sequence():
            continue
        greater_m, greater_m2 = pareto_comparison(m, m2, attributes)
        if greater_m == 0 and greater_m2>0:
            m.field["Pareto"] = 0
        if greater_m>0 and greater_m2==0:
            m2.field["Pareto"] = 0
    for m, m2 in combinations(set_of_markers, 2):
        if m.code()==m2.code() and m.sequence()==m2.sequence():
            m.field["Pareto"]=max(m.field["Pareto"],m2.field["Pareto"])
            m2.field["Pareto"]=max(m.field["Pareto"],m2.field["Pareto"])


def add_aggregative_score(set_of_markers):
    if len(set_of_markers)==0:
        return
    for m in set_of_markers:
        m.field["score"]=0.0
    list_of_attributes=[("Neighbour",5),("Subst",1), ("Mutability",0), ("Paired", 1), ("Intensity",2), ("Gamma", 1),("Spectra",1), ("Proh",1), ("CP",1), ("Exists",2), ("PTM_adequacy",1),("Delta_mz",1)]
    set_of_codes={m.code() for m in set_of_markers}
    set_of_good_markers=set()
    for code in set_of_codes:
        set_of_selected_markers={m for m in set_of_markers if m.code()==code}
        aggregative_score(set_of_selected_markers,  list_of_attributes)
        set_of_good_markers.update(select_good_scoring_candidates(set_of_selected_markers))
    return neighbour.remove_duplicates(set_of_good_markers)


# all markers are assumed to have the same code
# deprecated
"""def add_pareto_score(set_of_markers):
    if len(set_of_markers)==0:
        return []
    list_of_markers=list(set_of_markers)
    list_of_attributes=[{"s":m.field["Subst"], "m":m.field["Mutability"], "i":m.field["Intensity"], "g":m.field["Gamma"]} for m in list_of_markers]
    pareto_statut=apply_pareto(list_of_attributes)
    for i, m in enumerate(list_of_markers):
        m.field["Pareto"]=pareto_statut[i]
    return list_of_markers
"""

def compute_threshold(scores):
    threshold=scores[0]*0.6
    return threshold

# all markers are supposed to belong to the same species and codes
def select_good_scoring_candidates(set_of_markers):
    list_of_scores=[m.field["score"] for m in set_of_markers]
    list_of_scores.sort(reverse=True)
    threshold=compute_threshold(list_of_scores)
    for m in set_of_markers:
        m.field["keep"]=m.field["score"]>=threshold
        m.field["rk"]=list_of_scores.index(m.field["score"])
    return {m for m in set_of_markers if m.field["keep"]}
