"""
compute_masses.py                       

Everything that is related to the computation of peptide masses  

"""

from pyteomics import mass
import copy
import re


from src import utils as ut
from src import collagen as collagen
from src import message

HYDROXYPROLINE=15.994915 # proline (P)
DEAMIDATION=0.984016 # asparagine (N) and glutamine (Q)
CARBOXYLATION= 43.009 #
PHOSPHORYLATION= 79.9663 # serine (S), threonine (T), or tyrosine (Y)

def counting_matching_characters(sequence, set_of_PTM):
    return sum(char in set(set_of_PTM) for char in sequence)

def peptide_mass(sequence):
    """ mass of a single peptide, no PTM """
    if not ut.is_aa_sequence(sequence):
        return 0
    return mass.calculate_mass(sequence=sequence, ion_type='M', charge=1)
    # was (with molmass) https://pypi.org/project/molmass/
    # formula = Formula(sequence)
    # mass = formula.isotope.mass

def number_of_H(PTM_string):
    if PTM_string is None:
        return 0
    if 'H' in PTM_string:
        re_proline=re.compile('[0-9]*H')
        m = re_proline.search(PTM_string)
        proline=int(m.group().replace('H',''))
    else:
        proline=0
    return proline

def number_of_D(PTM_string):
    if PTM_string is None:
        return 0
    if 'D' in PTM_string:
        re_deamidation=re.compile('[0-9]*D')
        m = re_deamidation.search(PTM_string)
        deamidation=int(m.group().replace('D',''))
    else:
        deamidation=0
    return deamidation

def number_of_C(PTM_string):
    if PTM_string is None:
        return 0
    if 'C' in PTM_string:
        re_carboxylation=re.compile('[0-9]*C')
        m = re_carboxylation.search(PTM_string)
        carboxylation=int(m.group().replace('C',''))
    else:
        carboxylation=0
    return carboxylation

def PTM_mass(PTM_string):
    m=0
    # proline oxydation
    m+= number_of_H(PTM_string) * HYDROXYPROLINE
    # deamidation
    m+= number_of_D(PTM_string) * DEAMIDATION
    #carboxylation
    m+= number_of_C(PTM_string) * CARBOXYLATION
    return m

def peptide_mass_with_PTM(sequence, PTM_string):
    if not ut.is_aa_sequence(sequence) or PTM_string is None:
        return 0
    return peptide_mass(sequence)+PTM_mass(PTM_string)

#deprecated ?
def peptide_mass_with_proline(sequence, number_of_prolines):
    if not ut.is_aa_sequence(sequence):
        return 0
    mass = peptide_mass(sequence)
    if number_of_prolines>sequence.count('P'):
        number_of_prolines=sequence.count('P')
    mass+=HYDROXYPROLINE*number_of_prolines
    return mass
         
def peptide_mass_with_proline_range(sequence, min_P, max_P):
    """ compute the list of all masses, in the form of a pair (PTM,mass), for sequence corresponding for a given number of hydroxyprolines varying from min_P to max_P"""
    if not ut.is_aa_sequence(sequence):
        return []
    if (min_P, max_P)==(-1,-1):
        return [(" ", peptide_mass(sequence))]
    if max_P>sequence.count('P'):
        max_P=sequence.count('P')
    mass_list=[]
    mass = peptide_mass(sequence) + min_P*HYDROXYPROLINE
    for i in range (min_P, max_P+1):
        mass_list.append((str(i)+"H", mass))
        mass+=HYDROXYPROLINE
    return mass_list

def proline_range(sequence, more_hydroxyprolines=False):
    """ Estimates the minimal and maximal number of hydroxyprolines in a collagen sequence
        Returns (-1, -1) if the sequence does not follow the G-X-Y pattern"""
    strong_P, weak_P= collagen.P_pattern(sequence)
    if (strong_P, weak_P)==(-1,-1):
        return (-1, -1)
    if weak_P <4:
        min_P, max_P =strong_P, strong_P
    else:
        min_P, max_P = strong_P, strong_P+1
    if more_hydroxyprolines : # should be integrated to proline_range
        if min_P>0:
            min_P-=1
        if max_P<sequence.count('P'):
            max_P+=1
    return min_P, max_P
    

def update_PTM(new_sequence, source_sequence, set_of_source_PTM, deamidation):
    new_P_pattern = collagen.P_pattern(new_sequence)
    old_P_pattern = collagen.P_pattern(source_sequence)
    set_of_new_PTM= {str(new_P_pattern[0])+'H'}
    deamidated = 'Q' in new_sequence or 'N' in new_sequence
    for source_PTM in set_of_source_PTM:
        nb_H = number_of_H(source_PTM)
        if new_P_pattern[0] > old_P_pattern[0]:
             nb_H += 1
        elif new_P_pattern[0] < old_P_pattern[0]:
            nb_H =  max(0, nb_H - 1)
        new_PTM=str(nb_H)+'H'
        set_of_new_PTM.add(new_PTM)
        if deamidation and deamidated :
            set_of_new_PTM.add(new_PTM+'1D')
    return set_of_new_PTM

def compute_masses(m, more_hydroxyprolines, deamidation):
    sequence=m.sequence()
    if m.PTM() :
        m.field["Mass"]=peptide_mass_with_PTM(sequence, m.PTM())
        return {m}
    to_add = set()
    min_P, max_P=proline_range(sequence, more_hydroxyprolines)
    if (min_P, max_P)==(-1,-1):
        if "magic_number" in m.field:
            message.warning("Marker " + m.field[
                "magic_number"] + ": " + sequence + " is not a collagen sequence. Unable to compute hydroxyprolines")
        else:
            message.warning(sequence + " is not a collagen sequence. Unable to compute hydroxyprolines")
    mass_list = peptide_mass_with_proline_range(sequence, min_P, max_P)
    if deamidation and ('Q' in sequence or 'N' in sequence):
        mass_list_deamidation = [(ma[0] + "1D", ma[1] + DEAMIDATION) for ma in mass_list]
        mass_list = mass_list + mass_list_deamidation
    for ma in mass_list:
        new_marker = copy.deepcopy(m)
        new_marker.field["PTM"] = ma[0]
        new_marker.field["Mass"] = ma[1]
        to_add.add(new_marker)
    m.field["ToRemove"] = "yes"
    return to_add

def add_PTM_or_masses_to_markers(set_of_markers, more_hydroxyprolines=False, deamidation=False):
    """
    Compute PTM and masses when there are missing. Existing values are kept.
    Args:  
        a set of markers
    Return:
        a new set_of_markers, composed of the same markers with additional information about PTM and masses
    """
    to_add=set()
    for m in set_of_markers:
        if m.mass() is not None:
            to_add.add(m)
        if m.sequence() is None:
            continue
        if m.mass() is None :
            to_add.update(compute_masses(m, more_hydroxyprolines, deamidation))
    return to_add
    
     
def compatible_mass(sequence, PTM, mass, resolution):
    if mass is None:
        return True, PTM, None
    if PTM is not None:
        return ut.matching_masses(peptide_mass_with_PTM(sequence, PTM), mass, resolution), PTM, peptide_mass_with_PTM(sequence, PTM)
    # PTM is None:
    min_P, max_P=proline_range(sequence)
    if min_P>0: # test à intégrer à proline_range, ici et ailleurs
        min_P-=1
    if max_P<sequence.count('P'):
        max_P+=1
    mass_list = peptide_mass_with_proline_range(sequence, min_P, max_P)
    if ('Q' in sequence or 'N' in sequence):
        mass_list_deamidation=[(ma[0]+"1D", ma[1]+DEAMIDATION)
        for ma in mass_list]
        mass_list = mass_list + mass_list_deamidation
    for ma in mass_list:
        if ut.matching_masses(ma[1],mass, resolution):
            return True, ma[0], ma[1]
    return False, None, None
               
def add_deamidation(set_of_markers, set_of_codes=None):
    if set_of_codes is None :
        set_of_authorized_codes={m.code() for m in set_of_markers}
    elif len(set_of_codes)==0:
        return set()
    else:
        set_of_authorized_codes=set_of_codes
    set_of_new_markers=set()
    nothing_to_do={(m.sequence(), number_of_H(m.PTM()), m.taxid(), m.code()) for m in set_of_markers if m.sequence() is None or m.code() not in set_of_authorized_codes or (m.PTM() is not None and 'D' in m.PTM())}
    for m in set_of_markers:
        if  (m.sequence(), number_of_H(m.PTM()), m.taxid(), m.code()) not in nothing_to_do and (m.sequence() is not None and ('Q' in m.sequence() or 'N' in m.sequence())):
            new_marker=copy.deepcopy(m)
            if m.PTM() is None:
                new_marker.field["PTM"]='1D'
            else:
                new_marker.field["PTM"]=m.PTM()+'1D'
            new_marker.field["Mass"]=new_marker.mass()+DEAMIDATION
            new_marker.field["Comment"]=new_marker.comment()+ " + deamidation"
            set_of_new_markers.add(new_marker)
    return set_of_new_markers

