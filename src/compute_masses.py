"""
compute_masses.py                       

Everything that is related to the computation of peptide masses  

"""

from pyteomics import mass
import copy
import re
import sys

from src import markers 
from src import utils as ut


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


def PTM_mass(PTM_string):
    mass=0
    # proline oxydation
    proline="0"
    if 'H' in PTM_string:
        re_proline=re.compile('[0-9]*H')
        m = re_proline.search(PTM_string)
        proline=m.group().replace('H','')
        mass+= int(proline) * HYDROXYPROLINE 
    # deamidation
    deamidation="0"
    if 'D' in PTM_string:
        re_deamidation=re.compile('[0-9]*D')
        m = re_deamidation.search(PTM_string)
        deamidation=m.group().replace('D','')
        mass+= int(deamidation) * DEAMIDATION
    #carboxylation
    carboxylation="0"
    if 'C' in PTM_string:
        re_carboxylation=re.compile('[0-9]*C')
        m = re_carboxylation.search(PTM_string)
        carboxylation=m.group().replace('C','')
        mass+= int(carboxylation) * CARBOXYLATION
    return mass

def peptide_mass_with_PTM(sequence, PTM_string):
    if not ut.is_aa_sequence(sequence) or PTM_string is None:
        return 0
    return peptide_mass(sequence)+PTM_mass(PTM_string)

def peptide_mass_with_proline(sequence, number_of_prolines):
    if not ut.is_aa_sequence(sequence):
        return 0
    mass = peptide_mass(sequence)
    if number_of_prolines>sequence.count('P'):
        number_of_prolines=sequence.count('P')
    mass+=HYDROXYPROLINE*number_of_prolines
    return mass
         
def peptide_mass_with_proline_range(sequence, min_P, max_P):
    """ compute the list of all masses, in the form of a pair (PTM,mass), for sequence corresponding for a given number of oxyprolines varying from min_P to max_P""" 
    if not ut.is_aa_sequence(sequence):
        return []
    if max_P>sequence.count('P'):
        max_P=sequence.count('P')
    mass_list=[]
    mass = peptide_mass(sequence) + min_P*HYDROXYPROLINE
    for i in range (min_P, max_P+1):
        mass_list.append((str(i)+"H", mass))
        mass+=HYDROXYPROLINE
    return mass_list

# warning: may be incorrect if the peptide contains a motif of the form "GxGP"
def proline_range(sequence):
    """ Estimates the minimal and maximal number of oxyprolines in a sequence """
    #total number of P
    number_of_prolines=sequence.count('P')
    # number of P involved in the pattern "G.P"
    period=re.compile('G\wP')
    period_proline=len(period.findall(sequence))
    if number_of_prolines - period_proline<4:
        return period_proline, period_proline
    else:
        return period_proline, period_proline+1
    

def add_PTM_or_masses_to_markers(set_of_markers, more_hydroxyprolines=False, deamidation=False):
    """
    Compute PTM and masses when there are missing. Existing values are kept.

    Args:  
        a set of markers

    Return:
        a new set_of_markers, composed of the same markers with additional information about PTM and masses

    """
    set_of_new_markers=set()
    set_of_deprecated_markers=set()
    list_of_codes=[]
    for marker in set_of_markers:
        sequence=marker.sequence()
        if sequence is None:
            set_of_new_markers.add(marker)
            continue
        if marker.code() is None:
            if sequence not in list_of_codes:
                    list_of_codes.append(sequence)
            marker.field["Marker"]='m'+str(list_of_codes.index(sequence))
        if marker.mass() is None:
            if marker.PTM() is None:
                min_P, max_P=proline_range(sequence)
                if more_hydroxyprolines:
                    if min_P>0:
                        min_P-=1
                    if max_P<sequence.count('P'):
                        max_P+=1
                mass_list = peptide_mass_with_proline_range(sequence, min_P, max_P)
                if deamidation and ('Q' in sequence or 'N' in sequence):
                    mass_list_deamidation=[(ma[0]+"1D", ma[1]+DEAMIDATION) for ma in mass_list]
                    mass_list = mass_list + mass_list_deamidation
                for ma in mass_list:
                    new_marker=markers.Marker()
                    new_marker=copy.deepcopy(marker)
                    new_marker.field["PTM"]=ma[0]
                    new_marker.field["Mass"]=ma[1]
                    set_of_new_markers.add(new_marker)
                set_of_deprecated_markers.add(marker)
            else:
                mass=peptide_mass_with_PTM(sequence,marker.PTM())
                marker.field["Mass"]=mass
                    
    return set_of_markers.union(set_of_new_markers) - set_of_deprecated_markers
    
     
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
               
def add_deamidation(set_of_markers, set_of_codes):
    if len(set_of_codes)==0:
        set_of_codes={m.code() for m in set_of_markers}
    set_of_new_markers=set()
    for m in set_of_markers:
            if  m.code() in set_of_codes and (m.PTM() is None or 'D' not in m.PTM()) and (m.sequence() is None or ('Q' in m.sequence() or 'N' in m.sequence())):
                new_marker=copy.deepcopy(m)
                if m.PTM() is None:
                    new_marker.field["PTM"]='1D'
                else:
                    new_marker.field["PTM"]=m.PTM()+'1D'
                new_marker.field["Mass"]=new_marker.mass()+DEAMIDATION
                new_marker.field["Comment"]=new_marker.comment()+ " + deamidation"
                set_of_new_markers.add(new_marker)
    return set_of_new_markers
