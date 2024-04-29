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


OXYPROLINE=15.994915 # proline (P)
DEAMIDATION=0.984016 # asparagine (N) and glutamine (Q)
CARBOXYLATION= 43.009 #
PHOSPHORYLATION= 79.9663 # serine (S), threonine (T), or tyrosine (Y)


def couting_matching_characters(sequence, set_of_PTM):
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
    if 'O' in PTM_string:
        re_proline=re.compile('[0-9]*O')
        m = re_proline.search(PTM_string)
        proline=m.group().replace('O','')
        mass+= int(proline) * OXYPROLINE 
    # deamidation
    deamidation="0"
    if 'D' in PTM_string:
        re_deamidation=re.compile('[0-9]*D')
        m = re_deamidation.search(PTM_string)
        deamidation=m.group().replace('D','')
        mass+= int(deamidation) * DEAMIDATION
    carboxylation="0"
    if 'C' in PTM_string:
        re_carboxylation=re.compile('[0-9]*C')
        m = re_carboxylation.search(PTM_string)
        carboxylation=m.group().replace('C','')
        mass+= int(carboxylation) * CARBOXYLATION
    return mass

def peptide_mass_with_PTM(sequence, PTM_string):
    if not ut.is_aa_sequence(sequence):
        return 0
    return peptide_mass(sequence)+PTM_mass(PTM_string)

def peptide_mass_with_proline(sequence, number_of_prolines):
    if not ut.is_aa_sequence(sequence):
        return 0
    mass = peptide_mass(sequence)
    if number_of_prolines>sequence.count('O'):
        number_of_prolines=sequence.count('O')
    mass+=OXYPROLINE*number_of_prolines
    return mass
         
def peptide_mass_with_proline_range(sequence, min_P, max_P):
    """ compute the list of all masses, in the form of a pair (PTM,mass), for sequence corresponding for a given number of oxyprolines varying from min_P to max_P""" 
    if not ut.is_aa_sequence(sequence):
        return []
    if max_P>sequence.count('P'):
        max_P=sequence.count('P')
    mass_list=[]
    mass = peptide_mass(sequence) + min_P*OXYPROLINE
    for i in range (min_P, max_P+1):
        mass_list.append((str(i)+"O", mass))
        mass+=OXYPROLINE
    return mass_list
                
def proline_range(sequence):
    """ Estimates the minimal and maximal number of oxyprolines in a sequence """
    #total number of P
    proline=sequence.count('P')
    # number of P involved in the pattern "G.P"
    period=re.compile('G\wP')
    period_proline=len(period.findall(sequence))
    if proline - period_proline<3:
        return period_proline, period_proline
    else:
        return period_proline, period_proline+1
    

def add_PTM_or_masses_to_markers(set_of_markers, flexible=False):
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
        sequence=marker.sequence
        if marker.code==None:
            if sequence not in list_of_codes:
                    list_of_codes.append(sequence)
            marker.code='M'+str(list_of_codes.index(sequence))
        if marker.mass==None:
            if marker.ptm==None:
                min_P, max_P=proline_range(sequence)
                if flexible:
                    if min_P>0:
                        min_P-=1
                    if max_P<sequence.count('P'):
                        max_P+=1
                mass_list= peptide_mass_with_proline_range(sequence, min_P, max_P)
                for ma in mass_list:
                    new_marker=markers.Marker()
                    new_marker=copy.deepcopy(marker)
                    new_marker.ptm=ma[0]
                    new_marker.mass=ma[1]
                    set_of_new_markers.add(new_marker)
                set_of_deprecated_markers.add(marker)
            else:
                mass=peptide_mass_with_PTM(sequence,marker.ptm)
                marker.mass=mass
                    
    return set_of_markers.union(set_of_new_markers) - set_of_deprecated_markers
    
     
        
    
