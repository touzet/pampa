"""
compute_masses.py                       

Everything that is related to the computation of peptide masses  

External dependencies: module mass of pyteomics, copy, re

"""

from pyteomics import mass
import copy
import re
from markers import * 
from utils import * 

OXYPROLINE=15.994915
DEAMIDATION=0.984016

        
def peptide_mass(sequence):
    """ mass of a single peptide, no PTM """
    if not is_aa_sequence(sequence):
        return 0
    return mass.calculate_mass(sequence=sequence, ion_type='M', charge=1)
    # was (with molmass) https://pypi.org/project/molmass/
    # formula = Formula(sequence)
    # mass = formula.isotope.mass

def PTM_mass(PTM_string):
    mass=0
    # proline oxydation
    proline="0"
    if 'P' in ptm_string:
        re_proline=re.compile('[0-9]*P')
        m = re_proline.search(PTM_string)
        proline=m.group().replace('P','')
        mass+= int(proline) * OXYPROLINE 
    # deamidation
    deamidation="0"
    if 'D' in PTM_string:
        re_deamidation=re.compile('[0-9]*D')
        m = re_deamidation.search(PTM_string)
        deamidation=m.group().replace('D','')
        mass+= int(deamidation) * DEAMIDATION
    return mass

def peptide_mass_with_PTM(sequence, PTM_string):
    if not is_aa_sequence(sequence):
        return 0
    return peptide_mass(sequence)+PTM_mass(PTM_string)

def peptide_mass_with_proline(sequence, number_of_prolines):
    if not is_aa_sequence(sequence):
        return 0
    mass = peptide_mass(sequence)
    if number_of_prolines>sequence.count('P'):
        number_of_prolines=sequence.count('P')
    mass+=OXYPROLINE*number_of_prolines
    return mass
         
def peptide_mass_with_proline_range(sequence, min_P, max_P):
    """ compute the list of all masses, in the form of a pair (PTM,mass), for sequence corresponding for a given number of oxyprolines varying from min_P to max_P""" 
    if not is_aa_sequence(sequence):
        print(sequence)
        return []
    if max_P>sequence.count('P'):
        max_P=sequence.count('P')
    mass_list=[]
    mass = peptide_mass(sequence) + min_P*OXYPROLINE
    for i in range (min_P, max_P+1):
        mass_list.append((str(i)+"P", mass))
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
    

def add_PTM_or_masses_to_markers(set_of_markers):
    """
    Compute PTM and masses when there are missing. Existing values are kept.

    Args:  
        a set of markers

    Return:
        a new set_of_markers, composed of the same markers with additional information about PTM and masses

    """
    set_of_new_markers=set()
    set_of_deprecated_markers=set()
    for marker in set_of_markers:
        if len(marker.masses)==0:
            sequence=marker.sequence
            if len(marker.ptm)==0:
                min_P, max_P=proline_range(sequence)
                mass_list= peptide_mass_with_proline_range(sequence, min_P, max_P)
                for ma in mass_list:
                    new_marker=Marker()
                    new_marker=copy.deepcopy(marker)
                    new_marker.ptm=ma[0]
                    new_marker.masses={ma[1]}
                    set_of_new_markers.add(new_marker)
                set_of_deprecated_markers.add(marker)
            else:
                mass=peptide_mass_with_PTM(sequence)
                marker.masses={mass}
    return set_of_markers.union(set_of_new_markers) - set_of_deprecated_markers
    
        

