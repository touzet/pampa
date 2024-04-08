import copy
import os
from os import listdir
from os.path import join

from src import utils 
from src import mass_spectrum as ms


# input:
# list of markers, represented by a sorted list of masses (in increasing order)
# one spectrum
# indique pour chaque marqueur s'il y a un pic correspondant
# resultat: liste avec le nom du spectre si match 
def compare_set_of_markers_for_one_spectrum(mass_list, spectrum, resolution):
    spectrum.sort()  
    peak_list=[set() for m in mass_list] 
    current_j = 0
    for peak in spectrum:
        for j in range(current_j,len(mass_list)):
            difference = float(peak) - float(mass_list[j][0])
            if utils.matching_masses(peak, mass_list[j][0], resolution): # on a un match
                peak_list[j]={spectrum.name} 
            elif difference<0: # la masse du peptide est supérieure à la masse du pic  
                current_j = j # on reprend au même endroit où le pic d'avant s'est arrêté
                break # on ne va pas plus loin dans la liste des masses des peptides et on passe au pic suivant
    return peak_list



def compare_markers_with_spectra(set_of_markers, list_of_spectra, resolution):
    mass_list=[(m.mass,m) for m in set_of_markers]
    mass_list.sort(key=lambda x:x[0])
    peak_list_all_spectra=[set() for m in mass_list]
    for spectrum in list_of_spectra:
        peak_list= compare_set_of_markers_for_one_spectrum(mass_list, spectrum, resolution)
        peak_list_all_spectra=[a.union(b) for (a,b) in zip(peak_list_all_spectra, peak_list)]
    return {m[1]:b for (m,b) in zip(mass_list, peak_list_all_spectra)}

    
# les marqueurs sont supposés être de la même espèce que les spectres
def filter_set_of_markers(set_of_markers, list_of_spectra, resolution, min_nb_spectra=0):
    dict_markers= compare_markers_with_spectra(set_of_markers, list_of_spectra, resolution)
    set_of_confirmed_markers=set()
    for m in set_of_markers:
        found_spectra=dict_markers[m]
        if len(found_spectra)>=min_nb_spectra:
            m.comment=m.comment+" - "+ str(len(found_spectra))+"/"+str(len(list_of_spectra))+" spectra :" + str(found_spectra)
            set_of_confirmed_markers.add(m)
    return set_of_confirmed_markers           

def filter_set_of_markers_multispecies(set_of_markers, spectra_dir, resolution, min_ratio_spectra=0, min_ratio_species=0):
    dict_of_codes_ptm={(m.code,m.ptm):0 for m in set_of_markers}
    set_of_taxa={m.taxon_name for m in set_of_markers}
    dict_all_species={}
    dict_spectra={}
    missing_species=set()
    for taxon in set_of_taxa:
        taxon_dir= spectra_dir+"/"+(taxon.title()).replace(" ","") # TO REPLACE
        if not(os.path.isdir(taxon_dir)):
            missing_species.add(taxon)
            continue
        list_of_spectra=[]
        for f in listdir(taxon_dir):
            file_name= join(taxon_dir, f)
            spectrum=mass_spectrum.parser(file_name,f)
            list_of_spectra.append(spectrum)
        dict_spectra[taxon]=len(list_of_spectra)
        set_of_current_markers={m for m in set_of_markers if m.taxon_name==taxon}
        dict_current_species=compare_markers_with_spectra(set_of_current_markers, list_of_spectra, resolution)
        for m in set_of_current_markers:
            dict_all_species[m]=dict_current_species[m].copy()
            if len(dict_current_species[m])>=len(list_of_spectra)*min_ratio_spectra:
                dict_of_codes_ptm[(m.code,m.ptm)]+=1
        min_number_of_species=min_ratio_species*(len(set_of_taxa)-len(missing_species))
        set_of_filtered_markers={m for m in set_of_markers if dict_of_codes_ptm[(m.code,m.ptm)]>=min_number_of_species}
        for m in set_of_filtered_markers:
            if m in dict_all_species:
                m.comment=m.comment+" - "+ str(len(dict_all_species[m]))+"/"+str(dict_spectra[m.taxon_name])+" spectra :" + str(dict_all_species[m])
    return set_of_filtered_markers
