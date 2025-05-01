import copy
import os
from os import listdir
from os.path import join

from src import utils 
from src import mass_spectrum as ms
from src import compute_masses
from src import sequences
from src import markers


# input:
# list of markers, represented by a sorted list of masses (in increasing order)
# one spectrum
# output: list of matching peaks
def compare_set_of_markers_for_one_spectrum(mass_list, spectrum, resolution):
    spectrum.sort()  
    peak_list=[None for m in mass_list] 
    current_j = 0
    for peak in spectrum.peaks:
        for j in range(current_j,len(mass_list)):
            difference = float(peak.mass) - float(mass_list[j][0])
            if utils.matching_masses(peak.mass, mass_list[j][0], resolution): # on a un match
                peak_list[j]=peak.intensity 
            elif difference<0: # la masse du peptide est supérieure à la masse du pic  
                current_j = j # on reprend au même endroit où le pic d'avant s'est arrêté
                break # on ne va pas plus loin dans la liste des masses des peptides et on passe au pic suivant
    return peak_list


def compare_markers_with_spectra(set_of_markers, list_of_spectra, resolution):
    mass_list=[(m.mass(),m) for m in set_of_markers]
    mass_list.sort(key=lambda x:x[0])
    dict_markers={}
    dict_intensities={m[1]:[None for s in list_of_spectra] for m in mass_list} 
    for (i,spectrum) in enumerate(list_of_spectra):
        peak_list= compare_set_of_markers_for_one_spectrum(mass_list, spectrum, resolution)
        for j in range(0,len(mass_list)):
            if peak_list[j] is not None:
                utils.update_dictoset(dict_markers, mass_list[j][1], {spectrum.name})
                dict_intensities[mass_list[j][1]][i]=peak_list[j]
        
    return dict_markers, dict_intensities

    
# les marqueurs sont supposés être de la même espèce que les spectres
def filter_set_of_markers(set_of_markers, list_of_spectra, resolution, min_nb_spectra=0):
    dict_markers, dict_intensities = compare_markers_with_spectra(set_of_markers, list_of_spectra, resolution)
    set_of_confirmed_markers=set()
    for m in dict_markers:
        found_spectra=dict_markers[m]
        if len(found_spectra)>=min_nb_spectra:
            m=markers.update_comment(m, str(len(found_spectra))+"/"+str(len(list_of_spectra))+" spectra :" + str(found_spectra)+ " + ")
            m.field["Status"]="MS"
            set_of_confirmed_markers.add(m)
    
    return set_of_confirmed_markers



def filter_set_of_markers_multispecies(set_of_markers, spectra_dir, resolution, min_ratio_spectra=0, min_ratio_species=0):
    dict_of_codes_ptm={(m.code(),m.PTM()):0 for m in set_of_markers}
    set_of_taxa={m.taxon_name() for m in set_of_markers}
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
        set_of_current_markers={m for m in set_of_markers if m.taxon_name()==taxon}
        dict_current_species=compare_markers_with_spectra(set_of_current_markers, list_of_spectra, resolution)
        for m in set_of_current_markers:
            dict_all_species[m]=dict_current_species[m].copy()
            if len(dict_current_species[m])>=len(list_of_spectra)*min_ratio_spectra:
                dict_of_codes_ptm[(m.code(),m.PTM())]+=1
        min_number_of_species=min_ratio_species*(len(set_of_taxa)-len(missing_species))
        set_of_filtered_markers={m for m in set_of_markers if dict_of_codes_ptm[(m.code(),m.PTM())]>=min_number_of_species}
        for m in set_of_filtered_markers:
            if m in dict_all_species:
                m=markers.update_comment(m, " - "+ str(len(dict_all_species[m]))+"/"+str(dict_spectra[m.taxon_name])+" spectra :" + str(dict_all_species[m]))
    return set_of_filtered_markers



