"""
    assignment.py
"""

import math
import sys
from src import utils as ut
from src import markers

class Assignment(object):
    def __init__(self, spectrum_name="",  peak_to_taxid_name_ptm={}, combinations_of_peaks=[], results=set(), taxid_to_peak_name_ptm={}):
        self.spectrum_name = spectrum_name
        self.peak_to_taxid_name_ptm=peak_to_taxid_name_ptm
        self.combinations_of_peaks=combinations_of_peaks #list of combinations of peaks. Peak=(peak,name,ptm)
        self.results=results # list of 3-uplets (lca, score, set of 2-uplets (pvalue, taxid)), all taxid  have the same markers 
        self.taxid_to_peak_name_ptm=taxid_to_peak_name_ptm # index in combinations_of_peaks
    
def mcc(TP,TN, FP, FN):
    return (TN*TP-FP*FN)/math.sqrt((TN+FN)*(FP+TP)*(TN+FP)*(FN+TP))

def matching_masses(theoretical_peak, experimental_peak, resolution):
    delta= abs(theoretical_peak - experimental_peak)
    if resolution<1:
        return delta<=resolution
    else:
        return delta/theoretical_peak<=resolution/1000000
    
def assign_peaks_of_the_spectrum(spectrum, mass_taxid_name_list, resolution):  
    spectrum.sort()  # peaks are sorted according to their mass
    peak_to_taxid_name_ptm={} #key: peak, value: set of triplets (taxid,name of the marker,ptm)  
    current_j = 0
    for peak in spectrum:
        for j in range(current_j,len(mass_taxid_name_list)):
            diff_peak_pep = peak-mass_taxid_name_list[j][0]
            if matching_masses(mass_taxid_name_list[j][0], peak, resolution):  # there is a match 
                ut.update_dictoset(peak_to_taxid_name_ptm, peak, mass_taxid_name_list[j][1])
            elif diff_peak_pep<0: # peptide mass is greater than peak mass 
                current_j = j 
                break # next peak 

    return peak_to_taxid_name_ptm


def count_number_of_masses_per_taxid(list_of_taxid, mass_taxid_name_list):
    taxid_to_number_of_masses_dict={}
    for taxid in list_of_taxid:
        taxid_to_number_of_masses_dict[taxid]=0
    for (mass, set_of_pairs) in mass_taxid_name_list:
        set_of_taxid={x[0] for x in set_of_pairs}
        for taxid in list_of_taxid:
            if taxid in set_of_taxid:
                taxid_to_number_of_masses_dict[taxid]+=1
    return  taxid_to_number_of_masses_dict


def is_included(v,w):
   peaks_of_v={p[0] for  p in v}
   peaks_of_w={p[0] for p in w}
   return peaks_of_v<peaks_of_w

def assign_spectrum(spectrum, mass_taxid_name_list, resolution, taxonomy, threshold, allsolutions):
    assignment=Assignment()
    assignment.spectrum_name=spectrum.name
    assignment.peak_to_taxid_name_ptm=assign_peaks_of_the_spectrum(spectrum, mass_taxid_name_list, resolution)
   
    if len(assignment.peak_to_taxid_name_ptm)<2:
        return assignment 
    
    taxid_to_peak_name_ptm={}
    for peak in assignment.peak_to_taxid_name_ptm:
        for (taxid,name,ptm) in  assignment.peak_to_taxid_name_ptm[peak]:
            ut.update_dictoset(taxid_to_peak_name_ptm, taxid, {(peak, name,ptm)})

    max_number_of_markers=max({len(taxid_to_peak_name_ptm[taxid]) for taxid in taxid_to_peak_name_ptm})

    if max_number_of_markers<2:
        return assignment
 
    set_of_combinations=set()
    for v in taxid_to_peak_name_ptm.values():
        if len(v)>=max_number_of_markers*threshold/100:
            set_of_combinations.add(frozenset(v))
    if allsolutions:
        combinations_of_peaks=list(set_of_combinations)
    else:
        combinations_of_peaks=[]
        for v in set_of_combinations:
            is_optimal=True
            for w in  set_of_combinations: 
                if is_included(v,w):
                    is_optimal=False
            if is_optimal:
                combinations_of_peaks.append(v)
    
    for taxid in taxid_to_peak_name_ptm:
        try:
            i=combinations_of_peaks.index(taxid_to_peak_name_ptm[taxid])
            taxid_to_peak_name_ptm[taxid]=i
        except ValueError:
            taxid_to_peak_name_ptm[taxid]=-1
            
    dict_of_equivalent_taxid={}
    for taxid in taxid_to_peak_name_ptm:
        if taxid_to_peak_name_ptm[taxid]>=0:
             ut.update_dictoset(dict_of_equivalent_taxid, taxid_to_peak_name_ptm[taxid], {taxid})
    
    results=[]
    for i in range(len(combinations_of_peaks)):
        score=len(combinations_of_peaks[i])
        lca=taxonomy.hca(dict_of_equivalent_taxid[i])
        set_of_taxid=set()
        for taxid in  dict_of_equivalent_taxid[i]:
            set_of_taxid.add((0.1, taxid))
        results.append((lca, score, set_of_taxid))
            
    assignment.results=results
    assignment.taxid_to_peak_name_ptm=taxid_to_peak_name_ptm
    assignment.combinations_of_peaks=combinations_of_peaks
    return assignment 

def image(name, ptm):
    return name+ "-"+ptm


def create_main_result_file(output, list_of_assignments, taxonomy, B,list_of_marker_full_names):
    # Heading
    f1 = open(output, "w")
    #f1.write("Spectra:\t"+spectra_path+"\nMarkers:\t"+ marker_source + "\nResolution:\t"+str(resolution)+"\n" )
    s=""
    for (name,ptm) in  list_of_marker_full_names:
        s=s+"\t"+image(name,ptm)
    if taxonomy:
        s=s+"\tScore \t Assignment \t Rank \t Species\n"
    else:
        s=s+"\tScore \t Species\n"
    f1.write(s)

    # Body
    for i in range(len(list_of_assignments)):
        a=list_of_assignments[i]
        if len(a.results)==0: # no assignment
            s=a.spectrum_name+"\t"*len( list_of_marker_full_names)+"\t0\t None \n"
            f1.write(s)
            continue
        for j in range(len(a.results)): 
            name_to_peak={(name,ptm):peak for (peak,name,ptm) in a.combinations_of_peaks[j]}
            s=a.spectrum_name+"\t"
            for (name,ptm) in list_of_marker_full_names:
                if (name,ptm) in  name_to_peak:
                    s=s+str(name_to_peak[(name,ptm)])
                s=s+"\t"
            if taxonomy:
                lca=a.results[j][0]
                if lca!=None:
                    s=s+str(a.results[j][1])+"\t"+str(lca) + " ["+ B.taxidname[lca] + "]\t" + B.rank[lca] +"\t"
                else:
                    s=s+"\t None\t\t" 
            else:
                s=s+str(a.results[j][1])+"\t"
            for taxid in a.results[j][2]: 
                s=s+str(taxid)+" ["+B.taxidname[taxid[1]]+"] "
            s=s+"\n"
            f1.write(s)
    f1.close()

def create_detail_result_file(output2, list_of_useful_taxid,list_of_spectra, list_of_assignments, taxonomy, B,list_of_marker_full_names):
#
    f2 =open(output2, "w")
    
    #f2.write("Spectra:\t"+spectra_path+"\nMarkers:\t"+ marker_source + "\nResolution:\t"+str(resolution)+"\n" )
    #heading
    s="spectrum\t marker\t m/z \t intensity  "
    for taxid in  list_of_useful_taxid:
        s=s+"\t"+str(taxid)
    s=s+"\t Best Guess\n"
    f2.write(s)
    # Secondary file - body
    for i in range(len(list_of_spectra)):
        a=list_of_assignments[i]
        sp=list_of_spectra[i]
        mass_to_intensity={p.mass:p.intensity for p in sp.peaks}
        name_to_peak={}
        for peak in  a.peak_to_taxid_name_ptm:
            for (taxid, name, ptm) in  a.peak_to_taxid_name_ptm[peak]:
                ut.update_dictoset(name_to_peak, (name,ptm), {peak})
        for (name,ptm) in list_of_marker_full_names:
            if (name,ptm) in name_to_peak :
                # the peptide marker is present in the spectrum
                for peak in name_to_peak[(name,ptm)]:
                    found=False
                    #creation of a new line in the file
                    s=a.spectrum_name+"\t"+image(name,ptm)+"\t"+ str(peak)+"\t"+str(int(mass_to_intensity[peak]))+"\t"
                    for taxid in list_of_useful_taxid:
                        if taxid in a.taxid_to_peak_name_ptm and a.taxid_to_peak_name_ptm[taxid]>=0 and (taxid,name,ptm) in a.peak_to_taxid_name_ptm[peak] :
                            s=s+"X\t"
                            found=True
                        else:
                            s=s+"\t"
                    s=s+"\n"
                    if found :
                        f2.write(s)
        s=a.spectrum_name+"\t SCORE \t\t"
        score_to_taxid={}
        for taxid in list_of_useful_taxid:
            if taxid in a.taxid_to_peak_name_ptm and a.taxid_to_peak_name_ptm[taxid]>=0: # A revoir
                indice=a.taxid_to_peak_name_ptm[taxid]
                score=(a.results[indice])[1]
                s=s+"\t"+str(int(score))
                ut.update_dictoset(score_to_taxid, score, {taxid})
            else:
                s=s+"\t"
        s=s+"\t"
        for score in score_to_taxid:
            s=s+str(score)+" -> "
            for taxid in score_to_taxid[score]:
                s=s+taxid+" ["+B.taxidname[taxid]+"] "
        s=s+"\n"
        f2.write(s)

    f2.close()

def assign_all_spectra(list_of_spectra, set_of_markers, error, taxonomy, B, threshold, allsolutions, output,output2):

    # elements of the list are 2-uplets of the form  (mass, {(taxid, code, PTM)})
    mass_taxid_name_list=markers.sort_by_masses(set_of_markers)
    
    list_of_assignments=[]
    for spectrum in list_of_spectra:
        current_assignment = assign_spectrum(spectrum, mass_taxid_name_list, error, B, threshold,allsolutions)
        list_of_assignments.append(current_assignment)

    # selection of all useful taxid
    set_of_useful_taxid={z[1] for a in list_of_assignments for y in a.results for z in y[2] }
    list_of_useful_taxid=list(set_of_useful_taxid)
    list_of_useful_taxid.sort()

    # selection of all useful markers
    set_of_marker_full_names=set() #{x[1] for x in mass_taxid_name_list}
    for (x,y)  in  mass_taxid_name_list:
        set_of_marker_full_names.update(y)
    if len(set_of_marker_full_names)>40:
        # only markers with matches are displayed
        set_of_marker_full_names=set()
        for a in list_of_assignments:
            for y in a.combinations_of_peaks:
                set_of_marker_full_names.update(y)
    list_of_marker_full_names=list({(b,c) for (a,b,c) in set_of_marker_full_names})
    list_of_marker_full_names.sort()
   
    create_main_result_file(output, list_of_assignments, taxonomy, B,list_of_marker_full_names)
    create_detail_result_file(output2, list_of_useful_taxid,list_of_spectra, list_of_assignments, taxonomy, B,list_of_marker_full_names)

 
    
