"""
    assignment.py
"""

import math
from scipy.stats import binom
import sys
from src import utils 
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


def p_success(spectrum, resolution):
    cover=0
    max_cover=0
    for peak in spectrum:
        if peak-resolution<max_cover: #overlap
            cover=+peak+resolution-max_cover
            max_cover=peak+resolution
        else: #no overlap
            cover+=2*resolution
            max_cover=peak+resolution
    return cover/(max(spectrum)-min(spectrum)+2*resolution)

def pvalue(k,n,p_success):
    return binom.sf(k-1,n, p_success)
    

"""
def compute_pvalue(spectrum, set_of_markers, assignation, score_to_taxid, resolution):
    taxid_to_pvalue={}
    score_max=max(score_to_taxid)
    for taxid in score_to_taxid[score_max]:
        set_of_masses={m for marker in set_of_markers if marker.taxid==taxid for m in marker.masses}
        mass_range=max(set_of_masses)-min(set_of_masses)
        mass_number=len(set_of_masses)
        taxid_to_pvalue[taxid]=pvalue(score_max, spectrum,  mass_range, mass_number, resolution)
    for taxid in assignation.taxid_to_optimal:
        if not assignation.taxid_to_optimal[taxid]:
            taxid_to_pvalue[taxid]=1
    return taxid_to_pvalue
 """

def assign_peaks_of_the_spectrum(spectrum, mass_taxid_name_list, resolution):  
    spectrum.sort()  # peaks are sorted according to their mass
    peak_to_taxid_name_ptm={} #key: peak, value: set of triplets (taxid,name of the marker,ptm)  
    current_j = 0
    for peak in spectrum:
        for j in range(current_j,len(mass_taxid_name_list)):
            diff_peak_pep = peak-mass_taxid_name_list[j][0]
            if utils.matching_masses(mass_taxid_name_list[j][0], peak, resolution):  # there is a match 
                utils.update_dictoset(peak_to_taxid_name_ptm, peak, mass_taxid_name_list[j][1])
            elif diff_peak_pep<0: # peptide mass is greater than peak mass 
                current_j = j 
                break # next peak 
    return peak_to_taxid_name_ptm

def number_of_different_peaks(peaks):
    #number of peaks with different name+ptm
    return len({(name,ptm) for (mass,name,ptm) in peaks})

def is_included(v,w):
   peaks_of_v={p[0] for  p in v}
   peaks_of_w={p[0] for p in w}
   return peaks_of_v<peaks_of_w

def assign_spectrum(spectrum, mass_taxid_name_list, set_of_markers, resolution, taxonomy, threshold, allsolutions):
    assignment=Assignment()
    assignment.spectrum_name=spectrum.name
    assignment.peak_to_taxid_name_ptm=assign_peaks_of_the_spectrum(spectrum, mass_taxid_name_list, resolution)
   
    if len(assignment.peak_to_taxid_name_ptm)<4:
        return assignment 
    
    taxid_to_peak_name_ptm={}
    for peak in assignment.peak_to_taxid_name_ptm:
        for (taxid,name,ptm) in  assignment.peak_to_taxid_name_ptm[peak]:
            utils.update_dictoset(taxid_to_peak_name_ptm, taxid, {(peak, name,ptm)})
            
    max_number_of_markers=max({len(v) for v in taxid_to_peak_name_ptm.values()})

    if max_number_of_markers<4:
        return assignment 

    # Remove included taxid and equivalent taxid
    excluded_taxid=set()
    equivalent_taxid={}
    for taxid in taxid_to_peak_name_ptm:
        for taxid2 in taxid_to_peak_name_ptm:
            if taxid_to_peak_name_ptm[taxid]<taxid_to_peak_name_ptm[taxid2]:
                excluded_taxid.add(taxid)
                continue
            if  taxid_to_peak_name_ptm[taxid]==taxid_to_peak_name_ptm[taxid2]:
                if taxid<taxid2:
                    utils.update_dictoset(equivalent_taxid, taxid, {taxid2})
                else:
                    continue
    for taxid in equivalent_taxid:
        excluded_taxid.update(equivalent_taxid[taxid])            
    for taxid in excluded_taxid:
        del taxid_to_peak_name_ptm[taxid]

    taxid_to_score={taxid:number_of_different_peaks(taxid_to_peak_name_ptm[taxid]) for taxid in taxid_to_peak_name_ptm} 
        
    # Find assignment with best P-value
    taxid_to_pvalue={}
    p=p_success(spectrum, resolution)
    for taxid in taxid_to_peak_name_ptm:
        number_of_matching_masses=taxid_to_score[taxid]
        total_number_of_masses=len({m for m in set_of_markers if m.taxid==taxid})
        taxid_to_pvalue[taxid]=pvalue(number_of_matching_masses, total_number_of_masses, p)

    best_pvalue=1
    for taxid in taxid_to_pvalue:
        if taxid_to_pvalue[taxid]<best_pvalue:
            best_pvalue=taxid_to_pvalue[taxid]
            best_score=taxid_to_score[taxid]
    set_of_optimal_taxid={taxid for taxid in taxid_to_pvalue if taxid_to_pvalue[taxid]==best_pvalue or taxid_to_score[taxid]>=best_score*threshold/100}
    
    results=[] # each element is a 4-uplet (lca, pvalue, score, set of equivalent taxid)
    combinations_of_peaks=[] # each element is a set of 3-uplets (peak, name,ptm) for the taxid of results  
    for taxid in set_of_optimal_taxid:
        if taxid in equivalent_taxid:
            eq_taxid=equivalent_taxid[taxid].union({taxid})
        else:
            eq_taxid={taxid}
        lca=taxonomy.hca(eq_taxid)
        results.append((lca, taxid_to_pvalue[taxid], taxid_to_score[taxid],eq_taxid))
        combinations_of_peaks.append(taxid_to_peak_name_ptm[taxid])     
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
        s=s+"\tPvalue \t #peaks  \t Assignment \t Rank \t Species\n"
    else:
        s=s+"\tPvalue\t #peaks \t Species\n"
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
                    s=s+str(round(name_to_peak[(name,ptm)],3))
                s=s+"\t"
            if taxonomy:
                lca=a.results[j][0]
                if lca!=None:
                    s=s+"{:.2e}".format(a.results[j][1])+"\t"+str(a.results[j][2])+"\t"+str(lca) + " ["+ B.taxidname[lca] + "]\t" + B.rank[lca] +"\t"
                else:
                    s=s+"\t None\t\t" 
            else:
                s=s+"{:.2e}".format(a.results[j][1])+"\t"+str(a.results[j][2])+"\t"
            for taxid in a.results[j][3]: 
                s=s+str(taxid)+" ["+B.taxidname[taxid]+"] "
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
        name_ptm_to_mass={}
        mass_name_ptm_to_taxid={}
        taxid_to_score={}
        for j in range(len(a.results)):
            for taxid in a.results[j][3]:
                taxid_to_score[taxid]= "{:.2e}".format(a.results[j][1])
                for (mass,name,ptm) in a.combinations_of_peaks[j]:
                    utils.update_dictoset(mass_name_ptm_to_taxid, (mass,name,ptm),{taxid})
                    utils.update_dictoset(name_ptm_to_mass, (name,ptm), {mass})
        for (name,ptm) in list_of_marker_full_names:
            if (name,ptm) in name_ptm_to_mass :
                # the peptide marker is present in the spectrum
                for mass in name_ptm_to_mass[(name,ptm)]:
                    #creation of a new line in the file
                    found=False
                    s=a.spectrum_name+"\t"+image(name,ptm)+"\t"+ str(mass)+"\t"+str(int(mass_to_intensity[mass]))+"\t"
                    for taxid in list_of_useful_taxid:
                        if taxid in mass_name_ptm_to_taxid[(mass,name,ptm)] :
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
            if taxid in taxid_to_score:
                s=s+"\t"+taxid_to_score[taxid]
                utils.update_dictoset(score_to_taxid, taxid_to_score[taxid], {taxid})
            else:
                s=s+"\t"
        s=s+"\t"
        for score in score_to_taxid:
            s=s+score+" -> "
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
        current_assignment = assign_spectrum(spectrum, mass_taxid_name_list, set_of_markers, error, B, threshold,allsolutions)
        list_of_assignments.append(current_assignment)

    # selection of all useful taxid
    set_of_useful_taxid={z for a in list_of_assignments for y in a.results for z in y[3] }
    list_of_useful_taxid=list(set_of_useful_taxid)
    list_of_useful_taxid.sort()

    # selection of all useful markers
    set_of_marker_full_names=set() #{x[1] for x in mass_taxid_name_list}
    for (x,y) in  mass_taxid_name_list:
        set_of_marker_full_names.update(y)
    if len(set_of_marker_full_names)>40:
        # only markers with matches are displayed
        set_of_marker_full_names=set()
        for a in list_of_assignments:
            for y in a.combinations_of_peaks:
                set_of_marker_full_names.update(y)
    list_of_marker_full_names=list({(str(b),str(c)) for (a,b,c) in set_of_marker_full_names})
    list_of_marker_full_names.sort()
   
    create_main_result_file(output, list_of_assignments, taxonomy, B,list_of_marker_full_names)
    create_detail_result_file(output2, list_of_useful_taxid,list_of_spectra, list_of_assignments, taxonomy, B,list_of_marker_full_names)

 
    
