"""
    assignment.py
"""
import json
import math
from scipy.stats import binom
import sys
from src import utils 
from src import markers

class Annotated_peak(object):
    def __init__(self, mass=0, intensity=0, marker=None):
        self.mass = mass
        self.intensity = intensity  # list of (peak,name,ptm)
        self.marker = marker 
        
    def __eq__(self, other):
        return self.mass == other.mass and self.marker.code == other.marker.code

    def __hash__(self):
        return hash((self.mass, self.intensity, self.marker.code, self.marker.ptm))

class Taxon(object):
    def __init__(self, id=None, name=None):
        self.id = id
        self.name = name
        
    def __eq__(self, other):
        return self.id == other.id 

    def __hash__(self):
        return hash((self.id, self.name))
        
class Assignment(object):
    def __init__(self, spectrum_name="", spectrum_length=0, peaks=[], taxa=[], lca=None, lca_name=None, lca_rank=None, score=None, hca=None, hca_name=None, hca_rank=None, pvalue=None):
        self.spectrum_name = spectrum_name
        self.spectrum_length = spectrum_length
        self.peaks=peaks # list of Annotated_peaks
        self.taxa=taxa # list of objects Taxon, all Taxon have the same markers 
        self.lca = lca
        self.lca_name = lca_name
        self.lca_rank = lca_rank
        self.score = score # number of peaks in the assignment
        self.hca = hca
        self.hca_name = hca_name
        self.hca_rank = hca_rank
        self.pvalue = pvalue

        
def mcc(TP,TN, FP, FN):
    return (TN*TP-FP*FN)/math.sqrt((TN+FN)*(FP+TP)*(TN+FP)*(FN+TP))


def p_success(spectrum, resolution):
    cover=0
    max_cover=0
    if  resolution<1.1:
        delta=resolution # daltons
    else:
        delta=resolution/400  #ppm. Average m/z is supposed to be 2500
    for peak in spectrum:
        if peak.mass-resolution<max_cover: #overlap
            cover=+peak.mass+resolution-max_cover
            max_cover=peak.mass+resolution
        else: #no overlap
            cover+=2*resolution
            max_cover=peak.mass+resolution
    return cover/(max({peak.mass for peak in spectrum})-min({peak.mass for peak in spectrum})+2*resolution)

def pvalue_f(k,n,p_success):
    return binom.sf(k-1,n, p_success)
    
def assign_peaks_of_the_spectrum(spectrum, mass_markers_list, resolution):  
    spectrum.sort()  # peaks are sorted according to their mass
    peak_to_markers={} #key: peak, value: set of markers   
    current_j = 0
    for peak in spectrum:
        for j in range(current_j,len(mass_markers_list)):
            diff_peak_pep = peak.mass-mass_markers_list[j][0]
            if utils.matching_masses(mass_markers_list[j][0], peak.mass, resolution):  # there is a match 
                utils.update_dictoset(peak_to_markers, peak, mass_markers_list[j][1])
            elif diff_peak_pep<0: # peptide mass is greater than peak mass 
                current_j = j 
                break # next peak
    return peak_to_markers

def number_of_different_peaks(peaks):
    #number of peaks with different name+ptm
    return len({(p.marker.code,p.marker.ptm) for p in peaks})

def is_included(v,w):
   peaks_of_v={p.mass for  p in v}
   peaks_of_w={p.mass for p in w}
   return peaks_of_v<peaks_of_w

def assign_spectrum(spectrum, mass_markers_list, set_of_markers, resolution, taxonomy, threshold, allsolutions):
   # mass_taxid_name_list: contains the list of markers sorted by mass
    peak_to_markers=assign_peaks_of_the_spectrum(spectrum, mass_markers_list, resolution)
    if len(peak_to_markers)<4:
        return []

    taxid_to_annotated_peaks={}
    for peak in peak_to_markers:
        for m in peak_to_markers[peak]:
            utils.update_dictoset(taxid_to_annotated_peaks, m.taxid, {Annotated_peak(peak.mass, peak.intensity, m)})
            
    max_number_of_markers=max({len(v) for v in taxid_to_annotated_peaks.values()}) # incorrect, à revoir

    if max_number_of_markers<4:
        return [] 

    # Remove included taxid and equivalent taxid
    excluded_taxid=set()
    equivalent_taxid={}
    for taxid in taxid_to_annotated_peaks:
        for taxid2 in taxid_to_annotated_peaks:
            if  taxid_to_annotated_peaks[taxid]==taxid_to_annotated_peaks[taxid2]:
                if taxid < taxid2:
                    utils.update_dictoset(equivalent_taxid, taxid, {taxid2})
                else:
                    continue
            if taxid_to_annotated_peaks[taxid]<taxid_to_annotated_peaks[taxid2]:
                excluded_taxid.add(taxid)
                continue
    excluded_taxid.update({t for taxid in equivalent_taxid for t in equivalent_taxid[taxid]})
    for taxid in excluded_taxid:
        del taxid_to_annotated_peaks[taxid]

    taxid_to_score={taxid:number_of_different_peaks(taxid_to_annotated_peaks[taxid]) for taxid in taxid_to_annotated_peaks}
        
    # Find assignment with best P-value
    taxid_to_pvalue={}
    p=p_success(spectrum, resolution)
    for taxid in taxid_to_annotated_peaks:
        number_of_matching_masses=taxid_to_score[taxid]
        total_number_of_masses=len({m for m in set_of_markers if m.taxid==taxid})
        taxid_to_pvalue[taxid]=pvalue_f(number_of_matching_masses, total_number_of_masses, p)
    best_pvalue=1
    for taxid in taxid_to_pvalue:
        if taxid_to_pvalue[taxid]<best_pvalue:
            best_pvalue=taxid_to_pvalue[taxid]
            best_score=taxid_to_score[taxid]
    set_of_optimal_taxid={taxid for taxid in taxid_to_pvalue if taxid_to_pvalue[taxid]==best_pvalue or taxid_to_score[taxid]>=best_score*threshold/100}

    list_of_assignments=[]
    for taxid in set_of_optimal_taxid:
        if taxid in equivalent_taxid:
            eq_taxid=equivalent_taxid[taxid].union({taxid})
        else:
            eq_taxid={taxid}
        if taxonomy:
                lca=taxonomy.lca(eq_taxid)
        pvalue= taxid_to_pvalue[taxid]
        score=taxid_to_score[taxid]
        taxids=list({Taxon(taxid, None) for taxid in eq_taxid})
        set_of_peaks=taxid_to_annotated_peaks[taxid]     
        a=Assignment(spectrum.name, len(spectrum), list(set_of_peaks), taxids, lca, None, None, score, None, None, None, pvalue)
        list_of_assignments.append(a)
    return list_of_assignments 

    
def image(name, ptm):
    return name+ "-"+ptm



def create_json_result_file(output, list_of_assignments):
    json_file=open(output+".json", "w")
    list_of_serialized_objects=[]
    for a in list_of_assignments:
        a_dict = vars(a)
        serialized_list_of_taxa=[vars(t) for t in a_dict["taxa"]]
        serialized_list_of_peaks=[{"mass":p.mass, "intensity":p.intensity, "code":p.marker.code, "PTM":p.marker.ptm, "sequence":p.marker.sequence, "protein":p.marker.protein, "begin":p.marker.begin, "end":p.marker.end} for p in a_dict["peaks"]]
        a_dict["peaks"]=serialized_list_of_peaks
        a_dict["taxa"]= serialized_list_of_taxa=[vars(t) for t in a_dict["taxa"]]
        list_of_serialized_objects.append(a_dict)
    json.dump(list_of_serialized_objects, json_file, indent=4)
    json_file.close()
        
def create_main_result_file(output, list_of_assignments, taxonomy):
     
    set_of_marker_full_names={(m.marker.code, m.marker.ptm) for a in list_of_assignments for m in a.peaks}
    list_of_marker_full_names=list(set_of_marker_full_names)
    
    # Heading
    f1 = open(output+".tsv", "w")
    s=""
    for (name,ptm) in  list_of_marker_full_names:
        s=s+"\t"+image(name,ptm)
    if taxonomy:
        s=s+"\tPvalue \t #peaks  \t Assignment \t Rank \t Uncertainty \t Rank \t Species\n"
    else:
        s=s+"\tPvalue\t #peaks \t Species\n"
    f1.write(s)
    
    # Body
    for a in list_of_assignments:
        if len(a.taxa)==0: # no assignment A REVOIR A REVOIR
            s=a.spectrum_name+"\t"*len( list_of_marker_full_names)+"\t\t 0\t None \n"
            f1.write(s)
            continue
        name_to_peak={(p.marker.code,p.marker.ptm):p.mass for p in a.peaks}
        s=a.spectrum_name+"\t"
        for (name,ptm) in list_of_marker_full_names:
            if (name,ptm) in  name_to_peak:
                s=s+str(round(name_to_peak[(name,ptm)],3))
            s=s+"\t"
        if taxonomy:
            if a.lca!=None:
                s=s+"{:.2e}".format(a.pvalue)+"\t"+str(a.score)+"\t"+str(a.lca) + " ["+ a.lca_name+ "]\t" + a.lca_rank +"\t"+str(a.hca)+ " ["+ a.hca_name + "]\t" + a.hca_rank +"\t"
            else:
                s=s+"\t None\t\t\t\t" 
        else:
            s=s+"{:.2e}".format(a.pvalue)+"\t"+str(a.score)+"\t"
        for t in a.taxa: 
            s=s+str(t.id)+" ["+t.name+"] "
        s=s+"\n"
        f1.write(s)
    f1.close()

def create_detail_result_file(output_detail, list_of_assignments, taxonomy, B):

    # selection of all useful taxid
    set_of_useful_taxid={z.id for a in list_of_assignments for z in a.taxa}
    list_of_useful_taxid=list(set_of_useful_taxid)
    list_of_useful_taxid.sort()

    # assignments are clustered per spectrum
    dict_of_assignments={a.spectrum_name:{} for a in list_of_assignments}
    dict_of_pvalues={a.spectrum_name:{} for a in list_of_assignments}
    dict_of_equivalent_species={a.spectrum_name:[] for a in list_of_assignments}
    for a in list_of_assignments:
        dict_of_equivalent_species[a.spectrum_name].append((a.pvalue, a.taxa))
        for peak in a.peaks:
            utils.update_dictoset(dict_of_assignments[a.spectrum_name], peak, {z.id for z in a.taxa})
        for taxon in a.taxa:
            dict_of_pvalues[a.spectrum_name][taxon.id]=a.pvalue
            
    f2 =open(output_detail, "w")
    # heading
    s="spectrum\t marker\t m/z \t intensity  "
    for taxid in  list_of_useful_taxid:
        s=s+"\t"+str(taxid)
    s=s+"\t Best Guess\n"
    f2.write(s)
    # body
    for sp in dict_of_assignments:
        for peak in dict_of_assignments[sp]:
            s=sp+"\t"+image(peak.marker.code,peak.marker.ptm)+"\t"+ str(peak.mass)+"\t"+str(peak.intensity)+"\t"
            for taxid in list_of_useful_taxid:
                if taxid in dict_of_assignments[sp][peak] :
                    s=s+"X\t"
                else:
                    s=s+"\t"
            s=s+"\n"
            f2.write(s)
        s=sp+"\t SCORE \t\t"
        for taxid in list_of_useful_taxid:
            if taxid in dict_of_pvalues[sp]:
                s=s+"\t"+"{:.2e}".format(dict_of_pvalues[sp][taxid])
            else:
                s=s+"\t"
        s=s+"\t"
        for (pvalue, taxa) in dict_of_equivalent_species[sp]:
            s=s+"{:.2e}".format(pvalue)+": "
            for taxon in taxa:
                s=s+str(taxon.id)+"("+str(taxon.name)+") "
        f2.write(s+"\n")

    f2.close()

def assign_all_spectra(list_of_spectra, set_of_markers, error, taxonomy, B, threshold, allsolutions, output, output_detail):

    # elements of the list are 2-uplets of the form  (mass, {(taxid, code, PTM)})
    mass_markers_list=markers.sort_by_masses(set_of_markers)
    list_of_assignments=[]
    for spectrum in list_of_spectra:
        list_of_assignments.extend(assign_spectrum(spectrum, mass_markers_list, set_of_markers, error, B, threshold, allsolutions))

    # completion of assignments
    for a in list_of_assignments:
        list_of_taxa=[]
        for t in a.taxa:
            list_of_taxa.append(Taxon(t.id,B.name[t.id]))
        a.taxa = list_of_taxa                     
        if taxonomy and a.lca is not None:
            a.lca_rank = B.rank[a.lca]
            a.lca_name = B.name[a.lca]
            a.hca =  B.unary_ancestor(a.lca)
            if a.hca is not None:
                a.hca_rank = B.rank[a.hca]
                a.hca_name = B.name[a.hca]
   
    create_main_result_file(output, list_of_assignments, taxonomy)
    create_detail_result_file(output_detail, list_of_assignments, taxonomy,B)
    create_json_result_file(output, list_of_assignments)
