"""
    assignment.py
"""
import json
import math
from scipy.stats import binom
import csv
import sys
from src import utils 
from src import markers
from src import config

class Annotated_peak(object):
    def __init__(self, mass=0, intensity=0, marker=None):
        self.mass = mass
        self.intensity = intensity  # list of (peak,name,ptm)
        self.marker = marker 
        
    def __eq__(self, other):
        return self.mass == other.mass and self.marker.code() == other.marker.code()

    def __hash__(self):
        return hash((self.mass, self.intensity, self.marker.code(), self.marker.PTM()))

class Taxon(object):
    def __init__(self, id=None, name=None):
        self.id = id
        self.name = name
        
    def __eq__(self, other):
        return self.id == other.id 

    def __hash__(self):
        return hash((self.id, self.name))
        
class Assignment(object):
    def __init__(self, spectrum_name="", spectrum_length=0, peaks=[], taxa=[], lca=None, lca_name=None, lca_rank=None, score=None, hca=None, hca_name=None, hca_rank=None, pvalue=1):
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
    return len({(p.marker.code(),p.marker.PTM()) for p in peaks})

def is_included(v,w):
    # not used ?
   peaks_of_v={p.mass for  p in v}
   peaks_of_w={p.mass for p in w}
   return peaks_of_v<peaks_of_w


def is_better(v,w):
    masses_of_v={(a.mass):a.intensity for a in v}
    masses_of_w={(a.mass):a.intensity for a in w}
    
    IV={masses_of_v[k] for k  in masses_of_v.keys() - masses_of_w.keys()}
    IW ={masses_of_w[k] for k in masses_of_w.keys() - masses_of_v.keys()}
    if len(IV) == 0 or len(IW)==0:
        return False
    else:
        return  min(IV)>max(IW)
    
def no_assignment(spectrum):
    return [Assignment(spectrum.name)]
    
def assign_spectrum(spectrum, mass_markers_list, set_of_markers, resolution, taxonomy, threshold, allsolutions, minimum_number_of_peaks):
   # mass_taxid_name_list: contains the list of markers sorted by mass
    peak_to_markers=assign_peaks_of_the_spectrum(spectrum, mass_markers_list, resolution)
    if len(peak_to_markers)<minimum_number_of_peaks:
        return no_assignment(spectrum)

    taxid_to_annotated_peaks={}
    for peak in peak_to_markers:
        for m in peak_to_markers[peak]:
            utils.update_dictoset(taxid_to_annotated_peaks, m.taxid(), {Annotated_peak(peak.mass, peak.intensity, m)})
            
    max_number_of_markers=max({len(v) for v in taxid_to_annotated_peaks.values()}) # incorrect, à revoir

    if max_number_of_markers < minimum_number_of_peaks:
        return no_assignment(spectrum)

    # Remove included taxid and equivalent taxid
    excluded_taxid=set()
    equivalent_taxid={}
    for taxid in taxid_to_annotated_peaks:
        if len(taxid_to_annotated_peaks[taxid])<minimum_number_of_peaks:
            excluded_taxid.add(taxid)
            continue
        for taxid2 in taxid_to_annotated_peaks:
            if taxid < taxid2 and taxid_to_annotated_peaks[taxid]==taxid_to_annotated_peaks[taxid2]:
                    utils.update_dictoset(equivalent_taxid, taxid, {taxid2})
            if not allsolutions and taxid_to_annotated_peaks[taxid]<taxid_to_annotated_peaks[taxid2]:
                excluded_taxid.add(taxid)
           
    excluded_taxid.update({t for taxid in equivalent_taxid for t in equivalent_taxid[taxid]})
    for taxid in excluded_taxid:
        del taxid_to_annotated_peaks[taxid]

    taxid_to_score={taxid:number_of_different_peaks(taxid_to_annotated_peaks[taxid]) for taxid in taxid_to_annotated_peaks}
    for taxid in taxid_to_score:
        if taxid_to_score[taxid]<minimum_number_of_peaks:
            del taxid_to_annotated_peaks[taxid]
        
    # Find assignment with best P-value
    taxid_to_pvalue={}
    p=p_success(spectrum, resolution)
    for taxid in taxid_to_annotated_peaks:
        number_of_matching_masses=taxid_to_score[taxid]
        if number_of_matching_masses==0:
            taxid_to_pvalue[taxid]=1
        elif number_of_matching_masses<5:
            taxid_to_pvalue[taxid]=1/number_of_matching_masses
        else:
            total_number_of_masses=len({m.mass() for m in set_of_markers if m.taxid()==taxid})# à initialiser avant pour plus d'efficacite
            taxid_to_pvalue[taxid]=pvalue_f(number_of_matching_masses, total_number_of_masses, p)
    best_pvalue=1
    for taxid in taxid_to_pvalue:
        if taxid_to_pvalue[taxid]<best_pvalue:
            best_pvalue=taxid_to_pvalue[taxid]
            best_score=taxid_to_score[taxid]
    set_of_optimal_taxid={taxid for taxid in taxid_to_pvalue if taxid_to_pvalue[taxid]==best_pvalue or taxid_to_score[taxid]>=best_score*threshold/100}
    
    if not allsolutions and threshold==100:
        excluded_taxid=set()
        for taxid in set_of_optimal_taxid:
            for taxid2 in set_of_optimal_taxid:
                if is_better(taxid_to_annotated_peaks[taxid2], taxid_to_annotated_peaks[taxid]):
                    excluded_taxid.add(taxid)
        set_of_optimal_taxid= set_of_optimal_taxid - excluded_taxid

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
    if ptm is None:
        return name
    else:
        return name+ "-"+ ptm



def create_json_result_file(output, list_of_assignments):
    json_file=open(output+".json", "w")
    list_of_serialized_objects=[]
    for a in list_of_assignments:
        a_dict = vars(a)
        serialized_list_of_taxa=[vars(t) for t in a_dict["taxa"]]
        serialized_list_of_peaks=[{"mass":p.mass, "intensity":p.intensity, "code":p.marker.code(), "PTM":p.marker.PTM(), "sequence":p.marker.sequence(), "protein":p.marker.protein(), "begin":p.marker.begin(), "end":p.marker.end()} for p in a_dict["peaks"]]
        a_dict["peaks"]=serialized_list_of_peaks
        a_dict["taxa"]= serialized_list_of_taxa=[vars(t) for t in a_dict["taxa"]]
        list_of_serialized_objects.append(a_dict)
    json.dump(list_of_serialized_objects, json_file, indent=4)
    json_file.close()


def create_spectral_file(output, list_of_assignments, list_of_spectra):
# list_of_assignments and list_of_spectra should share the same indices !! 
      with open("spectra"+output+".tsv", mode='w', newline='') as outfile:
            writer = csv.writer(outfile, delimiter='\t')
            writer.writerow(['spectrum', 'm/z', 'intensity', 'marker'])
            for i, spectrum in enumerate(list_of_spectra):
                mass_dict={peak.mass:peak.marker for peak in list_of_assignments[i].peaks} 
                for peak in spectrum:
                    row=[spectrum.name, str(peak.mass), peak.intensity]
                    if  peak.mass in mass_dict:
                        row.append(mass_dict[peak.mass].code)
                    writer.writerow(row)
               
def extract_sorted_markers(list_of_assignments):
    set_of_marker_codes={m.marker.code() for a in list_of_assignments for m in a.peaks}
    list_of_marker_full_names=list({(m.marker.code(), m.marker.PTM()) for a in list_of_assignments for m in a.peaks})
    list_of_marker_codes=config.sort_headers(set_of_marker_codes)
    list_of_marker_full_names.sort(key=lambda x:(list_of_marker_codes.index(x[0]),x[1]))
    return list_of_marker_full_names
        
def create_main_result_file(output, list_of_assignments, taxonomy):
    list_of_marker_full_names= extract_sorted_markers(list_of_assignments)
    # Heading
    f1 = open(output+".tsv", "w")
    if taxonomy:
        s="Spectrum \t Assignment \t Rank \t Maximal clade  \t Rank \t Species \t Pvalue \t #peaks "
    else:
        s="Spectrum \t Species \tPvalue\t #peaks "
    for (name,ptm) in  list_of_marker_full_names:
        s=s+"\t"+image(name,ptm)
    s=s+"\n"
    f1.write(s)
    
    # Body
    for a in list_of_assignments:
        if len(a.taxa)==0: # no assignment
            if taxonomy:
                s=a.spectrum_name+" \t No assignment \t \t \t \t \t \t"
            else:
                s=a.spectrum_name+"\t No assignment \t \t "
            s=s+"\t"*len( list_of_marker_full_names)+"\n"
            f1.write(s)
            continue
        name_to_peak={(p.marker.code(),p.marker.PTM()):p.mass for p in a.peaks}
        s=a.spectrum_name+"\t"
        if taxonomy:
            if a.lca!=None:
                s=s+str(a.lca) + " ["+ a.lca_name+ "]\t" + a.lca_rank +"\t"+str(a.hca)+ " ["+ a.hca_name + "]\t" + a.hca_rank +"\t"
            else:
                s=s+" None\t\t\t\t"
        for t in a.taxa: 
            s=s+str(t.id)+" ["+t.name+"] "
        s=s+"\t {:.2e}".format(a.pvalue)+"\t"+str(a.score)
        for (name,ptm) in list_of_marker_full_names:
            s=s+"\t"
            if (name,ptm) in  name_to_peak:
                s=s+str(round(name_to_peak[(name,ptm)],3))
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
            s=sp+"\t"+image(peak.marker.code(),peak.marker.PTM())+"\t"+ str(peak.mass)+"\t"+str(peak.intensity)+"\t"
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
    minimum_number_of_peaks=int(config.parse_config_file()["minimum_number_of_peaks"])
    for spectrum in list_of_spectra:
        list_of_assignments.extend(assign_spectrum(spectrum, mass_markers_list, set_of_markers, error, B, threshold, allsolutions, minimum_number_of_peaks))

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
    #create_spectral_file(output, list_of_assignments, list_of_spectra)
    create_json_result_file(output, list_of_assignments)
    
