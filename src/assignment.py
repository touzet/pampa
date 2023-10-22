"""
    assignment.py
"""


from utils import *
import math


class Assignation(object):
    def __init__(self, spectrum_name="",  peak_to_taxid_name_ptm={}, taxid_to_score={}, taxid_to_optimal={}, taxid_to_peak_name_ptm={}, clusters={}, assigned_taxid=set(), score=0, lca=0):
        self.spectrum_name = spectrum_name
        self.peak_to_taxid_name_ptm=peak_to_taxid_name_ptm
        self.taxid_to_score=taxid_to_score
        self.taxid_to_optimal=taxid_to_optimal
        self.taxid_to_peak_name_ptm=taxid_to_peak_name_ptm
        self.clusters=clusters
        self.assigned_taxid = assigned_taxid
        self.score=0
        self.lca=lca


        
def mcc(TP,TN, FP, FN):
    return (TN*TP-FP*FN)/math.sqrt((TN+FN)*(FP+TP)*(TN+FP)*(FN+TP))


def assign_peaks_of_the_spectrum(spectrum, mass_taxid_name_list, resolution):  
    spectrum.sort()  # peaks are sorted according to their mass
    peak_to_taxid_name_ptm={} #key: peak, value: set of triplets (taxid,name,ptm)  
    current_j = 0
    for peak in spectrum:
        for j in range(current_j,len(mass_taxid_name_list)):
            diff_peak_pep = peak-mass_taxid_name_list[j][0]
            if abs(diff_peak_pep)<=resolution: # there is a match 
                update_dictoset(peak_to_taxid_name_ptm, peak, mass_taxid_name_list[j][1])
            elif diff_peak_pep<0: # peptide mass is gerater than peak mass 
                current_j = j 
                break # nex peak 

    assignation=Assignation()
    assignation.spectrum_name=spectrum.name
    assignation.peak_to_taxid_name_ptm= peak_to_taxid_name_ptm
    
    return assignation 


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

def assign_spectrum(spectrum, mass_taxid_name_list, list_of_taxid, resolution, taxonomy, threshold):
    assignation=assign_peaks_of_the_spectrum(spectrum, mass_taxid_name_list, resolution)
    taxid_to_score={}
    taxid_to_optimal={}
    taxid_to_peak_name_ptm={}
    if len(assignation.peak_to_taxid_name_ptm)==0:
        for taxid in list_of_taxid:
            taxid_to_score[taxid]=-1
        assignation.taxid_to_score=taxid_to_score
        assignation.taxid_to_peak_name_ptm=taxid_to_peak_name_ptm
        assignation.assigned_taxid=set()
        assignation.lca=taxonomy.root #lca of taxid with the best score
        assignation.score=-1
        return assignation

    taxid_to_number_of_masses_dict=count_number_of_masses_per_taxid(list_of_taxid, mass_taxid_name_list)

    for peak in assignation.peak_to_taxid_name_ptm:
        for (taxid,name,ptm) in  assignation.peak_to_taxid_name_ptm[peak]:
            update_dictoset(taxid_to_peak_name_ptm, taxid, {(peak, name,ptm)})

    score_max=-1
    for taxid in taxid_to_peak_name_ptm:
        score=len(taxid_to_peak_name_ptm[taxid])
        if score>=(score_max*threshold)/100:
            taxid_to_score[taxid]=score
            if score>score_max:
                score_max=score
                
    # suppression des espèces de mauvais score
    taxid_to_score={taxid: taxid_to_score[taxid] for taxid in taxid_to_score if taxid_to_score[taxid] >=(score_max*threshold)/100}
    
    # sléection des espèces de score maximum
    set_of_best_species={taxid for taxid in taxid_to_score.keys() if taxid_to_score[taxid]==score_max}
    
    # espèces sous-optimales et clusters
    clusters={} # key: one representative species, value: all species with same peaks
    for taxid in taxid_to_score:
        representant=taxid 
        taxid_to_optimal[taxid]=True
        for taxid2 in taxid_to_score:
            if taxid_to_peak_name_ptm[taxid]<taxid_to_peak_name_ptm[taxid2]:
                taxid_to_optimal[taxid]=False
        representant=taxid
        for taxid2 in clusters: 
            if taxid_to_peak_name_ptm[taxid]==taxid_to_peak_name_ptm[taxid2]:
                representant=taxid2
        update_dictoset(clusters, representant, {taxid})
        
    assignation.clusters=clusters
    assignation.taxid_to_optimal=taxid_to_optimal
    assignation.taxid_to_score=taxid_to_score
    assignation.taxid_to_peak_name_ptm= taxid_to_peak_name_ptm
    assignation.assigned_taxid=set_of_best_species
    if len(set_of_best_species)>0:
        assignation.lca=taxonomy.hca(set_of_best_species)
    else:
        assignation.lca=None
    assignation.score=score_max
    return assignation

def image(name, ptm):
    return name+ " - "+ptm

def assign_all_spectra(list_of_spectra, mass_taxid_name_list, list_of_taxid, list_of_marker_names, resolution, taxonomy, B, path_name, marker_source, threshold, outfile_name):
    set_of_marker_full_names=set()
    for (x,y)  in  mass_taxid_name_list:
        set_of_marker_full_names.update(y)
    list_of_marker_full_names=list({(b,c) for (a,b,c) in set_of_marker_full_names})
    list_of_marker_full_names.sort()
    list_of_assignations=[]
    for spectrum in list_of_spectra:
        current_assignation = assign_spectrum(spectrum, mass_taxid_name_list, list_of_taxid, resolution, B, threshold)
        list_of_assignations.append(current_assignation)
    # selection of all useful taxid
    set_of_useful_taxid=set()
    for a in list_of_assignations:
        set_of_useful_taxid.update((a.taxid_to_score).keys())
    list_of_useful_taxid=list(set_of_useful_taxid)
    list_of_useful_taxid.sort()

    # Primary file - heading
    f1 = open(outfile_name, "w+")
    f1.write("Spectra:\t"+path_name+"\nMarkers:\t"+ marker_source + "\nResolution:\t"+str(resolution)+"\n" )
    s=""
    for (name,ptm) in  list_of_marker_full_names:
        s=s+"\t"+image(name,ptm)
    if taxonomy:
        s=s+"\tScore \t Assignment \t Rank \t Species\n"
    else:
        s=s+"\tScore \t Species\n"
    f1.write(s)

    # Secondary file - heading
    f2 =open("detail_"+outfile_name, "w+")
    f2.write("Spectra:\t"+path_name+"\nMarkers:\t"+ marker_source + "\nResolution:\t"+str(resolution)+"\n" )
    s="\t \t m/z \t intensity  "
    for taxid in  list_of_useful_taxid:
        s=s+"\t"+str(taxid)
    s=s+"\t Best Guess\n"
    f2.write(s)

    
    for i in range(0,len(list_of_spectra)):
        a=list_of_assignations[i]
        sp=list_of_spectra[i]
        mass_to_intensity={p.mass:p.intensity for p in sp.peaks}
      
        # Primary file - body
        for taxid in a.clusters:
            if not a.taxid_to_optimal[taxid]:
                s=taxid+" no\n"
                continue
            name_to_peak={(name,ptm):peak for (peak,name,ptm) in a.taxid_to_peak_name_ptm[taxid]}
            s=a.spectrum_name+"\t"
            for (name,ptm) in list_of_marker_full_names:
                if (name,ptm) in  name_to_peak:
                    s=s+str(name_to_peak[(name,ptm)])
                s=s+"\t"
            if taxonomy:
                if a.lca!=None:
                    s=s+str(a.taxid_to_score[taxid])+"\t"+str(a.lca) + " ["+ B.taxidname[a.lca] + "]\t" + B.rank[a.lca] +"\t"
                else:
                    s=s+"\t None\t\t" 
            else:
                s=s+str(a.taxid_to_score[taxid])+"\t"
            for t2 in a.clusters[taxid]: 
                s=s+str(t2)+" ["+B.taxidname[t2]+"] "
            s=s+"\n"
            f1.write(s)
        
        # Secondary file - body
        name_to_peak={}
        for peak in  a.peak_to_taxid_name_ptm:
            for (taxid, name, ptm) in  a.peak_to_taxid_name_ptm[peak]:
                update_dictoset(name_to_peak, (name,ptm), {peak})
        for (name,ptm) in list_of_marker_full_names:
            if (name,ptm) in name_to_peak :
                # the peptide marker is present in the spectrum
                for peak in name_to_peak[(name,ptm)]:
                    found=False
                    #creation of a new line in the file
                    s=a.spectrum_name+"\t"+image(name,ptm)+"\t"+ str(peak)+"\t"+str(int(mass_to_intensity[peak]))+"\t"
                    for taxid in list_of_useful_taxid:
                        if  taxid in a.taxid_to_score and (taxid,name,ptm) in a.peak_to_taxid_name_ptm[peak] :
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
            if taxid in a.taxid_to_score:
                s=s+"\t"+str(int(a.taxid_to_score[taxid]))
                update_dictoset(score_to_taxid, a.taxid_to_score[taxid], {taxid})
            else:
                s=s+"\t"
        s=s+"\t"
        for score in score_to_taxid:
            s=s+str(score)+" -> "
            for taxid in score_to_taxid[score]:
                s=s+taxid+" ["+B.taxidname[taxid]+"] "
        s=s+"\n"
        f2.write(s)
        
    f1.close()
    f2.close()
    

 
    
