"""
   markers.py                             
"""

import copy
import sys
from utils import *

class Marker(object):

    def __init__(self, rank="", taxid="", taxon_name="", sequence="", ptm="", code="", masses=set(), protein="", seqid="", begin="", end="", comment=""):
        self.rank = rank
        self.taxid = taxid
        self.taxon_name = taxon_name
        self.sequence = sequence
        self.ptm = ptm
        self.code = code
        self.masses = masses
        self.protein = protein
        self.seqid = seqid
        self.begin = begin
        self.end = end
        self.comment = comment

    def __str__(self):
        return f"{self.rank}\t{self.taxid}\t{self.taxon_name.lstrip()}\t{self.sequence}\t{self.ptm}\t{self.code}\t"+image(self.masses)+f"\t{self.protein}\t{self.seqid}\t{self.begin}\t{self.end}\t{self.comment}"

    
def sort_by_masses(set_of_markers):
    taxid_name_ptm_to_mass_dict={}
    for m in set_of_markers:
        update_dictoset(taxid_name_ptm_to_mass_dict, (m.taxid,m.code,m.ptm), m.masses)
    mass_taxid_name_list=create_mass_Xid_list_from_dict(taxid_name_ptm_to_mass_dict)
    return mass_taxid_name_list
    
def build_markers_dictionnaries_from_set_of_markers(set_of_markers):
    set_of_codes=set()
    set_of_taxid=set()
    taxid_to_taxidname={}
    for current_marker in set_of_markers:
        taxid=current_marker.taxid
        taxid_to_taxidname[taxid]=current_marker.taxon_name
        set_of_taxid.add(taxid)
        set_of_codes.add(current_marker.code)
    mass_taxid_name_list=sort_by_masses(set_of_markers)
    return  mass_taxid_name_list, set_of_codes, set_of_taxid, taxid_to_taxidname


def build_markers_pepid_taxid_dictionnaries_from_set_of_markers(set_of_markers):
    taxid_to_pepid_dict={}
    set_of_taxid=set()
    taxid_to_taxidname={}
    number_of_pepid=0
    pepid_to_seq_list=[]
    seq_to_pepid_dict={}
    
    for current_marker in set_of_markers:
        taxid=current_marker.taxid
        taxid_to_taxidname[taxid]=current_marker.taxon_name
        set_of_taxid.add(taxid)
        seq=current_marker.sequence 
        if seq not in seq_to_pepid_dict:
            pepid=number_of_pepid
            number_of_pepid=number_of_pepid + 1
            seq_to_pepid_dict[seq]=pepid
            pepid_to_seq_list.append(seq)
        else:
            pepid=seq_to_pepid_dict[seq]
        if taxid in taxid_to_pepid_dict: # update_dictoset ?
            taxid_to_pepid_dict[taxid].add(pepid)
        else:
            taxid_to_pepid_dict[taxid]={pepid}

    pepid_to_taxid_list=create_dual_list(taxid_to_pepid_dict)
        
    return  pepid_to_taxid_list, pepid_to_seq_list, set_of_taxid, taxid_to_taxidname

def create_marker_landscape(marker_file_name, set_of_markers):
    """
    main file, with peak masses (one line per spectrum)
    """
    name_to_mass_dict={}
    mass_to_taxid_name_dict={}
    list_of_taxid=list({s.taxid for s in set_of_markers})
    list_of_taxid.sort()
    list_of_codes=list({s.code for s in set_of_markers})
    list_of_codes.sort()
    print(list_of_codes)
    for s in set_of_markers:
        for m in s.masses:
            update_dictoset(mass_to_taxid_name_dict, round(float(m),4), {(s.taxid, s.code)})
            update_dictoset(name_to_mass_dict, s.code,  {round(float(m),4)})
                
    f = open("landscape_"+marker_file_name, "w+")
    f.write(marker_file_name+"\n")
    s="\t"
    for taxid in  list_of_taxid:
        s=s+"\t"+taxid
    f.write(s+"\n")
    
    for name in list_of_codes:
        if name in name_to_mass_dict:
            mass_list=list(name_to_mass_dict[name])
            mass_list.sort()
            for mass in mass_list:
                #creation of a new line in the file
                s=name+"\t"+str(mass)+"\t"
                for taxid in list_of_taxid:
                    if (taxid,name) in mass_to_taxid_name_dict[mass]:
                         s=s+"X\t"
                    else:
                        s=s+"\t"
                f.write(s+"\n")
            f.write("- - -\t"*(len(list_of_taxid)+2)+"\n")
    f.close()



def is_included_taxid_sequence(set_of_markers1, set_of_markers2):
    """ 
    look if the set of pairs (taxid, sequence) for set_of_markers1 is included in the set of pairs
    (taxid, sequence) for set_of_markers2. 
    dict_result contains for each taxid the set of sequences present in set_of_markers1  
    and not in set_of_markers2
    """ 
    dict1={}
    dict2={}
    dict_result={}
    for m in set_of_markers1:
        update_dictoset(dict1, m.taxid, {m.sequence})
    for m in set_of_markers2:
        update_dictoset(dict2, m.taxid, {m.sequence})
    for taxid in dict1:
        if taxid not in dict2:
            dict_result[taxid]=dict1[taxid]
        else:
            diff_set=dict1[taxid].difference(dict2[taxid])
            if  len(diff_set)>0:
                dict_result[taxid]=diff_set
    return dict_result
            
    


def is_consistent(set_of_markers):
    """
    check the consistency of set_of_markers:
    one taxid is equivalent to one species and all identical sequences+PTM have the same mass
    """
    dict_taxid={}
    dict_taxonname={}
    dict_sequence={}
    for m in set_of_markers:
        taxon_name=(m.taxon_name).replace(" ","")
        taxon_name=taxon_name.lower()
        taxid=(m.taxid).replace(" ","")
        sequence=(m.sequence).replace(" ","")
        sequence=sequence.lower()
        ptm=(m.ptm).replace(" ","")
        ptm=ptm.replace(" ","")
        if len(taxid)!=0 and len(taxon_name)!=0:
            update_dictoset(dict_taxid, taxid, {taxon_name})
            update_dictoset(dict_taxonname, taxon_name, {taxid})
        if len(m.masses)>0 and len(sequence)>0:
                update_dictoset(dict_sequence, (sequence, ptm), frozenset(m.masses))
    for taxid, value in dict_taxid.items():
        if len(value)>1:
            print(" taxid", taxid, "->", value)
    for taxon, value in dict_taxonname.items():
        if len(value)>1:
            print(" taxon", taxon, "->", value)
    for sequence, value in dict_sequence.items():
        if len(value)>1:
            print(" sequence", sequence, "->", value)  



def check_set_of_markers(set_of_markers, file=""):
    """
    find equivalent species : species that cannot be distinguished because all peptides
    have identical masses.
    TO DO: deal with the case where a marker is in multiple copies (duplicated line)
    """
    dict_taxid_to_masses={}
    dict_taxid_to_name={}

    if len(file)>0:
        sys.stdout = open(file,'a+')
           
    for m in set_of_markers:
        dict_taxid_to_name[m.taxid]= m.taxon_name
        update_dictoset(dict_taxid_to_masses, m.taxid, m.masses)
    list_of_taxid=list(dict_taxid_to_masses.keys())

    dict_count={}
    for m in set_of_markers:
        increment_dictoset(dict_count, (m.protein, m.taxid))

    dict_t={}
    for taxid in list_of_taxid:
        s={(x[0],dict_count[(x[0],taxid)]) for x in dict_count if x[1]==taxid}
        dict_t[taxid]=s

    print("  Number of species : "+str(len(dict_taxid_to_name))+"\n")
        
    set_of_values={frozenset(x) for x in dict_t.values()}
    for v in set_of_values:
        print("  ", end="")
        for (x,y) in v:
            print(str(y)+" "+str(x)+" ", end="")
        print("")
        key_list=[str(t) for t in dict_t if dict_t[t]==v]
        for t in key_list:
            print("  "+dict_taxid_to_name[t]+" ["+str(t)+"]  ")
        print("")
    dict_equiv={}
    is_included=set()
   
    for (i,taxid1) in enumerate(list_of_taxid[0:len(list_of_taxid)-1:]):
        for (j,taxid2) in enumerate(list_of_taxid[i+1::]) :  
            if dict_taxid_to_masses[taxid1] == dict_taxid_to_masses[taxid2] and taxid2 not in dict_equiv:
                dict_equiv[taxid2]=taxid1
            elif  dict_taxid_to_masses[taxid1] < dict_taxid_to_masses[taxid2]:
                is_included.add((taxid1,taxid2))
            elif  dict_taxid_to_masses[taxid2] < dict_taxid_to_masses[taxid1]:
                is_included.add((taxid2,taxid1))
    dict_set={}
    for taxid in dict_equiv.values():
         dict_set[taxid]={t for t in dict_equiv if dict_equiv[t]==taxid}.union({taxid})

    print("---------------------------------")     
    print("  UNDISTINGUISHABLE SPECIES") 
    print("---------------------------------\n")

    existence=False
    for taxid in dict_set:
        if len(dict_set[taxid])>1:
            existence=True
            print("  "+str(dict_set[taxid]))
            #print("  - ", end="")
            #for t in dict_set[taxid]:
            #    print(dict_taxid_to_name[t]+", ", end=" ")        
    if not existence:
        print("  None")
    print("")
    print("---------------------------------")          
    print("  INCLUSION OF SPECIES")
    print("---------------------------------\n")   

    existence=False
    for taxid1 in list_of_taxid:
        s={taxid2 for taxid2 in list_of_taxid if (taxid1, taxid2)   in  is_included}
        if len(s)>0:
            existence=True
            print("  - "+dict_taxid_to_name[taxid1]+" is included in : ", end="")
            for taxid2 in s:
                print(dict_taxid_to_name[taxid2],end="  ")
            print("\n")
    if not existence:
        print("  None\n") 
 

def build_marker_dict_from_set_of_marker(set_of_markers):
    """
    dict_of_markers:
    key: peptide sequence
    value: set of 4-uplets (ptm, code, prot, masses)
    """
    dict_of_markers={}
    marker_count=0
    for m in set_of_markers:
        if len(m.sequence)==0 :
            continue
        model_marker=(m.ptm.strip(), m.code.strip(), (m.protein.strip()).upper(),m.taxid,m.begin)
        update_dictoset(dict_of_markers, m.sequence.strip(), {model_marker})
    return dict_of_markers


