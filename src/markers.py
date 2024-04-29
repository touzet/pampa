"""
   markers.py                             
"""

import copy
import sys
import re
import sys
from  src import utils as ut
from src import sequences
#from src import limit 

def pretty_print(s):
    if s==None:
        return ""
    else:
        return s
    
class Marker(object):

    def __init__(self, rank=None, taxid=None, taxon_name=None, sequence=None, ptm=None, code=None, mass=None, protein=None, helical=None, seqid=None, begin=None, end=None, comment=None):
        self.rank = rank
        self.taxid = taxid
        self.taxon_name = taxon_name
        self.sequence = sequence
        self.ptm = ptm
        self.code = code
        self.mass = mass
        self.protein = protein
        self.helical = helical
        self.seqid = seqid
        self.begin = begin #1-based notation
        self.end = end  #1-based notation
        self.comment = comment

    def __str__(self):
        return f"{pretty_print(self.rank)}\t{pretty_print(self.taxid)}\t{pretty_print(self.taxon_name)}\t{pretty_print(self.sequence)}\t{pretty_print(self.ptm)}\t{pretty_print(self.code)}\t{pretty_print(self.mass)}\t{pretty_print(self.protein)}\t{pretty_print(self.helical)}\t{pretty_print(self.seqid)}\t{pretty_print(self.begin)}\t{pretty_print(self.end)}\t{pretty_print(self.comment)}"

    
def sort_by_masses(set_of_markers):
    taxid_name_ptm_to_mass_dict={}
    for m in set_of_markers:
        ut.update_dictoset(taxid_name_ptm_to_mass_dict, (m.taxid,m.code,m.ptm), {float(m.mass)})
    mass_taxid_name_list=ut.create_mass_Xid_list_from_dict(taxid_name_ptm_to_mass_dict)
    return mass_taxid_name_list
    
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
        ut.update_dictoset(mass_to_taxid_name_dict, round(float(s.mass),4), {(s.taxid, s.code)})
        ut.update_dictoset(name_to_mass_dict, s.code, {round(float(s.mass),4)})
                
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
        ut.update_dictoset(dict1, m.taxid, {m.sequence})
    for m in set_of_markers2:
        ut.update_dictoset(dict2, m.taxid, {m.sequence})
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
    DEPRECATED
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
        if len(taxid)!=0 and len(taxon_name)!=0:
            ut.update_dictoset(dict_taxid, taxid, {taxon_name})
            ut.update_dictoset(dict_taxonname, taxon_name, {taxid})
        if len(m.mass)>0 and len(sequence)>0:
                ut.update_dictoset(dict_sequence, (sequence, ptm), {m.mass})
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
        ut.update_dictoset(dict_taxid_to_masses, m.taxid, {m.mass})
    list_of_taxid=list(dict_taxid_to_masses.keys())

    dict_count={}
    for m in set_of_markers:
        ut.increment_dictoset(dict_count, (m.protein, m.taxid))

    dict_t={}
    for taxid in list_of_taxid:
        s={(x[0],dict_count[(x[0],taxid)]) for x in dict_count if x[1]==taxid}
        dict_t[taxid]=s

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
            print("  - ", end="")
            for t in dict_set[taxid]:
                print(dict_taxid_to_name[t], end="   ")
            print("")
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
        ut.update_dictoset(dict_of_markers, m.sequence.strip(), {model_marker})
    return dict_of_markers


def remove_lost_taxid(set_of_markers, lost_taxid):
    set_of_markers={m for m in set_of_markers if m.taxid not in lost_taxid}
    return set_of_markers

def reduce(seq):
    seq=seq.replace(" ", "")
    seq=seq.lower()
    return  seq


def authorized_PTM(PTM_string, list_of_PTM):
    found_number=""
    for char in PTM_string:
        if char.isdigit():
            found_number += char
        else:
          if int(found_number)>0 and char not in list_of_PTM:
            return False
          found_number=""
    
    return True

  
def limit_markers(set_of_markers, limit_file):
    """
    build a new set of markers that is compliant with the specification of the limit_file (txt file, -l option) 
    """
    ## DEPRECATED ?
    new_set=set_of_markers
    file=open(limit_file).read().splitlines()
    for line in file:
        # each line defines a set of markers
        # the final set of markers is the union of all lines
        if line[:3]=="GN=":
            pattern = re.compile(r'GN=([\w\s,]+)')
            list_of_matches = pattern.findall(line)
            matches= list_of_matches[0].split(',')
            new_set={s for s in new_set if s.protein in matches}

        if line[:3]=="OS=":
            pattern = re.compile(r'OS=([\w\s,]+)')
            list_of_matches = pattern.findall(line)
            matches= list_of_matches[0].split(',')
            matches=[reduce(match) for match in matches]
            new_set={s for s in new_set if reduce(s.taxon_name) in matches}   

        if line[:6]=="SeqID=":
            pattern = re.compile(r'SeqID=([\w\s,]+)')
            list_of_matches = pattern.findall(line)
            matches= list_of_matches[0].split(',')
            matches=[reduce(match) for match in matches]
            new_set={s for s in new_set if reduce(s.seqid) in matches} 

        if line[:4]=="PTM=":
            pattern = re.compile(r'PTM=([\w\s,]+)')
            list_of_matches = pattern.findall(line)
            matches= list_of_matches[0].split(',')
            #matches=[reduce(match) for match in matches]
            new_set={s for s in new_set if authorized_PTM(s.ptm, matches)}
            
    summary_file.close()
    return new_set

def colinearity(set_of_markers):
    """
    print all markers as a matrix: marker name / position[mass] for the species
    """
    dict_taxon_name_to_taxid={m.taxon_name:m.taxid for m in set_of_markers}
    print("Total number of species: "+str(len( dict_taxon_name_to_taxid))+"\n")
    list_of_proteins=list({m.protein for m in set_of_markers if m.protein is not None})
    list_of_proteins.sort()
    for prot in list_of_proteins:
        print("\nGene : "+prot)
        list_of_seqid=list({(m.taxon_name, m.seqid) for m in set_of_markers if m.protein==prot})
        list_of_seqid.sort(key=lambda x: x[0])
        seqid_length=max({len(s[1]) for s in list_of_seqid})
        taxon_length=max({len(s[0]) for s in list_of_seqid})
        list_of_codes=list({str(m.code)+" - "+str(m.ptm) for m in set_of_markers if m.protein==prot})
        list_of_codes.sort()
        matrix = [["" for j in range(len(list_of_codes)+2)] for i in range(len(list_of_seqid)+1)]
        matrix_mass=[["" for j in range(len(list_of_codes)+2)] for i in range(len(list_of_seqid)+1)]
        for code in list_of_codes:
             matrix[0]=["",""]+list_of_codes
             matrix_mass[0]=["",""]+list_of_codes
        for (i,seq) in enumerate(list_of_seqid):
            matrix[i+1][0]= seq[0]
            matrix[i+1][1]= seq[1]
        for m in set_of_markers:
            if m.protein==prot:
                matrix[list_of_seqid.index((m.taxon_name, m.seqid))+1][list_of_codes.index(str(m.code)+" - "+str(m.ptm))+2]=str(m.begin)
                if m.mass!=None and len(str(m.mass))!=0:
                    matrix_mass[list_of_seqid.index((m.taxon_name, m.seqid))+1][list_of_codes.index(str(m.code)+" - "+str(m.ptm))+2]=str(round(float(m.mass),4))
        s=" "*taxon_length+"\t"+" "*seqid_length
        for code in list_of_codes:
            s=s+"\t"+"{:<15}".format(code)
        print(s)
        for i in range(len(list_of_seqid)):
            s="{:<{}}".format(list_of_seqid[i][0], taxon_length)+"\t"+"{:<{}}".format(list_of_seqid[i][1], seqid_length)
            for j in range(len(list_of_codes)):
                if len(matrix_mass[i+1][j+2])>0 : 
                    s=s+"\t"+"{:<15}".format(matrix[i+1][j+2]+"["+matrix_mass[i+1][j+2]+"]")
                else:
                    s=s+"\t"
            print(s)
    print("")


def sort_and_merge(set_of_markers):
    # merge all markers having the same taxid, sequence, PTM  and  mass into a single marker. The new seqId is the union of all taxid. What about the comment ?

    list_of_markers=[]
    for m in set_of_markers:
        m.taxid=pretty_print(m.taxid)
        m.code=pretty_print(m.code)
        m.sequence=pretty_print(m.sequence)
        m.ptm=pretty_print(m.ptm)
        m.seqid=pretty_print(m.seqid)
        if m.mass==None:
            m.mass=0
        list_of_markers.append(m)
    # TO DO: deal with empty mass 
    list_of_markers.sort(key= lambda m: (m.taxid,m.sequence,m.ptm,float(m.mass),m.seqid))
    for (i,m) in enumerate(list_of_markers):
        j=i+1
        while j<len(list_of_markers) and (list_of_markers[j].taxid, list_of_markers[j].sequence, list_of_markers[j].ptm,float(m.mass))==(m.taxid, m.sequence,m.ptm,float(m.mass)): 
            m.seqid=m.seqid+" "+list_of_markers[j].seqid
            if m.comment!=list_of_markers[j].comment:
                m.comment=m.comment+", "+list_of_markers[j].comment
            list_of_markers[j].taxid=None
            j=j+1

    list_of_markers=[m for m in list_of_markers if m.taxid is not None]
        
    return list_of_markers
        
def add_sequences_and_positions_to_markers(set_of_markers, set_of_sequences):
    set_of_marker_seqid={m.seqid for m in set_of_markers}
    dict_of_sequence_seqid={id:s for id in set_of_marker_seqid for s in set_of_sequences if s.seqid==id}
    set_of_missing_seqid=set_of_marker_seqid-dict_of_sequence_seqid.keys()
    if len(set_of_missing_seqid)>0:
        message.warning("No sequence found for SeqID "+str(set_of_missing_seqid)+". Ignored")
    for m in set_of_markers:
        if m.seqid in set_of_missing_seqid:
                continue
        s = dict_of_sequence_seqid[m.seqid]
        seq = s.sequence
        helical_start=sequences.helical_region(seq)[0]
        if m.sequence==None and m.begin and m.end:
            m.sequence = seq[int(m.begin)-1:int(m.end)+1]
        if m.begin==None:
            pos=seq.find(m.sequence)
            if pos > -1:
                m.begin=pos + 1
                if m.end==None:
                    m.end = pos+len(m.sequence)
        if m.helical==None:
            m.helical = int(m.begin) - helical_start + 1
               
    return set_of_markers
         
 
            

        
