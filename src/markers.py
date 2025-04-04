"""
   markers.py                             
"""

import copy
import sys
import re
import sys
from  src import utils as ut
from src import sequences
from src import compute_masses
from src import message
#from src import limit 

    
class Marker(object):
    def __init__(self, field={}):
        self.field = field

    def __str__(self):
        return str(self.field)
        
    def sequence(self):
        if "Sequence" in self.field:
            return self.field["Sequence"]
        else:
            return None
            
    def taxid(self):
        if "OX" in self.field:
            return self.field["OX"]
        else:
            return None
            
    def code(self):
        if "Marker" in self.field:
            return self.field["Marker"]
        else:
            return None
        
    def taxon_name(self):
        if "OS" in self.field:
            return self.field["OS"]
        else:
            return None

    def PTM(self):
        return self.field.get("PTM")

    def mass(self):
        if "Mass" not in self.field or float(self.field["Mass"])==0.0:
            return None
        else:
            return float(self.field["Mass"])
        
    def comment(self):
        if "Comment" in self.field:
            return self.field["Comment"]
        else:
            return ""
            
    def begin(self):
        return self.field.get("Begin")
        
    def end(self):
        return self.field.get("End")
        
    def helical(self):
        return self.field.get("Hel")
    
    def status(self):
        return self.field.get("Status")
            
    def protein(self):
        if "GN" in self.field:
            return self.field["GN"]
        else:
            return None
           
    def rank(self):
        return self.field.get("Rank")
     
    def seqid(self):
        return self.field.get("SeqID")
        
    def length(self):
        return self.field.get("Length")
    
    
def sort_by_masses(set_of_markers):
    mass_to_marker_list=[(m.mass(),m) for m in set_of_markers if m.mass() is not None]
    mass_to_marker_list.sort(key=lambda x: x[0])
    mass_markers_list=[]
    current_element=mass_to_marker_list[0]
    s={current_element[1]}
    for i in range(1,len(mass_to_marker_list)):
        new_element=mass_to_marker_list[i]
        if (new_element[0]==current_element[0]): # identical masses
            s.add(new_element[1])
        else: #new mass
            mass_markers_list.append((current_element[0],s))
            current_element=new_element
            s={new_element[1]}
    mass_markers_list.append((current_element[0],s))
    return mass_markers_list
    
def create_marker_landscape(marker_file_name, set_of_markers):
    """
    main file, with peak masses (one line per spectrum)
    """
    name_to_mass_dict={}
    mass_to_taxid_name_dict={}
    list_of_taxid=list({s.taxid() for s in set_of_markers})
    list_of_taxid.sort()
    list_of_codes=list({s.code() for s in set_of_markers})
    list_of_codes.sort()
    for s in set_of_markers:
        ut.update_dictoset(mass_to_taxid_name_dict, round(float(s.mass()),4), {(s.taxid(), s.code())})
        ut.update_dictoset(name_to_mass_dict, s.code, {round(float(s.mass()),4)})
                
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
        ut.update_dictoset(dict1, m.taxid(), {m.sequence()})
    for m in set_of_markers2:
        ut.update_dictoset(dict2, m.taxid(), {m.sequence()})
    for taxid in dict1:
        if taxid not in dict2:
            dict_result[taxid]=dict1[taxid]
        else:
            diff_set=dict1[taxid].difference(dict2[taxid])
            if  len(diff_set)>0:
                dict_result[taxid]=diff_set
    return dict_result
            

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
        dict_taxid_to_name[m.taxid()]= m.taxon_name()
        ut.update_dictoset(dict_taxid_to_masses, m.taxid(), {m.mass()})
    list_of_taxid=list(dict_taxid_to_masses.keys())

    dict_count={}
    for m in set_of_markers:
        ut.increment_dictoset(dict_count, (m.protein(), m.taxid()))

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
 


def remove_lost_taxid(set_of_markers, lost_taxid):
    set_of_markers={m for m in set_of_markers if m.taxid() not in lost_taxid}
    return set_of_markers


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

  
def colinearity(set_of_markers):
    """
    print all markers as a matrix: marker name / position[mass] for the species
    """
    dict_taxon_name_to_taxid={m.taxon_name():m.taxid() for m in set_of_markers}
    print("Total number of species: "+str(len( dict_taxon_name_to_taxid))+"\n")
    list_of_proteins=list({m.protein() for m in set_of_markers if m.protein() is not None})
    list_of_proteins.sort()
    for prot in list_of_proteins:
        print("\nGene : "+prot)
        list_of_seqid=list({(m.taxon_name(), m.seqid()) for m in set_of_markers if m.protein()==prot})
        list_of_seqid.sort(key=lambda x: x[0])
        seqid_length=max({len(ut.pretty_print(s[1])) for s in list_of_seqid})
        taxon_length=max({len(ut.pretty_print(s[0])) for s in list_of_seqid})
        list_of_codes=list({str(m.code())+" - "+str(m.PTM()) for m in set_of_markers if m.protein()==prot})
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
            if m.protein()==prot:
                matrix[list_of_seqid.index((m.taxon_name(), m.seqid()))+1][list_of_codes.index(str(m.code())+" - "+str(m.PTM()))+2]=str(m.begin())
                if m.mass()!=None :
                    matrix_mass[list_of_seqid.index((m.taxon_name(), m.seqid()))+1][list_of_codes.index(str(m.code())+" - "+str(m.PTM()))+2]=str(round(float(m.mass()),4))
        s=" "*taxon_length+"\t"+" "*seqid_length
        for code in list_of_codes:
            s=s+"\t"+"{:<15}".format(code)
        print(s)
        for i in range(len(list_of_seqid)):
            s="{:<{}}".format(ut.pretty_print(list_of_seqid[i][0]), taxon_length)+"\t"+"{:<{}}".format(ut.pretty_print(list_of_seqid[i][1]), seqid_length)
            for j in range(len(list_of_codes)):
                if len(matrix_mass[i+1][j+2])>0 : 
                    s=s+"\t"+"{:<15}".format(ut.pretty_print(matrix[i+1][j+2])+"["+ut.pretty_print(matrix_mass[i+1][j+2])+"]")
                else:
                    s=s+"\t"
            print(s)
    print("")


def str_union(s1, s2):
    if s1 is None:
        return s2
    if s2 is None:
        return s1
    s=set(s1.split()) | set(s2.split())
    return ' '.join(s)


def sort_and_merge(set_of_markers):
    # merge all markers having the same taxid, sequence, PTM  and  mass into a single marker. The new seqId is the union of all taxid. What about the comment ?
    if len(set_of_markers)<2:
        return list(set_of_markers)
    list_of_markers=list(set_of_markers)
    # TO DO: deal with empty mass 
    list_of_markers.sort(key= lambda m: (ut.none_str(m.taxid()),ut.none_str(m.sequence()),ut.none_str(m.PTM()),ut.none_float(m.mass()),ut.none_int(m.begin()),ut.none_str(m.seqid())))
    list_of_selected_markers=[]
    ref_marker=list_of_markers[0]
    for m in list_of_markers[1:]:
        if (ref_marker.taxid(), ref_marker.sequence(), ref_marker.PTM(), ref_marker.mass(),ref_marker.begin())==(m.taxid(), m.sequence(),m.PTM(), m.mass(), m.begin()):
            ref_marker.field["SeqID"]=str_union(ref_marker.seqid(),m.seqid())
        else:
            list_of_selected_markers.append(ref_marker)
            ref_marker=m
    list_of_selected_markers.append(ref_marker)
    return list_of_selected_markers
   

# compute the set of sequences that are compatibles with
# the seqid, the taxid, the taxon name or the protein of the marker m
def find_matching_sequences(m, set_of_sequences):
    matching_sequences=set_of_sequences
    sequence_fields={"OS", "OX", "SeqID", "GN"}
    for field in sequence_fields:
        if field in m.field:
            matching_sequences= {s for s in matching_sequences if ut.equiv(s.field[field],m.field[field])}
    return matching_sequences
    
def update_comment(m, comment):
    m.field["Comment"]=comment+ m.comment()
    return m
        
def post_comment(m, comment):
    m.field["Comment"]= m.comment() + comment
    return m

def supplement_marker(m, set_of_sequences):
    # This case is treated in find_sequence_from_mass
    if m.sequence() is None and m.begin() is None and m.helical() is None:
        return {m}
    matching_sequences=find_matching_sequences(m,set_of_sequences)
    # There is no matching sequence
    if len(matching_sequences)==0:
        message.warning("No protein sequence found for the marker " + str(m) + ".")
        m=update_comment(m, "No matching protein sequence. ")
        return {m}
    position=((m.length()) and (m.begin() or m.end() or m.helical())) or (m.begin() and m.end())
    if (m.seqid(), m.taxid(), m.taxon_name())!=(None, None, None):
        seq=next(iter(matching_sequences))
        m.field["OX"]=seq.taxid()
        m.field["OS"]=seq.taxon_name()
    if m.sequence() :
        return find_positions_from_sequence(m, matching_sequences)
    elif position :
        return find_sequence_from_positions(m, matching_sequences)
    else:
        return ({m})
    
def add_sequences_and_positions_to_markers(set_of_markers, set_of_sequences):
    set_of_new_markers=set()
    for m in set_of_markers:
        set_of_new_markers.update(supplement_marker(m, set_of_sequences))
    return set_of_new_markers

        
# use sequence fo find positions (m.sequence is not None)
def find_positions_from_sequence(m, matching_sequences):
    set_of_markers=set()
    if m.rank() is None:
        m.field["Rank"]="species"
    for seq in matching_sequences:
        pos=(seq.sequence()).find(m.sequence())
        if pos>-1:
            if m.length() is not None and m.length() !=len(m.sequence()):
                message.warning("Peptide length modified in marker "+str(m)+". ")
            if m.begin() is not None and m.begin()!=pos+1 :
                message.warning("Begin position modified in marker "+str(m)+". ")
            if m.end() is not None and m.end()!=pos + len(seq.sequence()) :
                message.warning("End position modified in marker "+str(m)+". ")
            helix=pos + 2 - sequences.helical_region(seq)[0]
            if m.helical() is not None and m.helical()!=helix:
                message.warning("Helical position modified in marker "+str(m)+". ")
            dict={x:m.field[x] for x in m.field}
            dict["OX"]=seq.taxid()
            dict["OS"]=seq.taxon_name()
            dict["GN"]=seq.protein()
            dict["Hel"]= pos - sequences.helical_region(seq)[0] +2
            dict["SeqID"]= seq.seqid()
            dict["Begin"]=pos+1
            dict["End"]= pos + len(m.sequence())
            dict["Length"]=len(m.sequence())
            m2=Marker(field=dict)
            m2=update_comment(m2, "Positions computed from peptide sequence. ")
            set_of_markers.add(m2)
    if len(set_of_markers)==0:
        m=update_comment(m, "No position found for peptide marker. ")
        message.warning("No position found for "+ m.sequence() +" in marker "+str(m))
        return {m}
    else:
        return set_of_markers

# find peptide sequence from its positions
def find_sequence_from_positions(m, matching_sequences):
    set_of_markers=set()
    for seq in matching_sequences:
        if m.begin() and m.length():
            begin=m.begin()
            length=m.length()
            end=m.begin()+m.length()-1
            helical=m.begin() - sequences.helical_region(seq)[0] + 1
        if m.begin() and m.end():
            begin=m.begin()
            end=m.end()
            length=m.end()-m.begin()+1
            helical=m.begin() - sequences.helical_region(seq)[0] + 1
        if m.helical() and m.length():
            length=m.length()
            helical=m.helical()
            begin=sequences.helical_region(seq)[0] +m.helical() -1
            end=begin+m.length()-1
        peptide_sequence=seq.sequence()[begin-1:end]
        dict={x:m.field[x] for x in m.field}
        dict["OX"]=seq.taxid()
        dict["OS"]=seq.taxon_name()
        dict["Sequence"]=peptide_sequence
        dict["GN"]=seq.protein()
        dict["Hel"]=helical
        dict["Length"]=length
        dict["SeqID"]=seq.seqid()
        dict["Begin"]=begin
        dict["End"]=end
        if m.status() is None:
            dict["Status"]="Genetic"
        m2=Marker(field=dict)
        m2=update_comment(m2,"Sequence deduced from positions. ")
        set_of_markers.add(m2)
    return set_of_markers
        
def check_masses_and_sequences(set_of_markers, resolution):
    for m in set_of_markers:
        if m.sequence() is None or m.mass() is None:
            continue
        m2=Marker(field={"Sequence":m.sequence()})
        candidate_markers=compute_masses.add_PTM_or_masses_to_markers({m2}, True, True)
        Found=False
        for cm in candidate_markers:
            if ut.matching_masses(cm.mass(), m.mass(), resolution):
                Found=True
                m.field["mass"]=cm.mass()
                m.field["PTM"]=cm.PTM()
                m=post_comment(m, " Adjusted mass (was "+str(m.mass())+").")
        if not Found:
            m=post_comment(m, " Inconsistent mass.")
    return set_of_markers

# use mass to find sequence
def find_sequences_from_mass(set_of_markers, set_of_sequences, resolution):
    set_of_target_markers={m for m in set_of_markers if (m.taxid() or m.taxon_name()) and m.mass() and m.sequence() is None}  # missing sequences for the set of markers
    target_mass_list=[(m.mass(),m) for m in set_of_target_markers]
    target_mass_list.sort(key=lambda x:x[0])
    set_of_target_sequences=set()
    for m in set_of_target_markers:
        set_of_target_sequences.update(find_matching_sequences(m, set_of_sequences))

    set_of_new_markers={m for m in set_of_markers if m.sequence()}
    set_of_denovo_markers = compute_masses.add_PTM_or_masses_to_markers(sequences.in_silico_digestion(set_of_target_sequences, 2, 10), True, True)
    denovo_mass_list=[(m.mass(),m) for m in set_of_denovo_markers] # masses of all tryptic peptides
    denovo_mass_list.sort(key=lambda x:x[0])
    founded_masses={y:0 for (x,y) in target_mass_list}  
    #for target_mass_list
    for mass in target_mass_list:
        for j in range(0,len(denovo_mass_list)):
            if ut.matching_masses(mass[0], denovo_mass_list[j][0], resolution) and (mass[1].taxid()==denovo_mass_list[j][1].taxid() or ut.equiv(mass[1].taxon_name(), denovo_mass_list[j][1].taxon_name())): # on a un match
                m2=mass[1]
                m= denovo_mass_list[j][1]
                founded_masses[m2]+=1
                dict={x:m.field[x] for x in m.field}
                if m2.code() is None:
                    dict["Marker"]=mass[1].code()
                else:
                    dict["Marker"]=m2.code()
                dict["Comment"]= m.comment() + "Sequence deduced from target mass. "
                new_marker=Marker(field=dict)
                set_of_new_markers.add(new_marker)
                continue
            if denovo_mass_list[j][0] - mass[0] >1:
                continue
        
    for m in founded_masses: 
        if founded_masses[m]==0:
            m=update_comment(m, "No matching peptide found. ")
            set_of_new_markers.add(m)
            message.warning(m.taxon_name()+": no matching peptide found for "+str(m.mass())+".")
        elif founded_masses[m]>1:
            message.warning(m.taxon_name()+": multiple peptides found for "+str(m.mass())+".")
                
    return set_of_new_markers


