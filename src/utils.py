"""
   utils.py
"""

def is_aa_sequence(sequence):
    """ test wether the sequence is an amino-acid sequence """
    if sequence==None or len(sequence)==0 :
        return False
    return all(aa in 'ACDEFGHIKLMNPQRSTVWY' for aa in sequence)

def matching_masses(theoretical_peak, experimental_peak, resolution):
    delta= abs(float(theoretical_peak) - float(experimental_peak))
    if resolution<2:
        return delta<=resolution #dalton
    else:
        return delta/theoretical_peak<=resolution/1000000 #ppm

def is_PTM(PTM_string, set_of_PTM):
    """ test wether the PTM expression is valid """
    if PTM_string==None:
        return True 
    found_number=""
    for char in PTM_string:
        if char.isdigit():
            found_number += char
        else:
          if len(found_number)==0 or (int(found_number)>0 and char not in set_of_PTM):
            return False
          found_number=""
    
    return True
    
def image(set_of_masses):
    result=""
    for mass in set_of_masses:
        result=result+" "+str(mass)
    return result




def update_dictoset(mydict, k,v):
    """
    mydict is a dictionnary whose values are sets (of elements).
    k is a key (new or not) and v is a set of elements that should all be added
    to the set of values attached to k
    """
    if k in mydict:
         mydict[k]=mydict[k].union(v)
    else:
        mydict[k]=v
        

def increment_dictoset(mydict, k):
    if k in mydict:
        mydict[k]=mydict[k]+1
    else:
        mydict[k]=1
                
 
def create_dual_list(mydict):
    """
    create a dual list for dictionary mydict.
    Values of mydict can be sets of IDs (pepids or shuffids)
    """
    l=set()
    for v in mydict.values():
        l.update(v)
    dual_list=[]
    for i in range(len(l)):
        dual_list.append(set())
    for k in mydict:
        for v in mydict[k]:
           dual_list[v].add(k)
    return dual_list

def create_dual_dict(mydict):
    dual_dict={}
    for k in mydict:
        for v in mydict[k]:
            update_dictoset(dual_dict, v,{k})
    return dual_dict
    

def create_mass_Xid_list_from_dict(Xid_to_mass_dict):
    """
    create a list sorted by masses 
    """
    tmp_mass_list=[] 
    for x, set_of_masses in Xid_to_mass_dict.items(): 
        for m in set_of_masses:
            tmp_mass_list.append((m,x))
    tmp_mass_list.sort(key=lambda x: x[0])
    mass_Xid_list=[]
    current_element=tmp_mass_list[0]
    s={current_element[1]}
    for i in range(1,len(tmp_mass_list)):
        new_element=tmp_mass_list[i]
        if (new_element[0]==current_element[0]): # identical masses
            s.add(new_element[1])
        else: #new mass
            mass_Xid_list.append((current_element[0],s))
            current_element=new_element
            s={new_element[1]}
    mass_Xid_list.append((current_element[0],s))
    return mass_Xid_list


def pretty_print_set(s):
    res_str=""
    for e in s:
        res_str=res_str+" "+str(e)
    return res_str
