"""
taxonomy.py                      

"""

import sys
#sys.path.append("src")
from src import utils 

class Taxonomy(object):
    def __init__(self, taxidname={}, rank={}, children={}, descendants={}, parent={},  root="0"):
        self.taxidname=taxidname
        # key (str): taxid 
        # value (str): name of the taxid
        self.rank=rank
         # key (str): taxid 
        # value (str): rank of the taxid
        self.children=children
        #key (str): taxid
        #value (set of str): set of children of the node taxid
        #   /!\  leaf nodes are not in this dict
        self.descendants=descendants
        #key (str): taxid 
        #value (set of str): set of descendants of the node taxid. Children are identified  by their taxid 
        self.parent=parent
        #key(str): taxid (str)
        #value (str): taxid of the parent of the node taxid
        self.root=root
        #taxid of the root (str)

    def __contains__(self, taxid):
        return taxid in self.taxidname.keys()

    def __len__(self):
        return (len( self.taxidname.keys()))
    
    def __iter__(self):
        return iter([x for x in self.taxidname.keys()])

    def is_leaf(self,taxid):
        return  (taxid not in self.children)
        
    def number_of_children(self, taxid):
        if taxid not in self.children:
            return 0
        else:
            return len(self.children[taxid])

    def init_parent(self):
         for taxid in self.children.keys():
             for t in  self.children[taxid]:
                self.parent[t]=taxid
        
    def init_root(self):     
        s=self.children.keys()
        if len(s)==0:
            self.root="0"
            self.children["0"]={c for c in self.rank.keys()}
            self.taxidname["0"]="None"
            self.rank["0"]=""
            return
        for taxid in self.children.keys():
            s=s - self.children[taxid]
        if len(s)==1:
            r=s.pop()
        else:
            r="0"
            self.children[r]=s
        self.root=r
        self.taxidname[r]="root"
        self.rank[r]="ROOT"
        

    def set_of_descendants_aux(self, taxid):      
        if self.number_of_children(taxid)==0: 
            self.descendants[taxid]={taxid}
        else:
            s={taxid}
            for t in self.children[taxid]:
                self.set_of_descendants_aux(t)
                s=s.union(self.descendants[t])
            self.descendants[taxid]=s
            
    def init_descendants(self):
          self.set_of_descendants_aux(self.root)

    def intersection(self, set_of_taxid):
        """ create a sub-taxonomy that contains only taxid from set_of_taxid, with their ancestors"""
        """ elements of set_of_taxid not present in the taxonomy are lost. """
        set_of_survivors=set()
        lost_taxid=set()
        for taxid in set_of_taxid:
            if taxid not in self.parent:
                lost_taxid.add(taxid)
            else:
                t=taxid
                while t in self.parent:
                    t=self.parent[t]
                    print(t)
                    set_of_survivors.add(t)
        set_of_survivors.update(set_of_taxid-lost_taxid)
        
        taxid_to_children_dict= {taxid: self.children[taxid] & set_of_survivors for taxid in set_of_survivors & self.children.keys()}
        taxid_to_name_dict = {taxid: self.taxidname[taxid]  for taxid in set_of_survivors} 
        taxid_to_rank_dict = {taxid: self.rank[taxid]  for taxid in set_of_survivors}

        B=Taxonomy()
        B.children = taxid_to_children_dict
        B.taxidname = taxid_to_name_dict
        B.rank=taxid_to_rank_dict
        B.init_root()
        B.init_descendants()
        B.init_parent()
        return B, lost_taxid

  
    def intersection_with_descendants(self, set_of_taxid):
        """ create a sub-taxonomy that contains all taxid from set_of_taxid, with their ancestors AND descendants """
        set_of_survivors={ m for m in set_of_taxid}
        for taxid in set_of_taxid :
            if self.number_of_children(taxid)>0:
                set_of_survivors=set_of_survivors.union(self.descendants[taxid])
            t=taxid
            while t in self.parent:
                t=self.parent[t]
                set_of_survivors.add(t)
        taxid_to_children_dict= {taxid: self.children[taxid] & set_of_survivors for taxid in set_of_survivors & self.children.keys()}
        taxid_to_name_dict = {taxid: self.taxidname[taxid]  for taxid in set_of_survivors} 
        taxid_to_rank_dict = {taxid: self.rank[taxid]  for taxid in set_of_survivors}

        B=Taxonomy()
        B.children = taxid_to_children_dict
        B.taxidname = taxid_to_name_dict
        B.rank=taxid_to_rank_dict
        B.init_root()
        B.init_descendants()
        B.init_parent()
        return B
    
    def lca(self, set_of_taxid):
        """ lowest common ancestor """
        if not set_of_taxid.issubset(self):
            return self.root
        ancestor=next(iter(set_of_taxid))
        while not set_of_taxid.issubset(self.descendants[ancestor]) and not ancestor==self.root :
            ancestor=self.parent[ancestor]
        return ancestor

    def hca(self, set_of_taxid):
        """ highest common ancestor: highest node whose descendants are set_of_taxid """
        hca=self.lca(set_of_taxid)
        if hca==self.root:
            return hca
        next_hca=self.parent[hca]
        while not next_hca==self.root and self.number_of_children(next_hca)==1:
            hca=next_hca
            next_hca=self.parent[next_hca]
        return hca
        
## end of class Taxonomy ##


def table_print_rec(t, taxid, rank):
    print("."*rank + " ["+t.rank[taxid]+"] "+ taxid +" "+ t.taxidname[taxid])
    if taxid not in t.children:
        return
    for c in t.children[taxid]:
        table_print_rec(t, c, rank+1)
        
def table_print(t, taxid=""):
    if not taxid:
        taxid=t.root
    table_print_rec(t,taxid,0)
     

def parse_taxonomy_simple_file(taxonomy_file):
    """
    taxonomy_file is a TSV file with 5 columns:
    Taxid | Common name	| Scientific name | Parent | Rank
    """
    taxonomy=Taxonomy()
    name_to_taxid_dict={} # key: (name, rank)
    taxid_to_children_dict={}
    taxid_to_name_dict={} 
    taxid_to_rank_dict={}
    taxid_to_parent_dict={}
    
    with open(taxonomy_file) as in_file:
        next(in_file)
        for line in in_file:
            columns = line.split("\t")
            taxid_to_name_dict.update({columns[0]:columns[2]})
            name_to_taxid_dict.update({columns[2]:columns[0]})
            taxid_to_rank_dict.update({columns[0]:columns[4].strip("\n")})
            taxid_to_parent_dict.update({columns[0]:columns[3]})
            
    for taxid in taxid_to_parent_dict.keys():
        utils.update_dictoset(taxid_to_children_dict,taxid_to_parent_dict[taxid],{taxid})

    taxonomy.children=taxid_to_children_dict
    taxonomy.taxidname=taxid_to_name_dict
    taxonomy.rank=taxid_to_rank_dict
    taxonomy.parent=taxid_to_parent_dict
    taxonomy.init_root()
    taxonomy.init_descendants()
    
    return taxonomy


 
def build_flat_taxonomy(set_of_markers):
    """ constructs a flat taxonomy (one root, all taxid are leaves) from a set of markers """
    taxidname={}
    rank={}
    for m in set_of_markers:
        taxidname[m.taxid]= m.taxon_name
        rank[m.taxid]=""
    t=Taxonomy(taxidname, rank)
    t.init_root()
    t.init_descendants()
    t.init_parent()
    return t

def search_taxid_from_taxon_name(taxon_name, taxonomy):
     for key, value in taxonomy.taxidname.items():
         if taxon_name == value:
             return key
    
