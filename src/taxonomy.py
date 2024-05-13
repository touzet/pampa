"""
taxonomy.py                      

"""

import sys

from src import utils
from src import message

class Taxonomy(object):
    def __init__(self, name={}, rank={}, children={}, descendants={}, parent={},  root=set()):
        self.name=name
        # key (str): taxid 
        # value (str): name of the taxid
        self.rank=rank
         # key (str): taxid 
        # value (str): rank of the taxid
        self.children=children
        #key (str): taxid
        #value (set of str): set of children of the node taxid
        #   /!\  leaf nodes are not in this dictionary
        self.descendants=descendants
        #key (str): taxid 
        #value (set of str): set of descendants of the node taxid. Children are identified  by their taxid 
        self.parent=parent
        #key(str): taxid (str)
        #value (str): taxid of the parent of the node taxid
        self.root=root
        #taxid of the root (str)

    def __contains__(self, taxid):
        return taxid in self.name.keys()

    def __len__(self):
        return (len( self.name.keys()))
    
    def __iter__(self):
        return iter([x for x in self.name.keys()])

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
        s=self.children.keys() # set of internal nodes
        if len(s)==0: # flat taxonomy
            self.root=self.name.keys()
            return
        for taxid in self.children.keys():
            s=s - self.children[taxid]
        self.root=s
        
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
      for taxid in self.root:
          self.set_of_descendants_aux(taxid)

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
                    set_of_survivors.add(t)
        set_of_survivors.update(set_of_taxid-lost_taxid)
        
        taxid_to_children_dict= {taxid: self.children[taxid] & set_of_survivors for taxid in set_of_survivors & self.children.keys()}
        taxid_to_name_dict = {taxid: self.name[taxid]  for taxid in set_of_survivors} 
        taxid_to_rank_dict = {taxid: self.rank[taxid]  for taxid in set_of_survivors}

        B=Taxonomy()
        B.children = taxid_to_children_dict
        B.name = taxid_to_name_dict
        B.rank=taxid_to_rank_dict
        B.init_root()
        B.init_descendants()
        B.init_parent()
        return B, lost_taxid

    # deprecated ?
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
        taxid_to_name_dict = {taxid: self.name[taxid]  for taxid in set_of_survivors} 
        taxid_to_rank_dict = {taxid: self.rank[taxid]  for taxid in set_of_survivors}

        B=Taxonomy()
        B.children = taxid_to_children_dict
        B.name = taxid_to_name_dict
        B.rank=taxid_to_rank_dict
        B.init_root()
        B.init_descendants()
        B.init_parent()
        return B
    
    def lca(self, set_of_taxid):
        """ lowest common ancestor """
        if not set_of_taxid.issubset(self):
            return None
        ancestor=next(iter(set_of_taxid))
        while not  ancestor in self.root and not set_of_taxid.issubset(self.descendants[ancestor])  :
            ancestor=self.parent[ancestor]
        if set_of_taxid.issubset(self.descendants[ancestor]):
            return ancestor
        else:
            return None

    def unary_ancestor(self, taxid):
        if taxid==None or taxid in  self.root:
            return None
        ancestor=taxid
        next_ancestor=self.parent[taxid]
        while not next_ancestor in self.root and self.number_of_children(next_ancestor)==1:
            ancestor=next_ancestor
            next_ancestor=self.parent[next_ancestor]
        if self.number_of_children(next_ancestor)==1:
            return next_ancestor
        else:
            return ancestor
        
    def hca(self, set_of_taxid):
        """ highest common ancestor: highest node whose descendant are set_of_taxid """
        hca=self.lca(set_of_taxid)
        return unary_ancestor(hca)
               
## end of class Taxonomy ##


def table_print_rec(t, taxid, rank):
    if taxid not in t.name:
        message.warning("Wrong taxid: "+ str(taxid))
        return
    print("."*rank + " ["+t.rank[taxid]+"] "+ taxid +" "+ t.name[taxid])
    if taxid not in t.children:
        # taxid is a leaf
        return
    for c in t.children[taxid]:
        table_print_rec(t, c, rank+1)
        
def table_print(t, taxid=None):
    if taxid:
        table_print(t,taxid,0)
    else:
        for tx in t.root: 
            table_print_rec(t,tx,0)
     

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
        for (i,line) in enumerate(in_file):
            columns = line.split("\t")
            if len(columns)<5 or len(columns[0])==0 :
                message.warning("File "+taxonomy_file+", line "+str(i+2)+": format error. Line is ignored")
            taxid_to_name_dict.update({utils.reduced(columns[0]):columns[2]})
            name_to_taxid_dict.update({columns[2]:utils.reduced(columns[0])})
            taxid_to_rank_dict.update({utils.reduced(columns[0]):columns[4].strip("\n")})
            taxid_to_parent_dict.update({utils.reduced(columns[0]):utils.reduced(columns[3])})

    excluded_taxid=taxid_to_parent_dict.values() - taxid_to_name_dict.keys()

    if len(excluded_taxid)>0:
        for taxid in excluded_taxid:
            message.warning("File "+taxonomy_file+": taxID "+str(taxid)+" is not documented. Ignored.")
        taxid_to_parent_dict={taxid:parent for taxid, parent  in  taxid_to_parent_dict.items() if parent not in excluded_taxid}
            
    for taxid in taxid_to_parent_dict.keys():
        utils.update_dictoset(taxid_to_children_dict,taxid_to_parent_dict[taxid],{taxid})

    taxonomy.children=taxid_to_children_dict
    taxonomy.name=taxid_to_name_dict
    taxonomy.rank=taxid_to_rank_dict
    taxonomy.parent=taxid_to_parent_dict
    taxonomy.init_root()
    taxonomy.init_descendants()
    
    return taxonomy

 
def build_flat_taxonomy(set_of_markers):
    """ constructs a flat taxonomy (the root is the set of species) from a set of markers """
    name={}
    rank={}
    for m in set_of_markers:
        name[m.taxid]= m.taxon_name
        rank[m.taxid]=""
    t=Taxonomy(name, rank)
    t.init_root()
    t.init_descendants()
    t.init_parent()
    return t

def search_taxid_from_taxon_name(taxon_name, taxonomy):
     for key, value in taxonomy.name.items():
         if taxon_name == value:
             return key

def create_taxonomy_file(taxonomy, outfile):
    file=open(outfile,"w")
    file.write("Taxon Id\tCommon name\tScientific name\tParent\tRank\n")
    for taxid in taxonomy:
        s=str(taxid)+"\t \t"+taxonomy.name[taxid]+"\t"
        if taxid not in taxonomy.root:
            s=s+str(taxonomy.parent.get(taxid))
        else:
            s=s+""
        s=s+" \t"+str(taxonomy.rank.get(taxid)+"\n")
        file.write(s)
    file.close()
