"""
taxonomy.py
"""

from src import markers, message, utils



class Taxonomy(object):
    def __init__(self, name={}, common_name={}, rank={}, children={}, descendants={}, parent={},  root=set()):
        self.name=name
        # key (str): taxid 
        # value (str): scientific name of the taxid
        self.common_name=common_name
        # key (str): taxid
        # value (str) : common name of the taxid
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
        s=self.name.keys()  # set of all nodes
        if len(s) == 0: # flat taxonomy
            self.root = self.name.keys()
            return
        for taxid in self.children.keys():
            s = s - self.children[taxid]
        self.root = s
        
    def set_of_descendants_aux(self, taxid):      
        if self.number_of_children(taxid)==0: 
            self.descendants[taxid]={taxid}
        else:
            s={taxid}
            for t in self.children[taxid]:
                self.set_of_descendants_aux(t)
                s.update(self.descendants[t])
            self.descendants[taxid]=s
            
    def init_descendants(self):
        self.descendants={}
        for taxid in self.root:
            self.set_of_descendants_aux(taxid)

    def intersection(self, set_of_taxid):
        """ create a sub-taxonomy that contains only taxid from set_of_taxid, with their ancestors"""
        """ elements of set_of_taxid not present in the taxonomy are lost. """
        """ B: new taxonomy"""
        """ lost_taxid: taxid from set_of_taxid that were not found in the taxonomy"""
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
        taxid_to_children= {}
        for taxid in set_of_survivors :
            if taxid in self.children and len(self.children[taxid] & set_of_survivors) > 0:
                taxid_to_children[taxid]= self.children[taxid] & set_of_survivors
        taxid_to_common_name = {taxid: self.common_name[taxid]  for taxid in set_of_survivors}
        taxid_to_name = {taxid: self.name[taxid]  for taxid in set_of_survivors}
        taxid_to_rank = {taxid: self.rank[taxid]  for taxid in set_of_survivors}
        B=Taxonomy()
        B.children = taxid_to_children
        B.name = taxid_to_name
        B.common_name=taxid_to_common_name
        B.rank=taxid_to_rank
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
        while not ancestor in self.root and not set_of_taxid.issubset(self.descendants[ancestor])  :
            ancestor=self.parent[ancestor]
        if set_of_taxid.issubset(self.descendants[ancestor]):
            return ancestor
        else:
            return None

    def unary_ancestor(self, taxid):
        if taxid is None or taxid in  self.root:
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
        return self.unary_ancestor(hca)
               
## end of class Taxonomy ##



def find_leaves(t, set_of_taxids):
    set_of_leaves={taxid for taxid in set_of_taxids if t.is_leaf(taxid)}
    return set_of_leaves

"""
t a été obtenue avec intersection_with_descendants auparavant
find the smallest subtree of t containing taxid + at least one more clade present in set_of_taxid
"""
def close_species(t, taxid, set_of_taxids):
    ancestor=t.unary_ancestor(taxid)
    ancestor=t.parent[ancestor] # TO DO: root test
    return (set_of_taxids.intersection(t.descendants[ancestor])).difference({taxid})

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
    if t is None:
        return
    if taxid:
        table_print_rec(t,taxid,0)
    else:
        for tx in t.root: 
            table_print_rec(t,tx,0)
     

def parse_taxonomy_simple_file(taxonomy_file):
    """
    taxonomy_file is a TSV file with 5 columns:
    Taxid | Common name	| Scientific name | Parent | Rank
    """
    if taxonomy_file is None:
        return None
    name_to_taxid_dict={} # key: (name, rank)
    taxid_to_children_dict={}
    taxid_to_name_dict={}
    taxid_to_common_name_dict={}
    taxid_to_rank_dict={}
    taxid_to_parent_dict={}
    with open(taxonomy_file) as in_file:
        next(in_file)
        for (i,line) in enumerate(in_file):
            columns = line.split("\t")
            if len(columns)<5 or len(columns[0])==0 :
                message.warning("File "+taxonomy_file+", line "+str(i+2)+": format error. Line is ignored")
            taxid_to_common_name_dict.update({utils.clean(columns[0]):columns[1].strip()})
            taxid_to_name_dict.update({utils.clean(columns[0]):columns[2].strip()})
            name_to_taxid_dict.update({columns[2].strip():utils.clean(columns[0])})
            taxid_to_rank_dict.update({utils.clean(columns[0]):columns[4].strip("\n")})
            taxid_to_parent_dict.update({utils.clean(columns[0]):utils.clean(columns[3])})
    excluded_taxid=taxid_to_parent_dict.values() - taxid_to_name_dict.keys()
    if len(excluded_taxid)>0:
        #for taxid in excluded_taxid:
        #    message.warning("File "+taxonomy_file+": taxID "+str(taxid)+" is not documented. Ignored.")
        taxid_to_parent_dict={taxid:parent for taxid, parent in taxid_to_parent_dict.items() if parent not in excluded_taxid}
    for taxid in taxid_to_parent_dict.keys():
        utils.update_dictoset(taxid_to_children_dict,taxid_to_parent_dict[taxid],{taxid})
    taxonomy=Taxonomy()
    taxonomy.children=taxid_to_children_dict
    taxonomy.name=taxid_to_name_dict
    taxonomy.common_name=taxid_to_common_name_dict
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
        name[m.taxid()]= m.taxon_name()
        rank[m.taxid()]=""
    t=Taxonomy(name, rank)
    t.init_root()
    t.init_descendants()
    t.init_parent()
    return t

def search_taxid_from_taxon_name(taxon_name, taxonomy):
    for key, value in taxonomy.name.items():
        if utils.equiv(taxon_name, value):
            return key
    return None

def search_taxid_from_common_name(common_name, taxonomy):
    for key, value in taxonomy.common_name.items():
        if utils.equiv(common_name, value):
            return key
    return None

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

def find_all_taxonomic_information_from_taxid(taxid, taxonomy):
    m=markers.Marker(field={"OX":taxid, "OS": taxonomy.name[taxid], "Rank": taxonomy.rank[taxid], "CommonName": taxonomy.common_name[taxid]})
    return m

# add taxid, taxon_name
def supplement_taxonomic_information(set_of_markers, taxo):
    taxa_with_missing_taxid={}
    for m in set_of_markers:
        if m.taxid() is None:
            if m.taxon_name() is None:
                None
            elif taxo is None:
                utils.update_dictoset(taxa_with_missing_taxid, m.taxon_name(), {m})
            else:
                m.field["OX"] = search_taxid_from_taxon_name(m.taxon_name(), taxo)
        else:
            if m.taxon_name() is None:
                if taxo is None:
                    None
                else:
                    m.field["OS"]=taxo.name[m.taxid()]
                    if taxo.common_name[m.taxid()] is not None and len(taxo.common_name[m.taxid()]) > 0:
                        m.field["Common Name"] = taxo.common_name[m.taxid()]
                    m.field["Rank"] = taxo.rank[m.taxid()]
            else:
                None
    for i,taxon in enumerate(list(taxa_with_missing_taxid.keys())):
        s={m.taxid() for m in set_of_markers if m.taxon_name()==taxon and m.taxid() is not None}
        if len(s)==0:
            new_taxid=''.join(word[0] for word in taxon.split())+str(i+1)
        else:
            new_taxid=s.pop()
        for m in taxa_with_missing_taxid[taxon]:
            m.field["OX"]= new_taxid
    return set_of_markers
    
def add_taxonomy_ranks(set_of_markers, t):
    if t is None:
        message.warning("No taxonomy provided. Unable to apply TAXONOMY completion.")
        return set_of_markers
    for m in set_of_markers:
        if m.taxid() not in t.common_name or m.taxid() not in t.rank or m.taxid() not in t.parent:
            message.warning("TaxID "+ str(m.taxid()) + " not found in TAXONOMY file.")
            FINISHED=True
        else:
            m.field["Common Name"]=t.common_name[m.taxid()]
            m.field["Rank"]=t.rank[m.taxid()]
            parent=t.parent[m.taxid()]
            FINISHED=False
        while not FINISHED:
            m.field[t.rank[parent]] = t.name[parent]
            if parent in t.parent:
                parent=t.parent[parent]
            else:
                FINISHED=True
    return set_of_markers

def find_closest_ID(target, set_of_taxids, taxo):
    node=target
    set_of_taxids={t for t in set_of_taxids if t is not None}
    if len(set_of_taxids)==0:
        return "."
    if target not in taxo:
        return ". "
    while node not in set_of_taxids and len(set_of_taxids & taxo.descendants[node])==0:
        if node in taxo.parent:
            node=taxo.parent[node]
        else:
            return ""
    return "["+taxo.rank[node]+"]. "

def update_taxonomy_and_set_of_markers(taxonomy, set_of_markers, set_of_mandatory_taxid=set()):
    if taxonomy is None:
        return build_flat_taxonomy(set_of_markers), set_of_markers
    set_of_taxid={m.taxid() for m in set_of_markers} | set_of_mandatory_taxid
    secondary_taxonomy, lost_taxid = taxonomy.intersection(set_of_taxid)
    set_of_markers = markers.remove_lost_taxid(set_of_markers, lost_taxid)
    for taxid in lost_taxid:
        message.warning("TaxID " + str(
            taxid) + " not found in taxonomy file. All markers associated to this TaxID are ignored.")
    return secondary_taxonomy, set_of_markers

def merge_taxonomy(set_of_markers, taxonomy):
    if taxonomy is None:
        return build_flat_taxonomy(set_of_markers)
    set_of_survivors = set()
    missing_taxid = set()
    set_of_taxid={m.taxid() for m in set_of_markers}
    for taxid in set_of_taxid:
        if taxid not in taxonomy.name:
            missing_taxid.add(taxid)
        else:
            t = taxid
            while t in taxonomy.parent:
                t = taxonomy.parent[t]
                set_of_survivors.add(t)
    set_of_survivors.update(set_of_taxid - missing_taxid)
    taxid_to_children = {}
    for taxid in set_of_survivors:
        if taxid in taxonomy.children and len(taxonomy.children[taxid] & set_of_survivors) > 0:
            taxid_to_children[taxid] = taxonomy.children[taxid] & set_of_survivors
    taxid_to_common_name = {taxid: taxonomy.common_name[taxid] for taxid in set_of_survivors}
    taxid_to_name = {taxid: taxonomy.name[taxid] for taxid in set_of_survivors}
    taxid_to_rank = {taxid: taxonomy.rank[taxid] for taxid in set_of_survivors}
    for m in set_of_markers:
        if m.taxid() in missing_taxid:
            taxid_to_name[m.taxid()] = m.taxon_name()
            taxid_to_rank[m.taxid()] = m.rank()
    if len(missing_taxid)>0:
        message.warning("Some TaxID were not found in the taxonomy file. \n" + str(missing_taxid)+"\nTaxonomic information ignored for associated markers.")
    B = Taxonomy()
    B.children = taxid_to_children
    B.name = taxid_to_name
    B.common_name = taxid_to_common_name
    B.rank = taxid_to_rank
    B.init_root()
    B.init_descendants()
    B.init_parent()
    print("R O O T:", str(B.root))
    return B



def find_taxid(target, taxonomy):
    if target in taxonomy:
        return target
    taxid = search_taxid_from_taxon_name(target, taxonomy)
    if taxid is not None:
        return taxid
    taxid = search_taxid_from_common_name(target, taxonomy)
    if taxid is not None:
        return taxid
    message.escape("The clade "+target+ " is not found in TAXONOMY file. Stopping execution.")