import json

def build_assignments_json(output, list_of_assignments, taxonomy, B):
# list of dictionnaries
list_of_lines=[]
     for i in range(len(list_of_assignments)):
        a=list_of_assignments[i]
        dict={}
        if len(a.results)==0: # no assignment
            dict[name]=a.spectrum_name
            list_of_lines.append(dict)
            continue
        for j in range(len(a.results)): 
            name_to_peak={(name,ptm):peak for (peak,name,ptm) in a.combinations_of_peaks[j]}
            s=a.spectrum_name+"\t"
            for (name,ptm) in list_of_marker_full_names:
                if (name,ptm) in  name_to_peak:
                    s=s+str(round(name_to_peak[(name,ptm)],3))
                s=s+"\t"
            if taxonomy:
                lca=a.results[j][0]
                if lca!=None:
                    hca=B.unary_ancestor(lca)
                    s=s+"{:.2e}".format(a.results[j][1])+"\t"+str(a.results[j][2])+"\t"+str(lca) + " ["+ B.name[lca] + "]\t" + B.rank[lca] +"\t"+str(hca)+ " ["+ B.name[hca] + "]\t" + B.rank[hca] +"\t"
                else:
                    s=s+"\t None\t\t\t\t" 
            else:
                s=s+"{:.2e}".format(a.results[j][1])+"\t"+str(a.results[j][2])+"\t"
            for taxid in a.results[j][3]: 
                s=s+str(taxid)+" ["+B.name[taxid]+"] "
            s=s+"\n"


       self.spectrum_name = spectrum_name
        self.peak_to_taxid_name_ptm=peak_to_taxid_name_ptm
        self.combinations_of_peaks=combinations_of_peaks #list of combinations of peaks. Peak=(peak,name,ptm)
        self.results=results # list of 3-uplets (lca, score, set of 2-uplets (pvalue, taxid)), all taxid  have the same markers 
        self.taxid_to_peak_name_ptm=taxid_to_peak_name_ptm # index in combinations_of_peaks
           

# Define the table with headings as a list of dictionaries
table = [
    {"Name": "John", "Age": 30, "City": "New York"},
    {"Name": "Alice", "Age": 25, "City": "Los Angeles"},
    {"Name": "Bob", "Age": 35, "City": "Chicago"}
]

# Specify the path to the JSON file
file_path = "table_with_headings.json"

# Write the table to the JSON file
with open(file_path, "w") as json_file:
    json.dump(table, json_file, indent=4)
