#!/usr/bin/env python3

import json

def parse_config_file(file_name):
    with open(file_name) as json_file:
        data = json.load(json_file)
    return data
# TO DO: add PTM definition here

def config_digestion(file_name):
    data=parse_config_file(file_name)
    return {k: data[k] for k in {"enzyme", "number_of_missed_cleavages",
        "min_peptide_length", "max_peptide_length"}}
        
def config_headers(file_name):
    data=parse_config_file(file_name)
    return data["taxonomy_ranks"]+data["peptide_table_order"]
    
def config_markers(file_name):
    data=parse_config_file(file_name)
    return data.get("marker_order")

def config_minimum_number_of_peaks(file_name):
    data=parse_config_file(file_name)
    nb=data.get("min_number_of_peaks")
    if nb is None:
        return 0
    else:
        return int(nb)

def config_selection_peaks(file_name):
    data=parse_config_file(file_name)
    nb=data.get("min_proportion_of_peaks")
    if nb is None:
        return 0.0
    else:
        return float(nb)

def config_taxonomy(file_name):
    return parse_config_file(file_name)["taxonomy"]

def config_peptide_table(file_name):
    return parse_config_file(file_name).get("peptide_table")

def config_matrices_and_co(file_name):
    return parse_config_file(file_name)["substitution_matrices"], parse_config_file(file_name)["gamma_matrices"], parse_config_file(file_name)["conserved"]