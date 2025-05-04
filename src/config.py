#!/usr/bin/env python3

import json

def parse_config_file(file_name="config.json"):
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
    return data["marker_order"]

def config_minimum_number_of_peaks(file_name):
    return int(parse_config_file(file_name)["min_number_of_peaks"])

def config_selection_peaks(file_name):
    return float(parse_config_file(file_name)["min_proportion_of_peaks"])
