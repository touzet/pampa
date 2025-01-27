#!/usr/bin/env python3

import json

def parse_config_file():
    with open('config.json') as json_file:
        data = json.load(json_file)
        return data
# TO DO: add PTM definition here
# TO DO: add warning when the file is not found

def sort_headers(set_of_headers):
    with open('config.json') as json_file:
        data = json.load(json_file)
    list_of_selected_headers=[]
    list_of_other_headers=[]
    list_of_config_headers=[]
    for key in data:
        if  isinstance(data[key], list):
            list_of_config_headers.extend(data[key])
    for element in set_of_headers:
        if element in list_of_config_headers:
            list_of_selected_headers.append((list_of_config_headers.index(element), element))
        else:
            list_of_other_headers.append(element)
    list_of_selected_headers.sort(key=lambda x:x[0])
    return [z[1] for z in list_of_selected_headers]+list_of_other_headers
