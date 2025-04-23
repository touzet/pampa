import os
import time
import sys

# local import
from src import markers


def print_set_of_sequences(set_of_sequences):
    for seq in set_of_sequences:
        print ("     "+seq.seqid()+"\t"+seq.protein()+ "\t "+seq.taxon_name()+" (TaxID:"+seq.taxid()+")")
    
def print_peptide_table(peptide_table, web):
    if web:
        peptide_table_dir, peptide_table_file = os.path.split(peptide_table)
        print("  Input peptide table: " + peptide_table_file)
    else:
        print("  Input peptide table: " + peptide_table)
    
def print_taxonomy_file(taxonomy, web):
    if taxonomy:
        if web:
            taxo_dir, taxo_file = os.path.split(taxonomy)
            print("  Taxonomy: "+taxo_file)
        else:
            print("  Taxonomy: "+taxonomy)
            
def print_digestion(config_digestion):
    print("     In silico digestion:")
    print("     - Enzyme: "+ config_digestion["enzyme"])
    print("     - Maximal number of missed cleavages: "+str(config_digestion["number_of_missed_cleavages"]))
    print("     - Minimal peptide length: "+str(config_digestion["min_peptide_length"]))
    print("     - Maximal peptide length: "+str(config_digestion["max_peptide_length"]))
            
def create_report_homology(peptide_table, set_of_markers, set_of_sequences, taxonomy, web):
    print("MODE : HOMOLOGY")
    print_peptide_table(peptide_table, web)
    print_taxonomy_file(taxonomy, web)
    print("---------------------------------")
    print("   INPUT SEQUENCES")
    print("---------------------------------")
    print_set_of_sequences(set_of_sequences)
    print("---------------------------------")
    print("   NEW PEPTIDE TABLE")
    print("---------------------------------\n")
    markers.colinearity(set_of_markers)
    markers.check_set_of_markers(set_of_markers)
    
def create_report_deamidation(peptide_table, set_of_codes):
    print("MODE : DEAMIDATION")
    print_peptide_table(peptide_table, web)
    if len(set_of_codes)==0:
        print("  Modified peptide markers: all")
    else:
        print("  Modified peptide markers: "+ str(set_of_codes))
   
   
def create_report_allpeptides(set_of_sequences, config_digestion):
    print("MODE : ALL PEPTIDES")
    print("")
    print("---------------------------------")
    print("   PARAMETERS")
    print("---------------------------------")
    print("")
    print_digestion(config_digestion)
    print("")
    print("---------------------------------")
    print("   INPUT SEQUENCES")
    print("---------------------------------")
    print("")
    print_set_of_sequences(set_of_sequences)
        
def create_report_supplement(peptide_table,list_of_markers=None, set_of_sequences=None):
    print("MODE : SUPPLEMENT\n")
    print("---------------------------------")
    print("   INPUT PEPTIDE TABLE           ")
    print("---------------------------------")
    print_peptide_table(peptide_table, web)
    if taxonomy:
        print("---------------------------------")
        print("   TAXONOMY                      ")
        print("---------------------------------")
        print_taxonomy_file(taxonomy, web)
    if set_of_sequences:
        print("---------------------------------")
        print("   INPUT SEQUENCES               ")
        print("---------------------------------")
        print_set_of_sequences(set_of_sequences)
        print("---------------------------------")
        print(" NEW PEPTIDE TABLE               ")
        print("---------------------------------")
        markers.colinearity(set(list_of_markers))
        markers.check_set_of_markers(set(list_of_markers))

def create_report_header(report):
    sys.stdout=open(report, 'w')
    print("============================================\n")
    print("                PAMPA CRAFT\n")
    print("============================================\n")
    print (time.ctime())
    print("")
    

def create_report_footer(output_dir, output, report):
    if os.path.getsize(os.path.join(output_dir,'warning.log')) > 0 and os.path.getsize(os.path.join(output_dir,'error.log'))==0:
        print("")
        print("---------------------------------")
        print("   WARNINGS")
        print("---------------------------------\n")
        with open(os.path.join(output_dir,'warning.log'), 'r') as file:
            for line in file:
                print("  "+line, end="")
    if os.path.getsize(os.path.join(output_dir,'error.log')) > 0:
        print("")
        with open(os.path.join(output_dir,'warning.log'), 'r') as file:
            for line in file:
                print("  "+line, end="")
        print("")
        print("---------------------------------")
        print("   * * *  NO OUTPUT FILE  * * *  ")
        print("---------------------------------")
    else:
        print("")
        print("---------------------------------")
        print("   OUTPUT FILES")
        print("---------------------------------\n")
        print("  Main result file (TSV) : "+output)
        print("  Report (this file)     : "+report)
    print("")
    sys.stdout = sys.__stdout__
    
