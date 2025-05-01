import os
import time
import sys

# local import
from src import markers


def print_set_of_sequences(set_of_sequences):
    for seq in set_of_sequences:
        print ("     "+seq.seqid()+"\t"+seq.protein()+ "\t "+seq.taxon_name()+" (TaxID:"+seq.taxid()+")")
    print("")
    
def print_peptide_table(peptide_table, web):
    if web:
        _, peptide_table_file = os.path.split(peptide_table)
        print("  Input peptide table: " + peptide_table_file)
    else:
        print("  Input peptide table: " + peptide_table)
    
def print_peptide_tables(peptide_table, web):
    if peptide_table is None:
        return
    if len(peptide_table)==1:
        print("---------------------------------")
        print("     INPUT PEPTIDE TABLE         ")
        print("---------------------------------")
        print("")
    else:
        print("---------------------------------")
        print("     INPUT PEPTIDE TABLES        ")
        print("---------------------------------")
        print("")
    if web:
        for pep in peptide_table:
            _, pep_file = os.path.split(pep)
            print("  "+pep_file, " ")
    else:
        for pep in peptide_table:
            print("  "+pep, " ")

    
def print_taxonomy_file(taxonomy, web):
    if taxonomy:
        if web:
            _, taxo_file = os.path.split(taxonomy)
            print("  Taxonomy: "+taxo_file)
        else:
            print("  Taxonomy: "+taxonomy)
            
            
def print_spectra_dir(spectra_dir, spectra_files, web):
    if spectra_dir:
        if web:
            _, taxo_file = os.path.split(taxonomy)
            print("  Taxonomy: "+taxo_file)
        else:
            print("  Taxonomy: "+taxonomy)
            

def print_spectra(spectra_dir, list_of_spectra, web):
    if not web:
        print("  Directory: "+spectra_dir)
    print("  "+str(len(list_of_spectra))+" files found\n")
    for f in list_of_spectra:
        print("  " + f.name + " (" + str(len(f)) + " peaks)")
    print("")
            
def print_digestion(config_digestion):
    print("     In silico digestion:")
    print("     - Enzyme: "+ config_digestion["enzyme"])
    print("     - Maximal number of missed cleavages: "+str(config_digestion["number_of_missed_cleavages"]))
    print("     - Minimal peptide length: "+str(config_digestion["min_peptide_length"]))
    print("     - Maximal peptide length: "+str(config_digestion["max_peptide_length"]))
    
def create_report_classify(spectra_dir, list_of_spectra, taxonomy, taxonomy_tree, peptide_table, fasta, fasta_dir, set_of_sequences, set_of_markers, limit, deamidation, error, neighbour, all, new_table, config_digestion, config_nb_of_peaks, web):

    # TO DO: display constraints

    print ("PAMPA CLASSIFY\n\n")
    print("---------------------------------")
    print("  MASS SPECTRA")
    print("---------------------------------\n")
    print_spectra(spectra_dir, list_of_spectra, web)
    print("")
    print("---------------------------------")
    print("  PEPTIDE MARKERS")
    print("---------------------------------\n")
    if peptide_table:
        print_peptide_tables(peptide_table, web)
        markers.short_colinearity(set_of_markers)
    else:
        print_fasta_file(fasta, web)
        print_fasta_dir(fasta, web)
        print_set_of_sequences(set_of_sequence)
        print_digestion(config_digestion)
        print("\n  New peptide table  : " + new_table)
   
    #print_limit(limit_file, web)
    print("")
  
    print("---------------------------------")
    print("   PARAMETERS")
    print("---------------------------------\n")
    print("  Minimum number of peaks :" + str(config_nb_of_peaks))
    print("  Near-optimal solutions  : ",end="")
    if neighbour==100:
        print("only solution with the highest number of matching peaks")
    else:
        print("up to "+str(neighbour)+"% matching peaks")
    print("  Selection of solutions  : ", end="")
    if not all:
        print ("peak intensity and inclusion selection")
    else:
        print ("all possible solutions")
    print("  Error margin tolerance  : "+str(error), end=" ")
    if error<1:
        print("Da")
    else:
        print("ppm")
    print("")
    if deamidation:
        print ("  Deamidation            : Yes")
    else:
        print ("  Deamidation            : None")
    #markers.colinearity(set_of_markers)
    #markers.check_set_of_markers(set_of_markers)
    if taxonomy :
        print("---------------------------------")
        print("  TAXONOMY")
        print("---------------------------------")
        print_taxonomy(taxonomy, web)
        ta.table_print(taxonomy_tree)

def create_report_homology(peptide_table, set_of_markers, set_of_sequences, taxonomy, web):
    print("PAMPA CRAFT, mode HOMOLOGY")
    print("")
    print_peptide_tables(peptide_table, web)
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
    
def create_report_deamidation(peptide_table, set_of_codes, web):
    print("PAMPA CRAFT, mode DEAMIDATION\n")
    if len(set_of_codes)==0:
        print("  Modified peptide markers: all\n")
    else:
        print("  Modified peptide markers: "+ str(set_of_codes)+"\n")
    print_peptide_tables(peptide_table, web)
    
   
def create_report_allpeptides(set_of_sequences, config_digestion):
    print("PAMPA CRAFT, mode ALL PEPTIDES")
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
    print("PAMPA CRAFT, mode SUPPLEMENT \n")
    print("---------------------------------")
    print("   INPUT PEPTIDE TABLE           ")
    print("---------------------------------")
    print_peptide_tables(peptide_table, web)
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

def create_report_header(command_line, report):
    sys.stdout=open(report, 'w')
    print("============================================\n")
    print("                    PAMPA  \n")
    print("============================================\n")
    print (time.ctime())
    print("")
    print(command_line)
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
    
