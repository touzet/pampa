#!/usr/bin/env python3

"""

   main_assign.py                              

"""


import argparse
import sys
import time
import os
from os import listdir
from os.path import join
import mass_spectrum 
from assignment import *
from peptide_table import *
from sequences import *


from taxonomy import*
from markers import *
from fasta_parsing import *
from compute_masses import *

def escape(message):
    print("\n "+message)
    sys.exit(" Stopping execution.\n\n Type <pampa assign -h> for more help.\n")


def check_and_update_parameters(spectra, taxonomy, peptide_table, fasta, directory, limit, error, neighbour, all, output, mammals):
    """
    Parameters checking and fixing
    """

    if spectra is None:
        escape("Missing parameter: spectra (-s). ")
    if not os.path.isdir(spectra):    
        escape("Directory "+spectra+" not found.")
    if error is None:
        escape("Missing parameter: error (-e)")
    if error<0:
        escape("Parameter error (-e) should be positive.")
    if output is None:
        escape("Missing parameter: output (-o)")

    if mammals :
        if taxonomy:
            print("Ignored option: -t "+taxonomy+"\n")
        if peptide_table:
            print("Ignored option: -p "+peptide_table+"\n")
        if fasta:
            print("Ignored option: -f "+fasta+"\n")
        if  limit:
          print("Ignored option: -l "+limit+"\n")
        taxonomy="../Taxonomy/taxonomy_mammals.tsv"
        peptide_table="../Peptide_tables/table_mammals_with_deamidation.tsv"

    if taxonomy:
        if not os.path.isfile(taxonomy):
            escape("File "+taxonomy+" not found.")

    if limit:
        if not os.path.isfile(limit):
            escape("File "+limit+" not found")
            
    if neighbour not in range(101):
        neighbour=100

    if all and neighbour is None:
        print ("Ignored parameter: -a (all). This parameter comes with -n (near-optimal solutions).")
        
    q = (peptide_table, fasta, directory)
    if not (q[0] or q[1] or q[2]):
        escape("Missing information for marker peptides (-p, -f or -d).")
    if (q[0] and q[1]) or (q[0] and q[2]) or (q[1] and q[2]) :
        escape("Options -p (peptide_table), -f (fasta) and -d (directory of fasta files) are mutually incompatible.")
    
    if q[0]:
        if not os.path.isfile(peptide_table):
           sys.exit("File "+peptide_table+" not found.")


    output_dir, output_file = os.path.split(output)
    # Ensure the output directory exists
    if len(output_dir)>0 and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    extension=output_file[-4:].lower()
    if extension!=".tsv":
        output_file=output_file+".tsv"
    else:
        output_file=output_file[:-4]+".tsv"

    output=os.path.join(output_dir, output_file)
    output2=os.path.join(output_dir, "detail_"+output_file)
    report_file="report_"+output_file.replace("tsv", "txt")
    report_path=os.path.join(output_dir, report_file)
        
    return (spectra, taxonomy, peptide_table, fasta, directory, limit, error, neighbour, all, output, output2, report_path, mammals)

        
def create_report(report_name, spectra, list_of_spectra, taxonomy, taxonomy_tree, lost_taxid, peptide_table, fasta, limit, directory, error, neighbour, output, output2,  set_of_markers):
    
    sys.stdout=open(report_name, 'w')
    print (time.ctime())
    print("")
    print("---------------------------------")
    print("   PARAMETERS")
    print("---------------------------------\n")
    if peptide_table:
        print("  Peptide table (markers) : "+ str(peptide_table))
    elif fasta:
        print("  Fasta file (markers)    : "+str(fasta))
    else:
        print("  Fasta file (markers)    : "+str(directory))
    print("  Taxonomy                : "+str(taxonomy))
    print("  Mass spectra            : "+spectra)
    print("  Near-optimal solutions  : ",end="")
    if neighbour==100:
        print("only optimal solutions")
    else:
        print("up to "+str(neighbour)+"% suboptimality")
    print("  Error margin            : "+str(error)+"Da/ppm")
    print("")
    print("---------------------------------")
    print("   OUTPUT FILES")
    print("---------------------------------\n")
    print("  Main result file         : "+output)
    print("  More details             : "+output2)
    if not peptide_table:
        print("  Table of peptides        : table_"+output) #/!\
    print("  Report (this file)       : " + report_name)
    print("")
    print("---------------------------------")
    print("  MASS SPECTRA")
    print("---------------------------------\n")
    print("  "+str(len(list_of_spectra))+" files found\n")
    for f in list_of_spectra:
        print("  "+f.name+ " ("+str(len(f))+" peaks)")
    print("")
    if peptide_table:
        print("---------------------------------")
        print("  PEPTIDE TABLE")
        print("---------------------------------\n")
    if fasta or directory:
       print("---------------------------------")
       print("  FASTA SEQUENCES")
       print("---------------------------------\n")
    colinearity(set_of_markers)
    check_set_of_markers(set_of_markers)
    if taxonomy is not None:
        print("---------------------------------")
        print("  TAXONOMY")
        print("---------------------------------")
        table_print(taxonomy_tree)
        if len(lost_taxid)>0:
            print ("\nTaxid not found : "+str(lost_taxid))
            print ("The associated markers are removed from the remaining of the analysis.")
    sys.stdout = sys.__stdout__
    
def main(spectra=None, taxonomy=None, peptide_table=None, fasta=None, directory=None, limit=None, error=None, neighbour=None, all=None, output=None, mammals=None):

    
    (spectra, taxonomy, peptide_table, fasta, directory, limit, error, neighbour, allsolutions, output, output2, report_path, mammals)=check_and_update_parameters(spectra, taxonomy, peptide_table, fasta, directory, limit, error, neighbour, all, output,mammals)

    # parsing spectra files
    list_of_spectra=[]
    spectra_dir=os.listdir(spectra)
    for f in spectra_dir:
        file_name= os.path.join(spectra, f)
        spectrum=mass_spectrum.parser(file_name,f)
        list_of_spectra.append(spectrum)
        list_of_spectra.sort(key=lambda x: x.name)

    # parsing models for organisms    
    if peptide_table:
        set_of_markers = parse_peptide_table(peptide_table)
        if limit:
            set_of_markers = extract_relevant_markers(set_of_markers, limit)
    if fasta :
        set_of_sequences= build_set_of_sequences_from_fasta_file(fasta,True)
        if limit:
            set_of_sequences=extract_relevant_sequences(set_of_sequences, limit)
        set_of_markers=add_PTM_or_masses_to_markers(in_silico_digestion(set_of_sequences,1,12,33, False))
        #build_peptide_table_from_set_of_markers(set_of_markers,"t able_"+output,"")
    if directory:
        set_of_sequences=build_set_of_sequences_from_fasta_dir(directory)
        if limit:
            set_of_sequences=extract_relevant_sequences(set_of_sequences, limit)
        set_of_markers=add_PTM_or_masses_to_markers(in_silico_digestion(set_of_sequences,1,12,33, False))
        #build_peptide_table_from_set_of_markers(set_of_markers,"table_"+output,"")

    
    # extraction of taxid, masses and names    
    mass_taxid_name_list, set_of_taxid, taxid_to_taxidname = build_marker_dictionnaries_from_set_of_markers(set_of_markers)
 
    # construction of the taxonomy  
    if taxonomy is not None:   
        primary_taxonomy=parse_taxonomy_simple_file(taxonomy)
        secondary_taxonomy, lost_taxid=primary_taxonomy.intersection(set_of_taxid)
        set_of_markers= remove_lost_taxid(set_of_markers, lost_taxid)
    else:
        secondary_taxonomy=build_flat_taxonomy(set_of_markers)

    create_report(report_path, spectra, list_of_spectra, taxonomy, secondary_taxonomy, lost_taxid, peptide_table, fasta,  directory, limit, error, neighbour, output, output2, set_of_markers)
        
    # species identification
    assign_all_spectra(list_of_spectra, mass_taxid_name_list, error, taxonomy, secondary_taxonomy, neighbour,allsolutions, output,output2)

    print("")
    print("   Job completed.")
    print("   All results are available in the three following files.") 
    print("")
    print(f"   - Assignments       : {output}")
    print(f"   - More details      : {output2}")
    print(f"   - Report on the run : {report_path}")
    print("")


