#!/usr/bin/env python3

import argparse
import sys
import time
import os
import logging

# local import
from src import mass_spectrum 
from src import assignment 
from src import peptide_table as pt
from src import sequences as seq
from src import taxonomy as ta 
from src import markers 
from src import fasta_parsing as fa
from src import compute_masses
from src import message
from src import limit as lmt
from src import config
from src import params_checker

def create_report(report_file, output_dir, spectra, final_list_of_spectra, taxonomy, taxonomy_tree, peptide_table, fasta, directory,  limit, deamidation, error, neighbour, all, output, output2, new_table, set_of_markers, WEB):

    # TO DO: display constraints

    report_name=os.path.join(output_dir, report_file)
    sys.stdout=open(report_name, 'w')
    print("============================================\n")
    print("              PAMPA CLASSIFY\n")
    print("============================================\n")
    print (time.ctime())
    if WEB:
        if peptide_table:
            peptide_table_file = ""
            for pt in peptide_table:
                peptide_table_dir, this_peptide_table_file = os.path.split(pt)
                peptide_table_file =  peptide_table_file + " "+ this_peptide_table_file
        if fasta:
            fasta_dir, fasta_file = os.path.split(fasta)
        if taxonomy:
          taxo_dir, taxo_file = os.path.split(taxonomy)
        else:
          taxo_file = None
        spectra_files = str([spectrum.name for spectrum in final_list_of_spectra])
    else:
        peptide_table_file=peptide_table
        fasta_file=fasta
        taxo_file=taxonomy
        spectra_files = str([spectrum.name for spectrum in final_list_of_spectra])
        spectra_files= spectra+"/ "+ spectra_files 

          
    print("")
    print("---------------------------------")
    print("   INPUT and PARAMETERS")
    print("---------------------------------\n")
    if peptide_table:
        print("  Peptide table (markers) : "+ str(peptide_table_file))
    elif fasta:
        print("  Fasta file (markers)    : "+str(fasta_file))
    else:
        print("  Fasta file (markers)    : "+str(directory))
    print("  Taxonomy                : "+str(taxo_file))
    print("  Mass spectra            : "+ str(spectra_files))
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
    if os.path.getsize(os.path.join(output_dir,'warning.log')) > 0:
        print("---------------------------------")
        print("   WARNING")
        print("---------------------------------\n")
        print("  Warnings were raised regarding your inputs.\n")
        with open(os.path.join(output_dir,'warning.log'), 'r') as file:
            for line in file:
                print("  - "+line, end="")
        print("")
    print("---------------------------------")
    print("   OUTPUT FILES")
    print("---------------------------------\n")
    print("  Main result file       : "+output+".tsv")
    print("  More details           : "+output2)
    if not peptide_table:
        print("  New peptide table      : " + new_table)
    print("  Report (this file)     : " + report_name)
    print("")
    print("---------------------------------")
    print("  MASS SPECTRA")
    print("---------------------------------\n")
    print("  "+str(len(final_list_of_spectra))+" files found\n")
    for f in final_list_of_spectra:
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
    markers.colinearity(set_of_markers)
    markers.check_set_of_markers(set_of_markers)
    if taxonomy is not None:
        print("---------------------------------")
        print("  TAXONOMY")
        print("---------------------------------")
        ta.table_print(taxonomy_tree)
    
    sys.stdout = sys.__stdout__
            

def main(spectra, taxonomy, peptide_table, fasta, directory, limit, deamidation, error, neighbour, allsolutions, output, mammals, web):
   
    try:

        (spectra, taxonomy, peptide_table, fasta, directory, limit, deamidation, error, neighbour, allsolutions, output, output2, output_dir, report_file, new_table)=params_checker.check_and_update_parameters_classify(spectra, taxonomy, peptide_table, fasta, directory, limit, deamidation, error, neighbour, allsolutions, output, mammals)


        # parsing spectra files
        final_list_of_spectra = []
        spectra_dir = os.listdir(spectra)
        for f in spectra_dir:
            file_name = os.path.join(spectra, f)
            list_of_spectra=mass_spectrum.parser(file_name,f)
            if len(list_of_spectra)>0:
                final_list_of_spectra.extend(list_of_spectra)
        if len(final_list_of_spectra)==0:
                message.escape("No valid spectra found.\n Please refer to the warning.log file for more detail.")
        final_list_of_spectra.sort(key=lambda x: x.name)
        
            
        # parsing taxonomy    
        if taxonomy is not None:   
            primary_taxonomy=ta.parse_taxonomy_simple_file(taxonomy)
        else:
            primary_taxonomy=None
            
        # limit file
        if limit:
            list_of_constraints=lmt.parse_limits(limit)
        else:
            list_of_constraints=[]
    
        # parsing models for organisms and applying limits    
        if peptide_table :
            set_of_markers, _ = pt.parse_peptide_tables(peptide_table, list_of_constraints, primary_taxonomy)
            set_of_markers=ta.supplement_taxonomic_information(set_of_markers, primary_taxonomy) # To check: see supplement
            set_of_markers=compute_masses.add_PTM_or_masses_to_markers(set_of_markers)
        if fasta or directory:
            set_of_sequences = fa.build_set_of_sequences(fasta, directory, list_of_constraints, primary_taxonomy)
            config_digestion=config.config_digestion()
            if len(set_of_sequences)==0:
                 message.escape("Fasta file(s): No valid sequences found.\nPlease refer to the warning.log file to trace back the errors.")
            set_of_markers=compute_masses.add_PTM_or_masses_to_markers(seq.in_silico_digestion(set_of_sequences,config_digestion))
        if len(set_of_markers)==0:
            message.escape("No valid peptide marker found.\nPlease refer to the warning.log file to trace back the errors.")

        if deamidation:
            set_of_codes=set()
            for constraint in list_of_constraints:
                if 'Deamidation' in constraint:
                    set_of_codes.update(constraint['Deamidation'])
            set_of_markers=set_of_markers.union(compute_masses.add_deamidation(set_of_markers, set_of_codes))
       
        set_of_markers=markers.sort_and_merge(set_of_markers)
       
        # construction of the secondary taxonomy and suppression of taxid not present in the taxonomy
        set_of_taxid={m.taxid() for m in set_of_markers}
        if primary_taxonomy :  
            secondary_taxonomy, lost_taxid=primary_taxonomy.intersection(set_of_taxid)
            set_of_markers= markers.remove_lost_taxid(set_of_markers, lost_taxid)
            for taxid in lost_taxid:
                message.warning("File "+taxonomy+": TaxID "+str(taxid)+" not found. All markers associated to this TaxID are ignored.") 
        else:
            secondary_taxonomy=ta.build_flat_taxonomy(set_of_markers)
        
        if fasta or directory:
            pt.build_peptide_table_from_set_of_markers(set_of_markers, new_table)
        
        create_report(report_file, output_dir, spectra, final_list_of_spectra, taxonomy, secondary_taxonomy, peptide_table, fasta,  directory, limit, deamidation, error, neighbour, all, output, output2, new_table, set_of_markers,web)
        
        # species identification
       
        assignment.assign_all_spectra(final_list_of_spectra, set_of_markers, error, taxonomy, secondary_taxonomy, neighbour, allsolutions, output, output2)

        if not web:
            print("")
            print("   Job completed.")
            print("   All results are available in the following files.") 
            print("")
            print(f"   - Assignments       : {output}.tsv")
            print(f"   - More detail       : {output2}")
            print(f"   - Report on the run : {os.path.join(output_dir, report_file)}")
            print("")
    # TO DO: add the new peptide table, if necessary

            if os.path.getsize(os.path.join(output_dir, "warning.log")) > 0:
                print("Warnings were raised regarding your inputs.")
                print("Please refer to the warning.log file for detail.")

    except message.InputError:
        if not web:
           print("\n   An error occured with your input. Stopping execution.")
           print("   Please refer to the warning.log files for more detail.")
        else:
           pass
 
if __name__ == "__main__":
    main()
