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




def check_and_update_parameters(spectra, taxonomy, peptide_table, fasta, directory, limit, error, neighbour, all, output, mammals, light):
    """
    Parameters checking and fixing. Configuration of loggers
    """

    if output is None:
        message.configure("")
        message.escape("Missing parameter: output (-o).")
        
    output_dir, output_file = os.path.split(output)
    
    if len(output_dir)>0 :
        # Ensure the output directory exists. If not, create it.
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    message.configure(output_dir)
            
    if spectra is None:
        message.escape("Missing parameter: spectra (-s).")
    if not os.path.isdir(spectra):    
        message.escape("Directory "+spectra+" not found.")
    if error is None:
        message.escape("Missing parameter: error (-e)")
    if error<0:
        message.escape("Parameter error (-e) should be positive.")
    if output is None:
        message.escape("Missing parameter: output (-o)")

    if mammals :
        if not os.path.isfile("Taxonomy/taxonomy_mammals.tsv"):
            message.escape("The taxonomy file has been deleted.")
        if not os.path.isfile("Peptide_tables/table_mammals_with_deamidation.tsv"):
            message.escape("The peptide table file has been deleted.")
        taxonomy="Taxonomy/taxonomy_mammals.tsv"
        peptide_table=["Peptide_tables/table_mammals_with_deamidation.tsv"]

    if taxonomy:
        if not os.path.isfile(taxonomy):
            message.warning("File "+taxonomy+" not found. Ignored.")
            taxonomy=None
        elif  os.path.getsize(taxonomy) == 0:
            message.warning("File "+taxonomy+" is empty. Ignored.")
            taxonomy=None
    if limit:
        if not os.path.isfile(limit):
            message.warning("File "+limit+" not found. No limit applied.")
        elif os.path.getsize(limit) == 0:
            message.warning("File "+limit+" is empty. No limit applied.")
    if neighbour not in range(101):
        neighbour=100
        message.warning("Parameter -n (neighbouring): value is 100")


    if all and neighbour is None:
        message.warning("Ignored parameter: -a (all). This parameter comes with -n (near-optimal solutions).")
        
    if not light:    
        q = (peptide_table, fasta, directory)
        if not (q[0] or q[1] or q[2]):
            message.escape("Missing information for marker peptides (-p, -f or -d).")
        if (q[0] and q[1]) or (q[0] and q[2]) or (q[1] and q[2]) :
            message.escape("Options -p (peptide_table), -f (fasta) and -d (directory of fasta files) are mutually incompatible.")

        if q[1]:
            if not os.path.isfile(fasta):
                message.escape("File "+fasta+" (-f) not found.")
            elif os.path.getsize(fasta) == 0:
                message.escape("File "+fasta+" is empty.")

        if q[2]: 
             if not os.path.isdir(directory):    
                message.escape("Directory "+directory+" (-d) not found.")
                
        if q[0]:
            for pep in peptide_table:
                if not os.path.isfile(pep):
                    message.escape("File "+pep+" not found.")
                elif os.path.getsize(pep) == 0:
                    message.escape("File "+pep+" is empty.")


   
    
    extension=output_file[-4:].lower()
    if extension==".tsv":
        output_file=output_file[:-4]

    output=os.path.join(output_dir, output_file)
    output2=os.path.join(output_dir, "detail_"+output_file+".tsv")
    report_file="report_"+output_file+".txt"
    

    if not light and not peptide_table:
        new_table=os.path.join(output_dir, "table_"+output_file)
    else:
        new_table=None

        
    return (spectra, taxonomy, peptide_table, fasta, directory, limit, error, neighbour, all, output, output2, output_dir, report_file, new_table)

        
def create_report(report_file, output_dir, spectra, list_of_spectra, taxonomy, taxonomy_tree, peptide_table, fasta, directory,  limit, error, neighbour, output, output2, new_table, set_of_markers, WEB):

    report_name=os.path.join(output_dir, report_file)
    sys.stdout=open(report_name, 'w')
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
        spectra_files = str([spectrum.name for spectrum in list_of_spectra])
    else:
        peptide_table_file=peptide_table
        fasta_file=fasta
        taxo_file=taxonomy
        spectra_files = str([spectrum.name for spectrum in list_of_spectra])
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
        print("only optimal solutions")
    else:
        print("up to "+str(neighbour)+"% suboptimality")
    print("  Error margin            : "+str(error), end=" ")
    if error<1:
        print("Da")
    else:
        print("ppm")
    print("")
    if os.path.getsize(os.path.join(output_dir,'warning.log')) > 0:
        print("---------------------------------")
        print("   WARNINGS")
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
    markers.colinearity(set_of_markers)
    markers.check_set_of_markers(set_of_markers)
    if taxonomy is not None:
        print("---------------------------------")
        print("  TAXONOMY")
        print("---------------------------------")
        ta.table_print(taxonomy_tree)
    
    sys.stdout = sys.__stdout__



def main(spectra, taxonomy, peptide_table, fasta, directory, limit, error, neighbour, allsolutions, output, mammals, light, web):
   
    try:

        (spectra, taxonomy, peptide_table, fasta, directory, limit, error, neighbour, allsolutions, output, output2, output_dir, report_file, new_table)=check_and_update_parameters(spectra, taxonomy, peptide_table, fasta, directory, limit, error, neighbour, all, output, mammals, light)

        # parsing spectra files
        list_of_spectra = []
        spectra_dir = os.listdir(spectra)
        for f in spectra_dir:
            file_name = os.path.join(spectra, f)
            spectrum=mass_spectrum.parser(file_name,f)
            if len(spectrum)>0:
                list_of_spectra.append(spectrum)
        if len(list_of_spectra)==0:
                message.escape("No valid spectra found.\n Please refer to the warning.log file for more details.")    
        list_of_spectra.sort(key=lambda x: x.name)
        
            
        # parsing taxonomy    
        if taxonomy is not None:   
            primary_taxonomy=ta.parse_taxonomy_simple_file(taxonomy)
        else:
            primary_taxonomy=None
    
        # parsing models for organisms and applying limits    
        if peptide_table :
            set_of_markers = pt.parse_peptide_tables(peptide_table, limit, primary_taxonomy)
            if not light:
                set_of_markers=compute_masses.add_PTM_or_masses_to_markers(set_of_markers)
        if fasta or directory:
            set_of_sequences = fa.build_set_of_sequences(fasta, directory, limit, primary_taxonomy)
            if len(set_of_sequences)==0:
                 message.escape("Fasta file(s): No valid sequences found.\nPlease refer to the warning.log file to trace back the errors.")
            set_of_markers = markers.sort_and_merge(compute_masses.add_PTM_or_masses_to_markers(seq.in_silico_digestion(set_of_sequences)))
            
        if len(set_of_markers)==0:
            message.escape("No valid peptide marker found.\nPlease refer to the warning.log file to trace back the errors.")
        
        # construction of the secondary taxonomy and suppression of taxid not present in the taxonomy
        set_of_taxid={m.taxid for m in set_of_markers}
        if primary_taxonomy is not None:  
            secondary_taxonomy, lost_taxid=primary_taxonomy.intersection(set_of_taxid)
            set_of_markers= markers.remove_lost_taxid(set_of_markers, lost_taxid)
            for taxid in lost_taxid:
                message.warning("File "+taxonomy+": TaxID "+str(taxid)+" not found. All markers associated to this TaxID are ignored.") 
        else:
            secondary_taxonomy=ta.build_flat_taxonomy(set_of_markers)
        
        if fasta or directory:
            pt.build_peptide_table_from_set_of_markers(set_of_markers, new_table)
        
        create_report(report_file, output_dir, spectra, list_of_spectra, taxonomy, secondary_taxonomy, peptide_table, fasta,  directory, limit, error, neighbour, output, output2, new_table, set_of_markers,web)
        
        # species identification
        assignment.assign_all_spectra(list_of_spectra, set_of_markers, error, taxonomy, secondary_taxonomy, neighbour, allsolutions, output, output2)

        if not web:
            print("")
            print("   Job completed.")
            print("   All results are available in the following files.") 
            print("")
            print(f"   - Assignments       : {output}.tsv")
            print(f"   - More details      : {output2}")
            print(f"   - Report on the run : {os.path.join(output_dir, report_file)}")
            print("")
    # TO DO: add the new peptide table, if necessary

            if os.path.getsize(os.path.join(output_dir, "warning.log")) > 0:
                print("Warnings were raised regarding your inputs.")
                print("Please refer to the warning.log file for details.")

    except message.InputError:
        if not web:
           print("\n   An error occured with your input. Stopping execution.")
           print("   Please refer to the error.log and warning.log files for more details.")
        else:
           pass
 
if __name__ == "__main__":
    main()
