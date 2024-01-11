#!/usr/bin/env python3

"""

   main_assign.py                              

"""


import argparse
import sys
import time
import os

from src import mass_spectrum 
from src import assignment 
from src import peptide_table as pt
from src import sequences as seq
from src import taxonomy as ta 
from src import markers 
from src import fasta_parsing as fa
from src import compute_masses


def escape(message):
    print("\n "+message)
    sys.exit(" Stopping execution.\n\n Type -h for more help.\n")

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
        for pep in peptide_table:
            if not os.path.isfile(pep):
                escape("File "+pep+" not found.")

    output_dir, output_file = os.path.split(output)
    # Ensure the output directory exists. If not, create it.
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

    if not peptide_table:
        new_table=os.path.join(output_dir, "table_"+output_file)
    else:
        new_table=None

    light=mammals
    
    return (spectra, taxonomy, peptide_table, fasta, directory, limit, error, neighbour, all, output, output2, report_path, new_table, light)

        
def create_report(report_name, spectra, list_of_spectra, taxonomy, taxonomy_tree, lost_taxid, peptide_table, fasta, directory,  limit, error, neighbour, output, output2, new_table, set_of_markers):
    
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
    print("  Error margin            : "+str(error), end=" ")
    if error<1:
        print("Da")
    else:
        print("ppm")
    print("")
    print("---------------------------------")
    print("   OUTPUT FILES")
    print("---------------------------------\n")
    print("  Main result file       : "+output)
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
        if len(lost_taxid)>0:
            print ("\nTaxid not found : "+str(lost_taxid))
            print ("The associated markers are removed from the remaining of the analysis.")
    sys.stdout = sys.__stdout__
    
def main(spectra=None, taxonomy=None, peptide_table=None, fasta=None, directory=None, limit=None, error=None, neighbour=None, all=None, output=None,  mammals=None, light=False):

    (spectra, taxonomy, peptide_table, fasta, directory, limit, error, neighbour, allsolutions, output, output2, report_path, new_table, light)=check_and_update_parameters(spectra, taxonomy, peptide_table, fasta, directory, limit, error, neighbour, all, output, mammals)

    # parsing spectra files
    list_of_spectra=[]
    spectra_dir=os.listdir(spectra)
    for f in spectra_dir:
        file_name= os.path.join(spectra, f)
        spectrum=mass_spectrum.parser(file_name,f)
        list_of_spectra.append(spectrum)
        list_of_spectra.sort(key=lambda x: x.name)

    # parsing taxonomy    
    if taxonomy is not None:   
        primary_taxonomy=ta.parse_taxonomy_simple_file(taxonomy)
    else:
        primary_taxonomy=None
        
    # parsing models for organisms and applying limits    
    if peptide_table :
        set_of_markers = pt.parse_peptide_tables(peptide_table, limit, primary_taxonomy)
        if len(set_of_markers)==0:
            sys.exit(" No marker peptide found. Stopping execution.\n")
        if not light:
            set_of_markers=compute_masses.add_PTM_or_masses_to_markers(set_of_markers)
    if fasta or directory:
        set_of_sequences = fa.build_set_of_sequences(fasta, directory, limit, primary_taxonomy)
        if len(set_of_sequences)==0:
            sys.exit(" No sequences found. Stopping execution.\n")
        set_of_markers = compute_masses.add_PTM_or_masses_to_markers(seq.in_silico_digestion(set_of_sequences))

    if len(set_of_markers)==0:
         sys.exit(" No marker peptide found. Stopping execution.\n")
        
    # construction of the secondary taxonomy and suppression of taxid not present in the taxonomy
    set_of_taxid={m.taxid for m in set_of_markers}
    if primary_taxonomy is not None:  
        secondary_taxonomy, lost_taxid=primary_taxonomy.intersection(set_of_taxid)
        set_of_markers= markers.remove_lost_taxid(set_of_markers, lost_taxid)
    else:
        secondary_taxonomy=ta.build_flat_taxonomy(set_of_markers)
        lost_taxid=set()

    if fasta or directory:
        pt.build_peptide_table_from_set_of_markers(set_of_markers, new_table)
        
    create_report(report_path, spectra, list_of_spectra, taxonomy, secondary_taxonomy, lost_taxid, peptide_table, fasta,  directory, limit, error, neighbour, output, output2, new_table, set_of_markers)
        
    # species identification
    assignment.assign_all_spectra(list_of_spectra, set_of_markers, error, taxonomy, secondary_taxonomy, neighbour, allsolutions, output, output2)

    print("")
    print("   Job completed.")
    print("   All results are available in the following files.") 
    print("")
    print(f"   - Assignments       : {output}")
    print(f"   - More details      : {output2}")
    print(f"   - Report on the run : {report_path}")
    print("")


