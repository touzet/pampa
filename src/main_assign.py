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


def check_parameters(SPECTRA, taxonomy, peptide_table, fasta, Entry, directory, resolution, outfile, amplitude):
    """
    Parameters checking and fixing
    """

    if SPECTRA is None:
        sys.exit("Missing parameter: SPECTRA")
    if not os.path.isdir(SPECTRA):    
        sys.exit("\n Directory "+SPECTRA+" not found. \n Stopping execution.\n")
    if resolution is None:
        sys.exit("\n Missing parameter: resolution. \n Stopping execution.\n")
    if outfile is None:
        sys.exit("\n Missing parameter: output.  \n Stopping execution.\n")
    

    if taxonomy:
        if not os.path.isfile(taxonomy):
            sys.exit("\n File "+taxonomy+" not found. Stopping execution.\n")

    if amplitude not in range(101):
        amplitude=100
        
    q = (peptide_table, fasta, Entry, directory)
    if not (q[0] or q[1] or q[2]):
        sys.exit("\n Missing information for marker peptides (-p, -f or -E).\n")
    if q[0] and q[1]:
        sys.exit("\n Options -p (peptide_table) and  -f (fasta) are mutually incompatible. \n")
    if q[0] and q[2]:
        sys.exit("\n Options -p (peptide_table) and  -E  are mutually incompatible. \n")
    if  q[1] and q[2]:
        print ("\n Options -f (fasta) and  -E  are mutually incompatible. \n")
        sys.exit(-1)
    if q[3] and not q[2]:
        sys.exit("\n Missing options: -E . \n")

    if q[0]:
        if not os.path.isfile(peptide_table):
           sys.exit("\n File "+peptide_table+" not found. Stopping execution.\n")

    if directory is None:
        dir='.'


def create_report(report_name,SPECTRA, list_of_spectra, taxonomy, taxonomy_tree, peptide_table, fasta, Entry, directory, resolution, outfile, amplitude, set_of_markers):
    
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
        print("  Fasta file (markers)    : "+str(Entry))
    print("  Taxonomy                : "+str(taxonomy))
    print("  Mass spectra            : "+SPECTRA)
    print("  Amplitude               : "+str(amplitude)+"%")
    print("  Resolution              : "+str(resolution))
    print("")
    print("---------------------------------")
    print("   OUTPUT FILES")
    print("---------------------------------\n")
    print("  Main result file         : "+outfile)
    print("  More details             : detail_"+outfile)
    if not peptide_table:
        print("  Table of peptides        : table_"+outfile)
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
    if fasta or Entry:
       print("---------------------------------")
       print("  FASTA SEQUENCES")
       print("---------------------------------\n")
    check_set_of_markers(set_of_markers)
    if taxonomy is not None:
        print("---------------------------------")
        print("  TAXONOMY")
        print("---------------------------------")
        table_print(taxonomy_tree)
    sys.stdout = sys.__stdout__
    
def main(SPECTRA=None, taxonomy=None, peptide_table=None, fasta=None, Entry=None, directory=None, resolution=None, outfile=None, amplitude=None):

    check_parameters(SPECTRA, taxonomy, peptide_table, fasta, Entry, directory, resolution, outfile, amplitude)
    
    extension=outfile[-4:].lower()
    if extension!=".tsv":
        outfile_name=outfil+".tsv"
    else:
        outfile=outfile[:-4]+".tsv"
    report_name="report_"+outfile.replace("tsv", "txt")

    list_of_spectra=[]
    spectra_dir=os.listdir(SPECTRA)
    for f in spectra_dir:
        file_name= join(SPECTRA, f)
        spectrum=mass_spectrum.parser(file_name,f)
        list_of_spectra.append(spectrum)
        list_of_spectra.sort(key=lambda x: x.name)

    if peptide_table:
        set_of_markers = parse_peptide_table(peptide_table)
    if fasta or Entry:
        if fasta:
           set_of_sequences= build_set_of_sequences_from_fasta_file(fasta,True)
        else: 
           set_of_sequences=build_set_of_sequences_from_TSV_fasta_dir(Entry,directory)
        set_of_markers=add_PTM_or_masses_to_markers(in_silico_digestion(set_of_sequences,1,12, False))
        build_peptide_table_from_set_of_markers(set_of_markers,"table_"+outfile,"")

    mass_taxid_name_list, set_of_marker_names, set_of_taxid, taxid_to_taxidname = build_markers_dictionnaries_from_set_of_markers(set_of_markers)
    list_of_marker_names=list(set_of_marker_names)
    list_of_marker_names.sort()
    list_of_taxid=list(set_of_taxid)
    list_of_taxid.sort()
   
    ## construction of the taxonomy  
    if taxonomy is not None:   
        primary_taxonomy=parse_taxonomy_simple_file(taxonomy)
        secondary_taxonomy=primary_taxonomy.intersection(set_of_taxid)
    else:
        secondary_taxonomy=build_flat_taxonomy(set_of_markers)

    create_report(report_name,SPECTRA, list_of_spectra, taxonomy, secondary_taxonomy, peptide_table, fasta, Entry, directory, resolution, outfile, amplitude, set_of_markers)
        
    ## species identification
    assign_all_spectra(list_of_spectra, mass_taxid_name_list, list_of_taxid, list_of_marker_names, resolution, taxonomy, secondary_taxonomy, SPECTRA, "XXXXXX", amplitude, outfile)

    print("")
    print("   Job completed.")
    print("   All results are available in the three following files.") 
    print("")
    print("   - Assignments       : "+outfile)
    print("   - More details      : detail_"+outfile)
    print("   - Report on the run : " + report_name)
    print("")


