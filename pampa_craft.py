#!/usr/bin/env python3

import argparse
import os
import time
import sys

# local import
from src import markers
from src import sequences as seq
from src import known_markers 
from src import fasta_parsing as fa
from src import peptide_table as pt
from src import compute_masses

def escape(message):
    print("\n "+message)
    sys.exit(" Stopping execution.\n\n Type <pampa_craft -h> for more help.\n")


def check_and_update_parameters(homology, denovo, fillin, peptide_table, fasta, directory, limit, output):
    """
    Parameters checking and fixing
    """

    if homology is None and denovo is None and fillin is None:
         escape("Missing parameter: --homology, --denovo or --fillin.")

    if (homology and denovo) or (homology and fillin) or (denovo and fillin):
        escape("Parameters --homology, --denovo and --fillin are mutually exclusive.")
    
    if output is None:
        escape("Missing parameter: output (-o).")
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
    report_file="report_"+output_file.replace("tsv", "txt")
    report=os.path.join(output_dir, report_file)

    if fillin or homology:
        if peptide_table is None:
            escape("Missing parameter: peptide_table (-p)")
        for pep in peptide_table:
            if not os.path.isfile(pep):
                escape("File "+pep+" not found.")
                
    if denovo and peptide_table:
        print("Warning: unused parameter, -p "+peptide_table)
                
    if limit:
        if not os.path.isfile(limit):
            escape("File "+limit+" not found (-l).")
            
    if denovo or homology:    
        q = (fasta, directory)
        if not (q[0] or q[1] ):
            escape("Missing target sequences (-f or -d).")
        if (q[0] and q[1])  :
            escape("Options -f (fasta file) and -d (directory of fasta files) are mutually exclusive.")
        if q[0]:
            if not os.path.isfile(fasta):
                escape("File "+fasta+" not found (-f).")
        if q[1]:
            if not os.path.isdir(directory):
                escape("Directory "+directory+" not found (-d).")
            
    return (homology, denovo, fillin, peptide_table, fasta, directory, limit, output, report)


def create_report(report, output, peptide_table, fasta,  directory, limit,  set_of_markers, set_of_new_markers):
    
    sys.stdout=open(report, 'w')
    print (time.ctime())
    print("")
    print("---------------------------------")
    print("   PARAMETERS")
    print("---------------------------------\n")
    if peptide_table:
        print("  Peptide table for model markers : ", end="")
        for file in peptide_table:
            print(file, end=" ")
        print("")
    if fasta:
        print("  Fasta file for new markers      : "+str(fasta))
    else:
        print("  Fasta file for new markers      : "+str(directory))
    if limit:
        print("  Limited to                      : "+ limit)
        
    print("")
    print("---------------------------------")
    print("   OUTPUT FILES")
    print("---------------------------------\n")
    print("  Main result file (TSV)   : "+output)
    print("  Report (this file)       : "+report)

    print("")
    print("---------------------------------")
    print("  NEW MARKERS FOR FASTA FILES ")
    print("---------------------------------\n")
    markers.colinearity(set_of_new_markers)
    markers.check_set_of_markers(set_of_new_markers)
    
    sys.stdout = sys.__stdout__
    

class CustomFormatter(argparse.HelpFormatter):
    def add_argument(self, action):
        pass
    def format_help(self):
        custom_paragraph = "\nUsage: pampa_craft  [-h] \n   --homology | --denovo | --fillin \n   [-f FASTA | -d DIRECTORY] [-p PEPTIDE_TABLE] [-l LIMIT] -o OUTPUT \n \nThis script is for the design of custom peptide tables.\n\n   --homology    Generate a new peptide table by homology\n   --denovo      De novo generation of a peptide table from Fasta sequences\n   --fillin      Add missing masses to an existing peptide table\n\nOptions coming with --homology\n   -p PEPTIDE_TABLE [PEPTIDE_TABLE]\n   Peptide table(s) that contain model peptide markers\n   -f FASTA      Fasta file for new species\n   -d DIRECTORY  Directory containing Fasta files for new species\n   -l LIMIT      Limit file that \n   -o OUTPUT     Path to the output file (new peptide table)\n\nOptions coming with --denovo\n   -f FASTA      Fasta file for new species\n   -d DIRECTORY  Directory containing Fasta files for new species\n   -l LIMIT      Limit file that \n   -o OUTPUT     Path to the output file (new peptide table)\n\nOptions coming with --fillin\n   -p PEPTIDE_TABLE\n   Peptide table for which masses should be completed.\n   -o OUTPUT     Path to the output file (new peptide table)\n"
        return custom_paragraph + super(CustomFormatter, self).format_help()


def main():
    parser = argparse.ArgumentParser(formatter_class=CustomFormatter, usage=argparse.SUPPRESS)
    parser.add_argument("--homology",  dest="homology", action='store_true', help="Generate a new table by homology.", required=False)  
    parser.add_argument("--denovo",  dest="denovo", action='store_true', help="De novo generation of a peptide table from Fasta sequences (specified with either -f or -d).", required=False)
    parser.add_argument("--fillin",  dest="fillin", action='store_true', help="Add missing masses to an existing peptide table (specified with -p).", required=False)
    parser.add_argument("-p", dest="peptide_table",nargs='+', help="Peptide table (TSV file). This option with required with --homolgy and --fillin.", type=str)
    parser.add_argument("-o", dest="output", help="Output path (should include the output file name)", type=str)
    parser.add_argument("-f", dest="fasta", help="FASTA file that contains new sequences", type=str)
    parser.add_argument("-d", dest="directory", help="Directory that contains FASTA files", type=str)
    parser.add_argument("-l", dest="limit",  help="Limit file (txt)", type=str)
  
    args = parser.parse_args()

    (homology, denovo, fillin, peptide_table, fasta, directory, limit, output,report)=check_and_update_parameters(args.homology, args.denovo, args.fillin, args.peptide_table, args.fasta, args.directory, args.limit, args.output)
    
    if homology:
        set_of_markers= pt.parse_peptide_tables(peptide_table, None, None)
        set_of_sequences = fa.build_set_of_sequences(fasta, directory, limit, None)      
        set_of_new_markers=known_markers.find_markers_all_sequences(set_of_sequences, set_of_markers)
        pt.build_peptide_table_from_set_of_markers(set_of_new_markers,output,"")

    if denovo:
        set_of_sequences = fa.build_set_of_sequences(fasta, directory, limit, None)
        if len(set_of_sequences)==0:
            sys.exit(" No sequences found. Stopping execution.\n")
        set_of_markers = compute_masses.add_PTM_or_masses_to_markers(seq.in_silico_digestion(set_of_sequences))
        if len(set_of_markers)==0:
            sys.exit(" No marker peptide found. Stopping execution.\n")
        pt.build_peptide_table_from_set_of_markers(set_of_markers, output)

    if fillin:
        set_of_markers = pt.parse_peptide_tables(peptide_table, None, None)
        set_of_new_markers=compute_masses.add_PTM_or_masses_to_markers(set_of_markers)
        pt.build_peptide_table_from_set_of_markers(set_of_new_markers,output)


if __name__ == "__main__":
    main()
