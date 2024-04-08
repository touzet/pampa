"""
main_build.py

building new markers and new peptide tables
""" 
import os
import time
import sys

# local import
from src import markers
from src import sequences 
from src import known_markers 
from src import fasta_parsing as fa
from src import peptide_table as pt

def escape(message):
    print("\n "+message)
    sys.exit(" Stopping execution.\n\n Type <pampa build -h> for more help.\n")


def check_and_update_parameters(peptide_table, fasta, directory, limit, output):
    """
    Parameters checking and fixing
    """
    
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

    if peptide_table is None:
        escape("Missing parameter: peptide_table (-p)")
    for pep in peptide_table:
            if not os.path.isfile(pep):
                escape("File "+pep+" not found.")    
    if limit:
        if not os.path.isfile(limit):
            escape("File "+limit+" not found (-l).")            
        
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
            
    return (peptide_table, fasta, directory, limit, output, report)


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
    


# a faire : calcul de la masse en option ?

def main(peptide_table=None, fasta=None, directory=None, limit=None, output=None):
    
    (peptide_table, fasta, directory, limit, output, report)=check_and_update_parameters(peptide_table, fasta, directory, limit, output)

    set_of_markers= pt.parse_peptide_tables(peptide_table, None, None)
    set_of_sequences = fa.build_set_of_sequences(fasta, directory, limit, None)      
    set_of_new_markers=known_markers.find_markers_all_sequences(set_of_sequences, set_of_markers)

    pt.build_peptide_table_from_set_of_markers(set_of_new_markers,output,"")

    create_report(report, output, peptide_table, fasta, directory, limit, set_of_markers, set_of_new_markers)

    print("")
    print("   Job completed.")
    print("   All results are available in the two following files.") 
    print("")
    print(f"   - New markers       : {output}")
    print(f"   - Report on the run : {report}")
    print("")
