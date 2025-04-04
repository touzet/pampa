#!/usr/bin/env python3

import argparse
import os
import time
import sys

# local import
from src import markers
from src import sequences as seq
from src import homology as homo
from src import fasta_parsing as fa
from src import peptide_table as pt
from src import compute_masses
from src import marker_filtering
from src import mass_spectrum
from src import message
from src import taxonomy as ta
from src import config
from src import supplement
from src import compute_masses
from src import limit as lmt

def check_and_update_parameters(homology, deamidation, allpeptides, fillin, selection, peptide_table, fasta, directory, spectra, limit, taxonomy, output):
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

    if deamidation is None and homology is None and allpeptides is None and fillin is None and spectra is None:
         message.escape("Missing parameter: --homology, --allpeptides, --fillin, --deamidation,  or -s.")

    if homology is None and allpeptides is None and fillin is None:
        if not peptide_table:
            message.escape("Missing parameter: -p (peptide table).")
        if fasta:
            message.warning("Unused parameter: -f (fasta file).")
        if directory:
            message.warning("Unused parameter: -d (directory for fasta files)")
            
    if (homology and allpeptides) or (homology and fillin) or (allpeptides and fillin):
        message.escape("Parameters --homology, --allpeptides and --fillin are mutually exclusive.")
    
    extension=output_file[-4:].lower()
    if extension!=".tsv":
        output_file=output_file+".tsv"
    else:
        output_file=output_file[:-4]+".tsv"
    output=os.path.join(output_dir, output_file)
    report_file="report_"+output_file.replace("tsv", "txt")
    report=os.path.join(output_dir, report_file)

    if fillin or homology or deamidation:
        if peptide_table is None:
            message.escape("Missing parameter: -p (peptide table)")
        for pep in peptide_table:
            if not os.path.isfile(pep):
                message.escape("File "+pep+" not found (-p).")
                
    if allpeptides and peptide_table:
        message.warning("Unused parameter: -p (peptide_table)")
                
    if limit:
        if not os.path.isfile(limit):
            message.escape("File "+limit+" not found (-l).")
        if os.path.getsize(limit) == 0:
            message.warning("File "+limit+" is empty.")
            
    if allpeptides or homology:
        q = (fasta, directory)
        if not (q[0] or q[1] ):
            message.escape("Missing target sequences (-f or -d).")
        if (q[0] and q[1])  :
            message.escape("Options -f (fasta file) and -d (directory of fasta files) are mutually exclusive.")
        if q[0]:
            if not os.path.isfile(fasta):
                message.escape("File "+fasta+" not found (-f).")
            if os.path.getsize(fasta) == 0:
                message.escape("File "+fasta+" is empty.")    
        if q[1]:
            if not os.path.isdir(directory):
                message.escape("Directory "+directory+" not found (-d).")
                
      # TO DO : add taxonomy
            
    return (homology, deamidation, allpeptides, fillin, selection, peptide_table, fasta, directory, spectra, limit, taxonomy, output, report, output_dir, report_file)

def create_report_homology(set_of_markers):
    print("  Mode : HOMOLOGY")
    print("  Input peptide table:" )
    print("---------------------------------")
    print(" NEW PEPTIDE TABLE")
    print("---------------------------------\n")
    markers.colinearity(set_of_markers)
    markers.check_set_of_markers(set_of_markers)
    
def create_report_deamidation(peptide_table, set_of_codes):
    print("  Mode : DEAMIDATION")
    print("  Peptide table: "+ str(peptide_table))
    if len(set_of_codes)==0:
        print("  Modified peptide markers: all")
    else:
        print("  Modified peptide markers: "+ str(set_of_codes))
        
def create_report_allpeptides(set_of_sequences, number_of_missed_cleavage, min_length, max_length):
    print("  Mode : ALL PEPTIDES")
    print("  In silico digestion:")
    print("     Enzyme: trypsine")
    print("     Maximal number of missed cleavages: "+str(number_of_missed_cleavage))
    print("     Minimal peptide length: "+str(min_length))
    print("     Maximal peptide length: "+str(max_length))

    print("---------------------------------")
    print("   INPUT SEQUENCES")
    print("---------------------------------")
    for seq in set_of_sequences:
        print ("     "+seq.seqid()+" "+seq.protein()+ " "+seq.taxon_name())
        
def create_report_supplement(peptide_table,list_of_markers=None, set_of_sequences=None):
    print("  Mode : SUPPLEMENT\n")
    print("---------------------------------")
    print("   INPUT PEPTIDE TABLE")
    print("---------------------------------")
    print(str(peptide_table)+"\n")
    if set_of_sequences:
        print("---------------------------------")
        print("   INPUT SEQUENCES")
        print("---------------------------------")
        for seq in set_of_sequences:
            print ("     "+str(seq.seqid())+" "+str(seq.protein())+ " "+str(seq.taxon_name()))
        print("---------------------------------")
        print(" NEW PEPTIDE TABLE")
        print("---------------------------------\n")
        markers.colinearity(set(list_of_markers))
        markers.check_set_of_markers(set(list_of_markers))

def create_report_header(report):
    sys.stdout=open(report, 'w')
    print("============================================\n")
    print("                PAMPA CRAFT\n")
    print("============================================\n")
    print (time.ctime())
    print("")
    print("---------------------------------")
    print("   PARAMETERS")
    print("---------------------------------\n")

def create_report_footer(output_dir, output, report):
    if os.path.getsize(os.path.join(output_dir,'warning.log')) > 0:
        print("---------------------------------")
        print("   WARNINGS")
        print("---------------------------------\n")
        print("  Warnings were raised regarding your inputs.\n")
        with open(os.path.join(output_dir,'warning.log'), 'r') as file:
            for line in file:
                print("  - "+line, end="")
        print("")
    print("")
    print("---------------------------------")
    print("   OUTPUT FILES")
    print("---------------------------------\n")
    print("  Main result file (TSV)   : "+output)
    print("  Report (this file)       : "+report)

    sys.stdout = sys.__stdout__
    

class CustomFormatter(argparse.HelpFormatter):
    def add_argument(self, action):
        pass
    def format_help(self):
        custom_paragraph = "\nUsage: pampa_craft  [-h] \n   --allpeptides| --deamidation | --fillin | --homology | --selection \n   [-f FASTA | -d DIRECTORY] [-p PEPTIDE_TABLE] [-s SPECTRA] [-l LIMIT] [-t TAXONOMY] -o OUTPUT \n \nThis script is for the design of custom peptide tables.\nIt should be invoked with one of the following pamaters:\n\n   --allpeptides    Generation of all tryptic peptides from FASTA sequences, possibly filtered by a set of MS spectra. \n   --deamidation    Addition of deamidation modifications to an existing peptide table\n   --fillin         Supplementing a partially filled peptide table: adding missing masses, positions, sequences etc.  \n   --homology       Construction of a new peptide table by homology\n   --selection      Filtration of markers of an existing peptide table with a set of MS spectra. \n\nOptions coming with --allpeptides\n   -f FASTA         Fasta file for new species\n   -d DIRECTORY     Directory containing Fasta files for new species\n   -l LIMIT         Limit file that applies constraints to the set of sequences (tokens GN, OS and OX). OPTIONAL\n   -s SPECTRA       Path to the spectra files. Authorized formats: cvs, mgd, mzML. OPTIONAL\n   -o OUTPUT        Path to the output file (new peptide table)\n\nOptions coming with --deamidation\n   -p PEPTIDE_TABLE Peptide table for which deamidation should be added.\n   -l LIMIT         Limit file to apply constraints on the set of markers affected by the modification (token Marker). OPTIONAL\n   -o OUTPUT        Path to  the output file (new peptide table)\n\nOptions coming with --fillin\n   -p PEPTIDE_TABLE Peptide table for which missing information should be completed.\n   -f FASTA         Fasta file for supplementary sequences.OPTIONAL\n   -d DIRECTORY     Directory containing Fasta files for supplementary sequences. OPTIONAL\n   -t TAXONOMY      Path to the taxonomy file needed to add missing taxonomic information. OPTIONAL\n   -o OUTPUT        Path to the output file (new peptide table)\n\nOptions coming with --homology\n   -p PEPTIDE_TABLE [PEPTIDE_TABLE]\n          Peptide table(s) that contain model peptide markers\n   -f FASTA         Fasta file for new species\n   -d DIRECTORY     Directory containing Fasta files for new species\n   -l LIMIT         Limit file that applies constraints on the set of sequences (tokens GN, OS and OX). OPTIONAL\n   -o OUTPUT        Path to the output file (new peptide table)\n\nOptions coming with --selection\n   -p PEPTIDE_TABLE Peptide table to be filtered.\n   -s SPECTRA       Path to the spectra files. Authorized formats: cvs, mgd, mzML.   \n   -o OUTPUT        Path to the output file (new peptide table)\n\n"
        return custom_paragraph + super(CustomFormatter, self).format_help()


def main():
    parser = argparse.ArgumentParser(formatter_class=CustomFormatter, usage=argparse.SUPPRESS)
    parser.add_argument("--homology",  dest="homology", action='store_true', help="Generate a new table by homology.", required=False)
    parser.add_argument("--deamidation", dest="deamidation", action='store_true', help="Add deamidation to marker masses. -l option can be used to specify the list of involved markers.")
    parser.add_argument("--allpeptides",  dest="allpeptides", action='store_true', help="Generation of all tryptic peptides a from FASTA sequences (specified with either -f or -d).", required=False)
    parser.add_argument("--selection", dest="selection", action='store_true', help="Selection of peptide markers from a set of spectra.", required=False)
    parser.add_argument("--fillin",  dest="fillin", action='store_true', help="Fill in missing information (such as masses, sequences...) to an existing peptide table (specified with -p).", required=False)
    parser.add_argument("-p", dest="peptide_table",nargs='+', help="Peptide table (TSV file). Required with --homology and --fillin.", type=str)
    parser.add_argument("-o", dest="output", help="Output path (should include the output file name)", type=str)
    parser.add_argument("-f", dest="fasta", help="FASTA file that contains new sequences.", type=str)
    parser.add_argument("-d", dest="directory", help="Directory that contains FASTA files.", type=str)
    parser.add_argument("-s", dest="spectra", help="Directory that contains spectra files (one spectrum per file) for marker filtering", type=str)
    parser.add_argument("-e", dest="resolution", help="Error margin for mass spectrum peaks. Recommended values: 0.01 for maldi FT and 0.1 for maldi TOF.", type=float)
    parser.add_argument("-l", dest="limit",  help="Limit file (txt)", type=str)
    parser.add_argument("-t", dest="taxonomy", help="Taxonomy file (TSV)", type=str)
    parser.add_argument("--web", dest="web",  action='store_true', help=argparse.SUPPRESS, required=False)
    args = parser.parse_args()

    try:
        # to do: add taxonomy 
        (homology, deamidation, allpeptides, fillin, selection, peptide_table, fasta, directory, spectra, limit, taxonomy, output, report, output_dir, report_file)=check_and_update_parameters(args.homology, args.deamidation, args.allpeptides, args.fillin, args.selection, args.peptide_table, args.fasta, args.directory, args.spectra, args.limit, args.taxonomy, args.output)
        
        create_report_header(report)
        list_of_constraints=lmt.parse_limits(limit)
        primary_taxonomy=ta.parse_taxonomy_simple_file(taxonomy)
    
        if homology:
            set_of_markers, _ = pt.parse_peptide_tables(peptide_table, None, None)
            set_of_sequences = fa.build_set_of_sequences(fasta, directory, list_of_constraints, primary_taxonomy)
            if taxonomy:
                primary_taxonomy=ta.parse_taxonomy_simple_file(taxonomy)
            list_of_markers=homo.find_markers_all_sequences(set_of_sequences, set_of_markers, primary_taxonomy)
            pt.build_peptide_table_from_set_of_markers(list_of_markers,output)
            create_report_homology(set_of_markers)
            
        if deamidation:
            set_of_codes=set()
            for constraint in list_of_constraints:
                if 'Deamidation' in constraint:
                    set_of_codes.update(constraint['Deamidation'])
            set_of_markers, list_of_headers=pt.parse_peptide_tables(peptide_table, None, None)
            set_of_new_markers=compute_masses.add_deamidation(set_of_markers, set_of_codes)
            pt.build_peptide_table_from_set_of_markers(set_of_markers.union(set_of_new_markers),output, list_of_headers)
            create_report_deamidation(peptide_table, set_of_codes)
            
        if allpeptides:
            set_of_sequences = fa.build_set_of_sequences(fasta, directory, list_of_constraints, None)
            if len(set_of_sequences)==0:
                message.escape("No valid sequences found.\n")
            data=config.parse_config_file()
            if spectra is None:
                set_of_new_markers = compute_masses.add_PTM_or_masses_to_markers(seq.in_silico_digestion(set_of_sequences, data["number_of_missed_cleavages"], data["min_peptide_length"], data["max_peptide_length"]))
                if len(set_of_new_markers)==0:
                    message.escape("No valid peptide markers found.\n")
            else:
                set_of_markers = compute_masses.add_PTM_or_masses_to_markers(seq.in_silico_digestion(set_of_sequences, data["number_of_missed_cleavages"], data["min_peptide_length"], data["max_peptide_length"]), True, True)
                if len(set_of_markers)==0:
                    message.escape("No valid peptide markers found.\n")
                final_list_of_spectra=[]
                for f in os.listdir(args.spectra):
                    file_name = os.path.join(spectra, f)
                    list_of_spectra=mass_spectrum.parser(file_name,f)
                    if len(list_of_spectra)>0:
                        final_list_of_spectra.extend(list_of_spectra)
                if len(final_list_of_spectra)==0:
                    message.escape("No valid spectra found.\n Please refer to the warning.log file for more detail.")
                minimal_number_of_spectra=max(1, 1*len(final_list_of_spectra)/5)
                set_of_new_markers=marker_filtering.filter_set_of_markers(set_of_markers, final_list_of_spectra, args.resolution, minimal_number_of_spectra)
            
            list_of_markers=markers.sort_and_merge(set_of_new_markers)
            pt.build_peptide_table_from_set_of_markers(list_of_markers,output)
            create_report_allpeptides(set_of_sequences, data["number_of_missed_cleavages"], data["min_peptide_length"], data["max_peptide_length"])

                     
        if selection: # not compatible with -f or -d
            set_of_markers, list_of_headers= pt.parse_peptide_tables(peptide_table, None, None)
            if len(set_of_markers)==0:
                message.escape("No valid markers found.\n")
            ####
            final_list_of_spectra=[]
            for f in os.listdir(args.spectra):
                file_name = os.path.join(spectra, f)
                list_of_spectra=mass_spectrum.parser(file_name,f)
                if len(list_of_spectra)>0:
                    final_list_of_spectra.extend(list_of_spectra)
            if len(final_list_of_spectra)==0:
                message.escape("No valid spectra found.\n Please refer to the warning.log file for more detail.")
            minimal_number_of_spectra=max(1, 2*len(final_list_of_spectra)/3)
            set_of_confirmed_markers=marker_filtering.filter_set_of_markers(set_of_markers, final_list_of_spectra, args.resolution, minimal_number_of_spectra)
            list_of_markers=markers.sort_and_merge(set_of_confirmed_markers)
            pt.build_peptide_table_from_set_of_markers(set_of_confirmed_markers,output, list_of_headers)

        if fillin:
            # to do: check that there is a single peptide table
            set_of_markers, list_of_headers = pt.parse_peptide_tables(peptide_table, list_of_constraints, None, False) # check list_of_constraints here.
            set_of_incomplete_markers, set_of_complete_markers, set_of_incomplete_fields  = supplement.search_for_incomplete_markers(set_of_markers, set(list_of_headers))
            if fasta or directory:
                set_of_sequences = fa.build_set_of_sequences(fasta, directory, list_of_constraints, None)
                set_of_incomplete_markers=markers.add_sequences_and_positions_to_markers(set_of_incomplete_markers, set_of_sequences)
                if args.resolution:
                    set_of_incomplete_markers=markers.check_masses_and_sequences(set_of_incomplete_markers, args.resolution)
                    set_of_incomplete_markers=markers.find_sequences_from_mass(set_of_incomplete_markers, set_of_sequences, args.resolution)
                set_of_incomplete_markers=supplement.add_digestion_status(set_of_incomplete_markers, set_of_sequences)
            set_of_new_markers=compute_masses.add_PTM_or_masses_to_markers(set_of_incomplete_markers)
            set_of_new_markers=ta.supplement_taxonomic_information(set_of_new_markers, primary_taxonomy)
            
            #list_of_markers=list(set_of_new_markers | set_of_complete_markers)
            list_of_markers=markers.sort_and_merge(set_of_new_markers | set_of_complete_markers)
            set_of_markers=ta.add_taxonomy_ranks(list_of_markers, primary_taxonomy)
            pt.build_peptide_table_from_set_of_markers(set_of_markers,output, list_of_headers)
            if fasta or directory:
                create_report_supplement(peptide_table, list_of_markers, set_of_sequences)
            else:
                create_report_supplement(peptide_table)

        create_report_footer(output_dir, output, report)
        
        print("")
        print("   Job completed.")
        print("   All results are available in the following files.")
        print("")
        print(f"   - New peptide table : {output}")
        print(f"   - Report on the run : {report}")
        print("")
    # TO DO: add the new peptide table, if necessary

        if os.path.getsize(os.path.join(output_dir, "warning.log")) > 0:
            print("Warnings were raised during the execution.")
            print("Please refer to the warning.log file for detail.")
        
    except message.InputError:
        if not args.web:
           print("\n   An error occured with your input. Stopping execution.")
           print("   Please refer to the warning.log file for detail.")
        else:
           pass

if __name__ == "__main__":
    main()
