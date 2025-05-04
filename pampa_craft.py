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
from src import config as conf
from src import supplement
from src import compute_masses
from src import limit as lmt
from src import params_checker
from src import report as rep

class CustomFormatter(argparse.HelpFormatter):
    def add_argument(self, action):
        pass
    def format_help(self):
        custom_paragraph = "\nUsage: pampa_craft  [-h] \n   --allpeptides| --deamidation | --fillin | --homology | --selection \n   [-f FASTA | -d DIRECTORY] [-p PEPTIDE_TABLE] [-s SPECTRA] [-l LIMIT] [-t TAXONOMY] -o OUTPUT \n \nThis script is for the design of custom peptide tables.\nIt should be invoked with one of the following pamaters:\n\n   --allpeptides    Generation of all tryptic peptides from FASTA sequences, possibly filtered by a set of MS spectra. \n   --deamidation    Addition of deamidation modifications to an existing peptide table\n   --fillin         Supplementing a partially filled peptide table: adding missing masses, positions, sequences etc.  \n   --homology       Construction of a new peptide table by homology\n   --selection      Filtration of markers of an existing peptide table with a set of MS spectra. \n\nOptions coming with --allpeptides\n   -f FASTA         Fasta file for new species\n   -d DIRECTORY     Directory containing Fasta files for new species\n   -l LIMIT         Limit file that applies constraints to the set of sequences (tokens GN, OS and OX). OPTIONAL\n   -s SPECTRA       Path to the spectra files. Authorized formats: cvs, mgd, mzML. OPTIONAL\n   -e ERROR         Error margin tolerance. OPTIONAL\n   -o OUTPUT        Path to the output file (new peptide table)\n\nOptions coming with --deamidation\n   -p PEPTIDE_TABLE Peptide table for which deamidation should be added.\n   -l LIMIT         Limit file to apply constraints on the set of markers affected by the modification (token Marker). OPTIONAL\n   -o OUTPUT        Path to  the output file (new peptide table)\n\nOptions coming with --fillin\n   -p PEPTIDE_TABLE Peptide table for which missing information should be completed.\n   -f FASTA         Fasta file for supplementary sequences.OPTIONAL\n   -d DIRECTORY     Directory containing Fasta files for supplementary sequences. OPTIONAL\n   -t TAXONOMY      Path to the taxonomy file needed to add missing taxonomic information. OPTIONAL\n   -o OUTPUT        Path to the output file (new peptide table)\n\nOptions coming with --homology\n   -p PEPTIDE_TABLE [PEPTIDE_TABLE]\n          Peptide table(s) that contain model peptide markers\n   -f FASTA         Fasta file for new species\n   -d DIRECTORY     Directory containing Fasta files for new species\n   -l LIMIT         Limit file that applies constraints on the set of sequences (tokens GN, OS and OX). OPTIONAL\n   -o OUTPUT        Path to the output file (new peptide table)\n\nOptions coming with --selection\n   -p PEPTIDE_TABLE Peptide table to be filtered.\n   -s SPECTRA       Path to the spectra files. Authorized formats: cvs, mgd, mzML.   \n   -e ERROR         Error margin tolerance   \n   -o OUTPUT        Path to the output file (new peptide table)\n\n"
        return custom_paragraph + super(CustomFormatter, self).format_help()


def main():
    parser = argparse.ArgumentParser(formatter_class=CustomFormatter, usage=argparse.SUPPRESS)
    parser.add_argument("--homology",  dest="homology", action='store_true', help="Generate a new table by homology.", required=False)
    parser.add_argument("--deamidation", dest="deamidation", action='store_true', help="Add deamidation to marker masses. -l option can be used to specify the list of involved markers.")
    parser.add_argument("--allpeptides",  dest="allpeptides", action='store_true', help="Generation of all tryptic peptides a from FASTA sequences (specified with either -f or -d).", required=False)
    parser.add_argument("--selection", dest="selection", action='store_true', help="Selection of peptide markers from a set of spectra.", required=False)
    parser.add_argument("--fillin",  dest="fillin", action='store_true', help="Fill in missing information (such as masses, sequences...) to an existing peptide table (specified with -p).", required=False)
    parser.add_argument("-p", dest="peptide_table", nargs='+', help="Peptide table (TSV file). Required with --homology and --fillin.", type=str)
    parser.add_argument("-o", dest="output", help="Output path (should include the output file name)", type=str)
    parser.add_argument("-f", dest="fasta", help="FASTA file that contains new sequences.", type=str)
    parser.add_argument("-d", dest="fasta_dir", help="Directory that contains FASTA files.", type=str)
    parser.add_argument("-s", dest="spectra", help="Directory that contains spectra files (one spectrum per file) for marker filtering", type=str)
    parser.add_argument("-e", dest="resolution", help="Error margin for mass spectrum peaks. Recommended values: 0.01 for MALDI FTICR and 0.1 for MALDI TOF.", type=float)
    parser.add_argument("-l", dest="limit",  help="Limit file (txt)", type=str)
    parser.add_argument("-t", dest="taxonomy", help="Taxonomy file (TSV)", type=str)
    parser.add_argument("--web", dest="web",  action='store_true', help=argparse.SUPPRESS, required=False)
    parser.add_argument("-c", dest="config", help="Config file (json). Default is config.json", type=str, required=False)
    args = parser.parse_args()

    try:
        output_dir=""
        output=""
        report="report.txt"
        web=args.web
        output_dir, output_file, report_file, _ , _= params_checker.logger_and_outputdir_configuration(args.output, " ".join(sys.argv))
        output=os.path.join(output_dir, output_file)
        report=os.path.join(output_dir, report_file)
        rep.create_report_header(" ".join(sys.argv), report)
        (homology, deamidation, allpeptides, fillin, selection, peptide_table, fasta, directory, spectra, error, limit, taxonomy, config) = params_checker.check_and_update_parameters_craft(args.homology, args.deamidation, args.allpeptides, args.fillin, args.selection, args.peptide_table, args.fasta, args.fasta_dir, args.spectra, args.resolution, args.limit, args.taxonomy, args.config)
        
        list_of_constraints = lmt.parse_limits(limit)
        primary_taxonomy = ta.parse_taxonomy_simple_file(taxonomy)
        config_markers = conf.config_markers(config)
        config_headers = conf.config_headers(config)
    
        if homology:
            set_of_markers, _ = pt.parse_peptide_tables(peptide_table, None, None)
            set_of_sequences = fa.build_set_of_sequences(fasta, directory, list_of_constraints, primary_taxonomy)
            config_digestion=conf.config_digestion(config)
            list_of_markers=homo.find_markers_all_sequences(set_of_sequences, set_of_markers, primary_taxonomy, config_digestion)
            list_of_markers=ta.add_taxonomy_ranks(list_of_markers, primary_taxonomy)
            pt.build_peptide_table_from_set_of_markers(list_of_markers,output, config_headers, config_markers)
            rep.create_report_homology(peptide_table, list_of_markers, set_of_sequences, taxonomy, web)
            
        if deamidation:
            set_of_codes=set()
            for constraint in list_of_constraints:
                if 'Deamidation' in constraint:
                    set_of_codes.update(constraint['Deamidation'])
            set_of_markers, list_of_headers=pt.parse_peptide_tables(peptide_table, None, None)
            set_of_new_markers=compute_masses.add_deamidation(set_of_markers, set_of_codes)
            pt.build_peptide_table_from_set_of_markers(set_of_markers.union(set_of_new_markers),output, list_of_headers, config_markers)
            rep.create_report_deamidation(peptide_table, set_of_codes, web)
            
        if allpeptides:
            set_of_sequences = fa.build_set_of_sequences(fasta, directory, list_of_constraints, None)
            if len(set_of_sequences)==0:
                message.escape("No valid sequences found.\n")
            config_digestion=conf.config_digestion(config)
            if spectra is None:
                rep.create_report_allpeptides(set_of_sequences, config_digestion, web)
                set_of_new_markers = compute_masses.add_PTM_or_masses_to_markers(seq.in_silico_digestion(set_of_sequences, config_digestion))
                if len(set_of_new_markers)==0:
                    message.escape("No valid peptide markers found.\n")
            else:
                set_of_markers = compute_masses.add_PTM_or_masses_to_markers(seq.in_silico_digestion(set_of_sequences, config_digestion), True, True)
                if len(set_of_markers)==0:
                    message.escape("No valid peptide markers found.\n")
                final_list_of_spectra=[]
                for f in os.listdir(args.spectra):
                    file_name = os.path.join(spectra, f)
                    list_of_spectra=mass_spectrum.parser(file_name,f)
                    if len(list_of_spectra)>0:
                        final_list_of_spectra.extend(list_of_spectra)
                if len(final_list_of_spectra)==0:
                    message.escape("No valid spectra found.\n Please refer to the warning.log file or the report file for more detail.")
                config_selection=conf.config_selection_peaks(config)
                rep.create_report_allpeptides(set_of_sequences, config_digestion, web, spectra, final_list_of_spectra, error, config_selection)
                minimal_number_of_spectra=max(1, len(final_list_of_spectra)*config_selection)
                set_of_new_markers=marker_filtering.filter_set_of_markers(set_of_markers, final_list_of_spectra, error, minimal_number_of_spectra)
            
            list_of_markers=markers.sort_and_merge(set_of_new_markers)
            pt.build_peptide_table_from_set_of_markers(list_of_markers,output, config_headers)
            
            
        if selection: # not compatible with -f or -d
            # TO DO: add deamidation ?
            set_of_markers, list_of_headers= pt.parse_peptide_tables(peptide_table, None, None)
            set_of_markers=compute_masses.add_PTM_or_masses_to_markers(set_of_markers)
            if len(set_of_markers)==0:
                message.escape("No valid peptide markers found.\n")
            final_list_of_spectra=[]
            for f in os.listdir(spectra):
                file_name = os.path.join(spectra, f)
                list_of_spectra=mass_spectrum.parser(file_name,f)
                if len(list_of_spectra)>0:
                    final_list_of_spectra.extend(list_of_spectra)
            if len(final_list_of_spectra)==0:
                message.escape("No valid spectra found.\n Please refer to the warning.log file or the report file for more detail.")
            config_selection=conf.config_selection_peaks(config)
            minimal_number_of_spectra=max(1, len(final_list_of_spectra)*config_selection)
            set_of_confirmed_markers=marker_filtering.filter_set_of_markers(set_of_markers, final_list_of_spectra, error, minimal_number_of_spectra)
            list_of_markers=markers.sort_and_merge(set_of_confirmed_markers)
            pt.build_peptide_table_from_set_of_markers(set_of_confirmed_markers,output, list_of_headers, config_markers)
            rep.create_report_selection(spectra, final_list_of_spectra, peptide_table, set_of_markers, config_selection, error, web)

        if fillin:
            # to do: check that there is a single peptide table
            set_of_markers, list_of_headers = pt.parse_peptide_tables(peptide_table, list_of_constraints, None, False) # check list_of_constraints here.
            set_of_incomplete_markers, set_of_complete_markers, set_of_incomplete_fields  = supplement.search_for_incomplete_markers(set_of_markers, set(list_of_headers))
            if fasta or directory:
                set_of_sequences = fa.build_set_of_sequences(fasta, directory, list_of_constraints, None)
                config_digestion=conf.config_digestion(config)
                set_of_incomplete_markers=markers.add_sequences_and_positions_to_markers(set_of_incomplete_markers, set_of_sequences)
                if error:
                    set_of_incomplete_markers=markers.check_masses_and_sequences(set_of_incomplete_markers, error)
                    set_of_incomplete_markers=markers.find_sequences_from_mass(set_of_incomplete_markers, set_of_sequences, error)
                set_of_incomplete_markers=supplement.add_digestion_status(set_of_incomplete_markers, set_of_sequences, config_digestion)
            set_of_new_markers=supplement.add_length(compute_masses.add_PTM_or_masses_to_markers(set_of_incomplete_markers))
            set_of_new_markers=ta.supplement_taxonomic_information(set_of_new_markers, primary_taxonomy)
            
            #list_of_markers=list(set_of_new_markers | set_of_complete_markers)
            list_of_markers=markers.sort_and_merge(set_of_new_markers | set_of_complete_markers)
            set_of_markers=ta.add_taxonomy_ranks(list_of_markers, primary_taxonomy)
            pt.build_peptide_table_from_set_of_markers(set_of_markers,output, list_of_headers, config_markers)
            if fasta or directory:
                rep.create_report_supplement(peptide_table, web, taxonomy, list_of_markers, set_of_sequences)
            else:
                rep.create_report_supplement(peptide_table, web, taxonomy)

        rep.create_report_footer(output_dir, output, report)
        
        if not web:
            print("")
            print("Job completed.")
            print("All results are available in the following files.")
            print("")
            print(f"   - New peptide table : {output}")
            print(f"   - Report on the run : {report}")
            print("")

            if os.path.getsize(os.path.join(output_dir, "warning.log")) > 0:
                print("Warnings were raised during the execution.")
                print("Please refer to the warning.log file or the report file for detail.")
        
    except message.InputError:
        rep.create_report_footer(output_dir, output, report)
        if not web:
           print("\n   An error occured with your input. Stopping execution.")
           print("   Please refer to the warning.log file or the "+report+" file for more detail.")
        else:
           pass

if __name__ == "__main__":
    main()
