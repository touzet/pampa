#!/usr/bin/env python3

import os

# local import
from src import mass_spectrum as ms
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
from src import report as rep
            

def main(command_line, spectra, taxonomy, peptide_table, fasta, fasta_dir, limit, deamidation, error, neighbour, allsolutions, output, mammals, placentals, birds, web, config_file, isotopes):
   
    try:
        output_dir=""
        report="report.txt"
        output_dir, output_file, report_file, detail_file, output_json = params_checker.logger_and_outputdir_configuration(output, command_line)
        output=os.path.join(output_dir, output_file)
        report=os.path.join(output_dir, report_file)
        detail=os.path.join(output_dir, detail_file)
        jsonf=os.path.join(output_dir, output_json)
        rep.create_report_header(command_line, report, web)
        
        (spectra, taxonomy, peptide_table, fasta, fasta_dir, limit, deamidation, error, neighbour, allsolutions, config_file) = params_checker.check_and_update_parameters_classify(spectra, taxonomy, peptide_table, fasta, fasta_dir, limit, deamidation, error, neighbour, allsolutions, mammals, placentals, birds, config_file)

        # parsing taxonomy
        full_taxonomy = ta.parse_taxonomy_simple_file(taxonomy)
        # parsing spectra files
        list_of_spectra=ms.parse_spectra_files(spectra)
        # parsing limit file
        list_of_constraints=lmt.parse_limits(limit)
        set_of_codes_for_deamidation = lmt.deamidated_codes(list_of_constraints, deamidation)
        # parsing models for organisms and applying limits    
        if peptide_table :
            set_of_markers, _ = pt.parse_peptide_tables(peptide_table, list_of_constraints, full_taxonomy)
            ta.supplement_taxonomic_information(set_of_markers, full_taxonomy) # To check: see supplement
            set_of_markers=compute_masses.add_PTM_or_masses_to_markers(set_of_markers)
            config_digestion=None
            set_of_sequences=None
            new_table=None
        if fasta or fasta_dir:
            set_of_sequences = fa.build_set_of_sequences(fasta, fasta_dir, list_of_constraints, full_taxonomy)
            config_digestion=config.config_digestion(config_file)
            if len(set_of_sequences)==0:
                 message.escape("Fasta file(s): No valid sequences found.\nPlease refer to the warning.log file to trace back the errors.")
            set_of_markers=compute_masses.add_PTM_or_masses_to_markers(seq.in_silico_digestion(set_of_sequences,config_digestion))
        if len(set_of_markers)==0:
            message.escape("No valid peptide marker found.\nPlease refer to the warning.log file to trace back the errors.")
        set_of_markers.update(compute_masses.add_deamidation(set_of_markers, set_of_codes_for_deamidation))
        set_of_markers=markers.sort_and_merge(set_of_markers)
       
        # parsing taxonomy and filtering markers accordingly
        final_taxonomy = ta.merge_taxonomy(set_of_markers, full_taxonomy)

        if fasta or fasta_dir:
            config_headers = config.config_headers(config_file)
            new_table=os.path.join(output_dir, "table_"+output_file)
            pt.build_peptide_table_from_set_of_markers(set_of_markers, new_table, config_headers)
                
        config_nb_of_peaks=config.config_minimum_number_of_peaks(config_file)
        config_markers=config.config_markers(config_file)

        # report creation
        rep.create_report_classify(spectra, list_of_spectra, taxonomy, final_taxonomy, peptide_table, fasta, fasta_dir, set_of_sequences, set_of_markers, limit, list_of_constraints, deamidation, error, neighbour, all, new_table, config_digestion, config_nb_of_peaks, web)
        rep.create_report_footer(output_dir, output, report, web)

        # species identification
        assignment.assign_all_spectra(list_of_spectra, set_of_markers, error, taxonomy, final_taxonomy, neighbour, allsolutions, config_nb_of_peaks, config_markers, output, detail, jsonf, isotopes)
        
        if not web:
            print("")
            print("   Job completed.")
            print("   All results are available in the following files.") 
            print("")
            print(f"   - Assignments       : {output}")
            print(f"   - More detail       : {detail}")
            print(f"   - Report on the run : {report}")
            print("")
    # TO DO: add the new peptide table, if necessary

            if os.path.getsize(os.path.join(output_dir, "warning.log")) > 0:
                print("Warnings were raised during execution.")
                print("Please refer to the warning.log file or the report file for details.")


    except message.InputError:
        rep.create_report_footer(output_dir, output, report, web)
        if not web:
           print("\n   An error occurred with your input. Stopping execution.")
           print("   Please refer to the warning.log file or the "+report+" file for more detail.")
        else:
           pass
 
if __name__ == "__main__":
    main()
