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
from src import report as rep
            

def main(command_line, spectra, taxonomy, peptide_table, fasta, fasta_dir, limit, deamidation, error, neighbour, allsolutions, output, mammals, web, config_file):
   
    try:
        output_dir=""
        report="report.txt"
        output_dir, output_file, report_file, detail_file, output_json = params_checker.logger_and_outputdir_configuration(output, command_line)
        output=os.path.join(output_dir, output_file)
        report=os.path.join(output_dir, report_file)
        detail=os.path.join(output_dir, detail_file)
        jsonf=os.path.join(output_dir, output_json)
        rep.create_report_header(command_line, report)
        
        (spectra, taxonomy, peptide_table, fasta, fasta_dir, limit, deamidation, error, neighbour, allsolutions, config_file) = params_checker.check_and_update_parameters_classify(spectra, taxonomy, peptide_table, fasta, fasta_dir, limit, deamidation, error, neighbour, allsolutions, mammals, config_file)


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
        primary_taxonomy=ta.parse_taxonomy_simple_file(taxonomy)
            
        # limit file
        list_of_constraints=lmt.parse_limits(limit)
        
        # parsing models for organisms and applying limits    
        if peptide_table :
            set_of_markers, _ = pt.parse_peptide_tables(peptide_table, list_of_constraints, primary_taxonomy)
            set_of_markers=ta.supplement_taxonomic_information(set_of_markers, primary_taxonomy) # To check: see supplement
            set_of_markers=compute_masses.add_PTM_or_masses_to_markers(set_of_markers)
            config_digestion=None
            set_of_sequences=None
        if fasta or fasta_dir:
            set_of_sequences = fa.build_set_of_sequences(fasta, fasta_dir, list_of_constraints, primary_taxonomy)
            config_digestion=config.config_digestion(config_file)
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
            secondary_taxonomy, lost_taxid = primary_taxonomy.intersection(set_of_taxid)
            set_of_markers= markers.remove_lost_taxid(set_of_markers, lost_taxid)
            for taxid in lost_taxid:
                message.warning("File "+taxonomy+": TaxID "+str(taxid)+" not found. All markers associated to this TaxID are ignored.") 
        else:
            secondary_taxonomy=ta.build_flat_taxonomy(set_of_markers)
            
        if fasta or fasta_dir:
            config_headers = config.config_headers(config_file)
            new_table=os.path.join(output_dir, "table_"+output_file)
            pt.build_peptide_table_from_set_of_markers(set_of_markers, new_table, config_headers)
                
        config_nb_of_peaks=config.config_minimum_number_of_peaks(config_file)
        config_markers=config.config_markers(config_file)
        
        rep.create_report_classify(spectra, final_list_of_spectra, taxonomy, secondary_taxonomy, peptide_table, fasta, fasta_dir, set_of_sequences, set_of_markers, limit, deamidation, error, neighbour, all, new_table, config_digestion, config_nb_of_peaks, web)
        
        rep.create_report_footer(output_dir, output, report)

        # species identification
       
        assignment.assign_all_spectra(final_list_of_spectra, set_of_markers, error, taxonomy, secondary_taxonomy, neighbour, allsolutions, config_nb_of_peaks, config_markers, output, detail, jsonf)
        
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
                print("Warnings were raised regarding your inputs.")
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
