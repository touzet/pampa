#!/usr/bin/env python3

import argparse
import os
import sys

# local import
from src import markers
from src import sequences as seq
from src import homology as homo
from src import fasta_parsing as fa
from src import peptide_table as pt
from src import marker_filtering
from src import mass_spectrum as ms
from src import message
from src import taxonomy as taxo
from src import config as conf
from src import supplement
from src import compute_masses
from src import limit as lmt
from src import params_checker
from src import matrix
from src import gamma
from src import reconstruct
from src import report as rep
from src import conserved_markers
from src import neighbour
from src import pareto


class CustomFormatter(argparse.HelpFormatter):
    def add_argument(self, action):
        pass
    def format_help(self):
        custom_paragraph = ("\nUsage: pampa_craft  [-h] \n   --"
                            "allpeptides| --deamidation | --fillin | --homology | --selection | --reconstruction \n   [-f FASTA | -d DIRECTORY] [-p PEPTIDE_TABLE] [-s SPECTRA] [-l LIMIT] [-t TAXONOMY] -o OUTPUT \n \nThis script is for the design of custom peptide tables.\nIt should be invoked with one of the following parameters:\n\n   --allpeptides    Generation of all tryptic peptides from FASTA sequences, possibly filtered by a set of MS spectra. \n   --deamidation    Addition of deamidation modifications to an existing peptide table\n   --fillin         Supplementing a partially filled peptide table: adding missing masses, positions, sequences etc.  \n   --homology       Construction of a new peptide table by homology\n   --selection      Filtration of markers of an existing peptide table with a set of MS spectra. \n   --reconstruction Inference of peptide sequences from MS spectra without genetic data\n\nOptions coming with --allpeptides\n   -f FASTA         Fasta file for new species\n   -d DIRECTORY     Directory containing Fasta files for new species\n   -l LIMIT         Limit file that applies constraints to the set of sequences (tokens GN, OS and OX). OPTIONAL\n   -s SPECTRA       Path to the spectra files. Authorized formats: cvs, mgd, mzML. OPTIONAL\n   -e ERROR         Error margin tolerance. OPTIONAL\n   -o OUTPUT        Path to the output file (new peptide table)\n\nOptions coming with --deamidation\n   -p PEPTIDE_TABLE Peptide table for which deamidation should be added.\n   -l LIMIT         Limit file to apply constraints on the set of markers affected by the modification (token Marker). OPTIONAL\n   -o OUTPUT        Path to  the output file (new peptide table)\n\nOptions coming with --fillin\n   -p PEPTIDE_TABLE Peptide table for which missing information should be completed.\n   -f FASTA         Fasta file for supplementary sequences.OPTIONAL\n   -d DIRECTORY     Directory containing Fasta files for supplementary sequences. OPTIONAL\n   -t TAXONOMY      Path to the taxonomy file needed to add missing taxonomic information. OPTIONAL\n   -o OUTPUT        Path to the output file (new peptide table)\n\nOptions coming with --homology\n   -p PEPTIDE_TABLE [PEPTIDE_TABLE]\n          Peptide table(s) that contain model peptide markers\n   -f FASTA         Fasta file for new species\n   -d DIRECTORY     Directory containing Fasta files for new species\n   -l LIMIT         Limit file that applies constraints on the set of sequences (tokens GN, OS and OX). OPTIONAL\n   -o OUTPUT        Path to the output file (new peptide table)\n\nOptions coming with --selection\n   -p PEPTIDE_TABLE Peptide table to be filtered.\n   -s SPECTRA       Path to the spectra files. Authorized formats: cvs, mgd, mzML.   \n   -e ERROR         Error margin tolerance   \n   -o OUTPUT        Path to the output file (new peptide table)\n\n")
        return custom_paragraph + super(CustomFormatter, self).format_help()


def main():
    parser = argparse.ArgumentParser(formatter_class=CustomFormatter, usage=argparse.SUPPRESS)
    parser.add_argument("--homology",  dest="homology", action='store_true', help="Generate a new table by homology.", required=False)
    parser.add_argument("--deamidation", dest="deamidation", action='store_true', help="Add deamidation to marker masses. -l option can be used to specify the list of involved markers.")
    parser.add_argument("--allpeptides",  dest="allpeptides", action='store_true', help="Generation of all tryptic peptides a from FASTA sequences (specified with either -f or -d).", required=False)
    parser.add_argument("--selection", dest="selection", action='store_true', help="Selection of peptide markers from a set of spectra.", required=False)
    parser.add_argument("--fillin",  dest="fillin", action='store_true', help="Fill in missing information (such as masses, sequences...) to an existing peptide table (specified with -p).", required=False)
    parser.add_argument("--reconstruction",  dest="reconstruction", action='store_true', help="Reconstruct peptide sequences from MS spectra", required=False)
    parser.add_argument("--placentals", dest="placentals", action='store_true', help="Use placental taxonomic model",required=False)
    parser.add_argument("--mammals", dest="mammals", action='store_true', help="Use mammal taxonomic model", required=False)
    parser.add_argument("--birds", dest="birds", action='store_true', help="Use bird taxonomic model", required=False)
    parser.add_argument("--custom", dest="custom", action='store_true', help="Use custom, user-defined taxonomic model",required=False)
    parser.add_argument("-p", dest="peptide_table", nargs='+', help="Peptide table (TSV file). Required with --homology and --fillin.", type=str)
    parser.add_argument("-o", dest="output", help="Output path (should include the output file name)", type=str)
    parser.add_argument("-f", dest="fasta", help="FASTA file that contains new sequences.", type=str)
    parser.add_argument("-d", dest="fasta_dir", help="Directory that contains FASTA files.", type=str)
    parser.add_argument("-s", dest="spectra", help="Directory that contains spectra files (one spectrum per file) for marker filtering", type=str)
    parser.add_argument("-e", dest="resolution", help="Error margin for mass spectrum peaks. Recommended values: 0.01 for MALDI FTICR and 0.1 for MALDI TOF.", type=float)
    parser.add_argument("-l", dest="limit",  help="Limit file (txt)", type=str)
    parser.add_argument("-t", dest="taxonomy", help="Taxonomy file (TSV)", type=str)
    parser.add_argument("-x", dest="target", help="Target species (or clade) for peptide sequence reconstruction: taxid or scientific name.", type=str)
    parser.add_argument("-X", dest="targetfile", help="Target species File (TSV)  for peptide sequence reconstruction: name of spectral_file taxid.", type=str)
    parser.add_argument("--web", dest="web",  action='store_true', help=argparse.SUPPRESS, required=False)
    parser.add_argument("-c", dest="config", help="Config file (json). Default is config.json", type=str, required=False)
    args = parser.parse_args()

    output_dir = ""
    output = ""
    report = "report.txt"
    web = args.web
    try:
        output_dir, output_file, report_file, _ , _= params_checker.logger_and_outputdir_configuration(args.output, " ".join(sys.argv))
        output=os.path.join(output_dir, output_file)
        report=os.path.join(output_dir, report_file)
        rep.create_report_header(" ".join(sys.argv), report, web)
        (homology, deamidation, allpeptides, fillin, selection, reconstruction, peptide_table, fasta, fasta_dir,
         spectra, error, limit, taxonomy, config, placentals, birds, custom,
         target,targetfile) = params_checker.check_and_update_parameters_craft(args.homology, args.deamidation, args.allpeptides,
                                                                    args.fillin, args.selection, args.reconstruction,
                                                                    args.peptide_table, args.fasta, args.fasta_dir,
                                                                    args.spectra, args.resolution, args.limit,
                                                                    args.taxonomy, args.config, args.placentals,
                                                                    args.birds, args.mammals, args.custom, args.target, args.targetfile)

        list_of_constraints = lmt.parse_limits(limit)
        set_of_codes_for_deamidation = lmt.deamidated_codes(list_of_constraints, deamidation)
        full_taxonomy = taxo.parse_taxonomy_simple_file(taxonomy)
        config_markers = conf.config_markers(config)
        config_headers = conf.config_headers(config)
        set_of_sequences=set()

        if homology:
            set_of_markers, list_of_headers = pt.parse_peptide_tables(peptide_table, None, None)
            headers={header.lower() for header in list_of_headers}
            set_of_sequences = fa.build_set_of_sequences(fasta, fasta_dir, list_of_constraints, full_taxonomy)
            config_digestion=conf.config_digestion(config)
            list_of_new_markers=homo.find_markers_all_sequences(set_of_sequences, set_of_markers, full_taxonomy, config_digestion)
            if taxonomy:
                taxo.supplement_taxonomic_information(list_of_new_markers, full_taxonomy)
                taxo.add_taxonomy_ranks(list_of_new_markers, full_taxonomy, headers)
            supplement.add_marker_names(list_of_new_markers)
            pt.build_peptide_table_from_set_of_markers(list_of_new_markers,output, config_headers, config_markers)
            rep.create_report_homology(peptide_table, set_of_markers, list_of_new_markers, fasta, fasta_dir, set_of_sequences, taxonomy, config_digestion, limit, list_of_constraints, deamidation, web)
            deamidation=False

        if allpeptides:
            set_of_sequences = fa.build_set_of_sequences(fasta, fasta_dir, list_of_constraints, full_taxonomy)
            if len(set_of_sequences)==0:
                message.escape("No valid sequences found.\n")
            config_digestion=conf.config_digestion(config)
            if spectra is None:
                rep.create_report_allpeptides(fasta, fasta_dir, set_of_sequences, config_digestion, limit, list_of_constraints, deamidation, web)
                set_of_new_markers = compute_masses.add_PTM_or_masses_to_markers(seq.in_silico_digestion(set_of_sequences, config_digestion))
                if len(set_of_new_markers)==0:
                    message.escape("No valid peptide markers found.\n")
                set_of_new_markers.update(compute_masses.add_deamidation(set_of_new_markers, set_of_codes_for_deamidation))
                supplement.add_marker_comment(set_of_new_markers, "In silico digestion.")
            else:
                set_of_markers = compute_masses.add_PTM_or_masses_to_markers(seq.in_silico_digestion(set_of_sequences, config_digestion), True, True)
                set_of_markers.update(compute_masses.add_deamidation(set_of_markers, set_of_codes_for_deamidation))
                if len(set_of_markers)==0:
                    message.escape("No valid peptide markers found.\n")
                list_of_spectra=ms.parse_spectra_files(spectra)
                config_selection=conf.config_selection_peaks(config)
                rep.create_report_allpeptides(fasta, fasta_dir, set_of_sequences, config_digestion, limit, list_of_constraints, deamidation, web, spectra, list_of_spectra, error, config_selection)
                minimal_number_of_spectra=max(1, len(list_of_spectra)*config_selection)
                set_of_new_markers=marker_filtering.filter_set_of_markers(set_of_markers, list_of_spectra, error, minimal_number_of_spectra)
                supplement.add_marker_comment(set_of_new_markers, "In silico digestion, and MALDI filtering.")
            list_of_markers=markers.sort_and_merge(set_of_new_markers)
            #supplement.add_marker_names(list_of_markers)
            pt.build_peptide_table_from_set_of_markers(list_of_markers,output, config_headers)
            deamidation=False

        if selection: # not compatible with -f or -d
            # TO DO: add deamidation ?
            set_of_markers, list_of_headers= pt.parse_peptide_tables(peptide_table, None, None)
            set_of_markers=compute_masses.add_PTM_or_masses_to_markers(set_of_markers)
            if len(set_of_markers)==0:
                message.escape("No valid peptide markers found.\n")
            list_of_spectra=ms.parse_spectra_files(spectra)
            config_selection=conf.config_selection_peaks(config)
            minimal_number_of_spectra=max(1, len(list_of_spectra)*config_selection)
            set_of_confirmed_markers=marker_filtering.filter_set_of_markers(set_of_markers, list_of_spectra, error, minimal_number_of_spectra)
            #list_of_markers=markers.sort_and_merge(set_of_confirmed_markers)
            pt.build_peptide_table_from_set_of_markers(set_of_confirmed_markers,output, list_of_headers, config_markers)
            rep.create_report_selection(spectra, list_of_spectra, peptide_table, set_of_markers, config_selection, error, web)

        if fillin:
            # to do: check that there is a single peptide table
            set_of_markers, user_headers = pt.parse_peptide_tables(peptide_table, list_of_constraints, None, False) # check list_of_constraints here.
            #supplement.add_magic_numbers(set_of_markers)
            headers={header.lower() for header in user_headers}
            rep.print_markers_with_hidden_numbers(set_of_markers, peptide_table, web)
            if fasta or fasta_dir:
                set_of_sequences = fa.build_set_of_sequences(fasta, fasta_dir, list_of_constraints, None)
                set_of_taxid = {seq.taxid() for seq in set_of_sequences}
                if full_taxonomy is not None:
                    small_taxonomy, _ = full_taxonomy.intersection(set_of_taxid)
                else:
                    small_taxonomy = None
                config_digestion=conf.config_digestion(config)
                taxonomy_ranks=["Family", "family", "Order", "order", "Genus", "genus"] #config_taxonomy
                set_of_markers=markers.supplement_markers(set_of_markers, set_of_sequences, error, small_taxonomy, taxonomy_ranks, config_digestion)
                supplement.add_digestion_status(set_of_markers, set_of_sequences, config_digestion)

            set_of_markers=compute_masses.add_PTM_or_masses_to_markers(set_of_markers)
            supplement.add_length(set_of_markers)
            supplement.add_marker_names(set_of_markers)
            taxo.supplement_taxonomic_information(set_of_markers, full_taxonomy)
            supplement.add_taxid(set_of_markers, set_of_sequences, full_taxonomy)
            #list_of_markers=list(set_of_new_markers | set_of_complete_markers)
            list_of_markers=markers.sort_and_merge(set_of_markers)
            taxo.add_taxonomy_ranks(list_of_markers, full_taxonomy, headers)
            new_headers=pt.merge_headers(user_headers, config_headers, set_of_markers) # !! list_of_markers
            pt.build_peptide_table_from_set_of_markers(list_of_markers, output, new_headers)
            if fasta or fasta_dir:
                rep.create_report_supplement(peptide_table, fasta, fasta_dir, set_of_markers, web, taxonomy, set_of_sequences)
            else:
                rep.create_report_supplement(peptide_table, fasta, fasta_dir, set_of_markers, web, taxonomy)

        if reconstruction:
            # initialisation des spectres
            list_of_all_spectra = ms.parse_spectra_files(spectra)
            list_of_target_taxid = reconstruct.initialize_target_taxids(args.target, args.targetfile,
                                                                        list_of_all_spectra, full_taxonomy)
            config_digestion = conf.config_digestion(config)
            set_of_markers, headers = pt.parse_peptide_tables(peptide_table, None, None)
            set_of_markers.update(compute_masses.add_deamidation(set_of_markers, set_of_codes_for_deamidation))
            markers.find_overlapping_markers(set_of_markers)
            set_of_all_new_markers = set()
            set_of_all_matching_markers = set()
            set_of_taxid = {m.taxid() for m in set_of_markers}
            config_substitution, config_gamma, config_conserved = conf.config_matrices_and_co(config)
            matrices = matrix.load_collagen_matrices(config_substitution)
            gamma_marker = gamma.initialize_gamma_from_csv_files(config_gamma, set_of_markers)
            dict_of_variable_positions=neighbour.compute_variable_positions(set_of_markers, gamma_marker)
            for target_taxid in list_of_target_taxid.keys():
                small_taxonomy, set_of_small_markers = taxo.update_taxonomy_and_set_of_markers(full_taxonomy,
                                                                                             set_of_markers,
                                                                                             full_taxonomy.descendants[
                                                                                                 target_taxid])
                set_of_close_species = taxo.close_species(small_taxonomy, target_taxid, set_of_taxid)
                set_of_CPs = conserved_markers.extract_CPs_from_file(target_taxid, config_conserved, full_taxonomy)
                list_of_spectra = list_of_target_taxid[target_taxid]
                for spectrum in list_of_spectra:
                    spectrum.add_median_intensity()
                config_selection = conf.config_selection_peaks(config)
                min_nb_spectra = len(list_of_spectra) * 0.1  # config_selection  # A reporter ailleurs
                # début des calculs
                set_of_close_markers = {m for m in set_of_markers if
                                           m.taxid() in set_of_close_species}  # reconstruct.build_candidate_markers(target_taxid, set_of_close_species, set_of_markers)

                # else:
                #   set_of_codes=lmt.deamidated_codes(list_of_constraints, set_of_selected_markers)
                # set_of_selected_markers.update(compute_masses.add_deamidation(set_of_selected_markers, set_of_codes))
                set_of_matching_markers, set_of_orphan_markers, list_of_orphan_spectra = reconstruct.find_orphan_peaks_and_peptides(
                    target_taxid, list_of_spectra, set_of_close_markers, error)
                # set_of_new_markers = set_of_matching_markers
                set_of_all_matching_markers.update(set_of_matching_markers)
                # set_of_new_markers=neighbour.find_candidate_sequences_for_orphan_markers(set_of_orphan_markers, target_taxid, list_of_orphan_spectra, error, min_nb_spectra, matrices, gamma_marker, set_of_small_markers, set_of_codes, config_digestion, set_of_CPs)
                set_of_neighbour_markers = neighbour.find_candidate_sequences_for_orphan_markers(set_of_orphan_markers,
                                                                                                 target_taxid,
                                                                                                 list_of_orphan_spectra,
                                                                                                 error, min_nb_spectra,
                                                                                                 matrices, gamma_marker,
                                                                                                 set_of_small_markers,
                                                                                                 dict_of_variable_positions,
                                                                                                 set_of_codes_for_deamidation,
                                                                                                 config_digestion,
                                                                                                 set_of_CPs)
                set_of_new_markers = pareto.add_aggregative_score(set_of_matching_markers | set_of_neighbour_markers)
                set_of_all_new_markers.update(set_of_new_markers)
            taxo.supplement_taxonomic_information(set_of_all_new_markers, full_taxonomy)
            # ta.supplement_taxonomic_information(set_of_all_matching_markers, full_taxonomy)
            pt.build_peptide_table_from_set_of_markers(set_of_all_new_markers, output, headers)
            rep.create_report_reconstruction(peptide_table, set_of_markers, target,  targetfile, spectra, list_of_target_taxid,list_of_all_spectra, taxonomy, full_taxonomy, web, set_of_close_species, config_substitution, config_gamma, config_conserved )
            deamidation = False

        if deamidation:
            set_of_markers, list_of_headers=pt.parse_peptide_tables(peptide_table, None, None)
            set_of_markers.update(compute_masses.add_deamidation(set_of_markers, set_of_codes_for_deamidation))
            pt.build_peptide_table_from_set_of_markers(set_of_markers,output, list_of_headers, config_markers)
            rep.create_report_deamidation(peptide_table, set_of_markers, set_of_codes_for_deamidation, web)

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
                print("Warnings were raised during execution.")
                print("Please refer to the warning.log file or the report file for details.\n\n")
        
    except message.InputError:
        rep.create_report_footer(output_dir, output, report)
        if not web:
           print("\n   An error occurred with your input. Stopping execution.")
           print(f"   Please refer to the warning.log file or the {report} file for more detail.\n\n")

if __name__ == "__main__":
    main()
