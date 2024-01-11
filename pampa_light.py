"""
pampa_ligth.py

light version of pampa:
assignment with preset options, coming with predefined peptide tables and predefined taxonomies

"""

#!/usr/bin/env python3


import argparse
from src import main_assign


def main():
    parser = argparse.ArgumentParser(usage="\n \n pampa_light [-h] <general options> <preset options for organism selection> <non-mandatory advanced options>")
    group2=parser.add_argument_group('\ngeneral options')
    group2.add_argument("-s", dest="spectra", help="Path to the spectra files (one spectrum per file). Authorized formats: cvs, mgd, mzML.", type=str, required=True)
    group2.add_argument("-e", dest="error", help="Maximal error margin for the observation (in Dalton or ppm). Recommended values: 0.02 for maldi FT and 0.1 for maldi TOF.", type=float, required=True) 
    group2.add_argument("-o", dest="output", help="Output path (should include the output file name)", type=str, required=True)

    group = parser.add_argument_group('preset options for organism selection')
    group.add_argument("--mammals", help="Preset option: use pre-computed peptide table and taxonomy for mammals.", action='store_true')

    group3=parser.add_argument_group('\nnon-mandatory advanced options: organism and markers selection')
    group3.add_argument("-p", dest="peptide_table", help="Peptide table (TSV file)", type=str)
    group3.add_argument("-f", dest="fasta", help="Fasta sequences", type=str)
    group3.add_argument("-d", dest="directory",  help="Directory where to find  Fasta files.", type=str)
    group3.add_argument("-l", dest="limit",  help="Limit the set of peptides or fasta sequences to organisms, molecules or sequence ID specified in this file  (TSV files)", type=str, required=False)
    group3.add_argument("-t", dest="taxonomy", help="Taxonomy (TSV file), optional.", type=str)

    group4=parser.add_argument_group('\nnon-mandatory advanced options: suboptimal solutions')
    group4.add_argument("-n", dest="neighbour", help="Also provide near-optimal solutions within a specified percentage margin, ranging between 0 and 100. Default is 100. Only suboptimal solutions are provided.", type=int, required=False, default=100)
    group4.add_argument("-a", dest="all", action='store_true', help="Provide all solutions, and not only suboptimal solutions, within the percentage margin specified with parameter -n. Default is False.", required=False)  
  
    args = parser.parse_args()

    main_assign.main(args.spectra, args.taxonomy, args.peptide_table, args.fasta, args.directory, args.limit, args.error, args.neighbour, args.all, args.output, args.mammals, True)  
    
if __name__ == "__main__":
    main()
