#!/usr/bin/env python3

import argparse
import sys
import time
import os

# local import

from src import classify


def main():
    parser = argparse.ArgumentParser(#formatter_class=CustomFormatter,
                                         usage="pampa_classify [-h]\n   -s SPECTRA \n   -e ERROR \n   -o OUTPUT FILE \n   -p PEPTIDE_TABLE [PEPTIDE_TABLE] | -f FASTA  | -d DIRECTORY\n   [-l LIMIT]\n   [-t TAXONOMY]\n   [-n NEIGHBOURING] [-a]", description="This script is for species identification.")


    group2=parser.add_argument_group('\nGeneral options')
    group2.add_argument("-s", dest="spectra", help="Path to the spectra files (one spectrum per file). Authorized formats: cvs, mgd, mzML.", type=str, required=True)
    group2.add_argument("-e", dest="error", help="Maximal error margin for the observation (in Dalton or ppm). Recommended values: 0.02 for maldi FT and 0.1 for maldi TOF.", type=float) 
    group2.add_argument("-o", dest="output", help="Output path (should include the output file name)", type=str)
    group3=parser.add_argument_group('\nOptions for organism selection')
    group3.add_argument("-p", dest="peptide_table",nargs='+', help="Peptide table(s) (TSV file(s))", type=str)
    group3.add_argument("-f", dest="fasta", help="Fasta sequences", type=str)
    group3.add_argument("-d", dest="directory",  help="Directory where to find  Fasta files.", type=str)
    group3.add_argument("-l", dest="limit",  help="Limit the set of peptides or fasta sequences to organisms, molecules or sequence ID specified in this file  (TSV files)", type=str, required=False)
    group3.add_argument("-t", dest="taxonomy", help="Taxonomy (TSV file), optional.", type=str)
    group4=parser.add_argument_group('\nOptions for suboptimal solutions')
    group4.add_argument("-n", dest="neighbour", help="Provide near-optimal solutions within a specified percentage margin, ranging between 0 and 100. Default is 100. With -n 100, only optimal solutions are provided.", type=int, required=False, default=100)
    group4.add_argument("-a", dest="all", action='store_true', help="Provide all solutions within the percentage margin specified with option  -n, and not only suboptimal solutions.  Default is False.", required=False)  
    group4.add_argument("--web", dest="web",  action='store_true', help=argparse.SUPPRESS, required=False)
    
    args = parser.parse_args()
    
    classify.main(args.spectra, args.taxonomy, args.peptide_table, args.fasta, args.directory, args.limit, args.error, args.neighbour, args.all, args.output,None, False, args.web)
 
if __name__ == "__main__":
    main()
