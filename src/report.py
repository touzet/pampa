import os
import time
import sys


# local import
from src import markers
from src import utils

def print_title(title):
    print("-"*30)
    print("   "+title)
    print("-"*30)
    print("")
    
def print_set_of_sequences(set_of_sequences):
    for seq in set_of_sequences:
        print ("  "+seq.seqid()+"\t"+seq.protein()+ "\t "+seq.taxon_name()+" (TaxID:"+seq.taxid()+")")
    print("")
    
def print_set_of_markers(set_of_markers):
    list_of_species=list({(m.taxon_name(),m.taxid())for m in set_of_markers})
    list_of_species.sort(key=lambda x: x[0])
    print("  Total number of species : "+str(len(list_of_species))+"\n")
    taxon_length=max({len(utils.pretty_print(s[0])+utils.pretty_print(s[1])) for s in list_of_species})
    taxid_length=max({len(utils.pretty_print(s[1])) for s in list_of_species})
    if taxon_length < 35:
        taxon_length=max({len(utils.pretty_print(s[0])) for s in list_of_species})
    else:
        taxon_length=35-taxid_length
    list_of_codes=list({str(m.code())+"-"+str(m.PTM()) for m in set_of_markers})
    list_of_codes.sort()
    matrix = [["" for j in range(len(list_of_codes)+2)] for i in range(len(list_of_species)+1)]
    for code in list_of_codes:
        matrix[0]=["",""]+list_of_codes
    for (i,sp) in enumerate(list_of_species):
        matrix[i+1][0]= list_of_species[i][0]
        matrix[i+1][1]= list_of_species[i][1]
    for m in set_of_markers:
        matrix[list_of_species.index((m.taxon_name(), m.taxid()))+1][list_of_codes.index(str(m.code())+"-"+str(m.PTM()))+2]=str(round(float(m.mass()),1)).rjust(6)
    for i in range(len(list_of_species)):
        species=utils.pretty_print(list_of_species[i][0])
        s= "  "+(species[:taxon_length] if len(species) > taxon_length else species)
        s= (s +" (taxID:"+utils.pretty_print(list_of_species[i][1])+")").ljust(taxon_length+taxid_length+10)
        for j in range(len(list_of_codes)):
            if len(matrix[i+1][j+2])>0 :
                s=s+matrix[i+1][j+2]+" "
        print(s)
    print("")

def print_peptide_table(peptide_table, web):
    if web:
        _, peptide_table_file = os.path.split(peptide_table)
        print("  Input peptide table: " + peptide_table_file)
    else:
        print("  Input peptide table: " + peptide_table)
    print("")
    
def print_peptide_tables(peptide_table, web):
    if peptide_table is None:
        return
    print ("  Petide table file(s) : ", end="")
    if web:
        for pep in peptide_table:
            _, pep_file = os.path.split(pep)
            print(pep_file, end=" ")
    else:
        for pep in peptide_table:
            print(pep, end=" ")
    print("")

def print_file(file_path, title, web):
    if file_path:
        if web:
            _, file_name = os.path.split(file_path)
            print("  "+title+" : "+file_name)
        else:
            print("  "+title+" : "+file_path)
      
def print_fasta(fasta, fasta_dir, web):
    if fasta:
        print_file(fasta, "Fasta file", web)
    else:
        if not web:
            print("  Fasta directory: "+fasta_dir)
        
            
def print_spectra(spectra_dir, list_of_spectra, web):
    if not web:
        print("  Directory: "+spectra_dir+"\n")
    print("  "+str(len(list_of_spectra))+" spectral files found\n")
    for f in list_of_spectra:
        print("  " + f.name + " (" + str(len(f)) + " peaks)")
    print("")
        
def print_error(error):
    if error is None:
        return
    print("  Error margin tolerance : "+str(error), end=" ")
    if error<1:
        print("Da")
    else:
        print("ppm")
            
def print_limit(list_of_constraints):
    for d in list_of_constraints:
        print("  ",end="")
        for key in d:
            print (utils.restitute_field(key)+" : "+str(d[key]), end=" ")
    print("\n")
        
def print_digestion(config_digestion):
    print("  In silico digestion :")
    print("     - Enzyme: "+ config_digestion["enzyme"])
    print("     - Maximal number of missed cleavages : "+str(config_digestion["number_of_missed_cleavages"]))
    print("     - Minimal peptide length : "+str(config_digestion["min_peptide_length"]))
    print("     - Maximal peptide length : "+str(config_digestion["max_peptide_length"]))
    
def create_report_classify(spectra_dir, list_of_spectra, taxonomy, taxonomy_tree, peptide_table, fasta, fasta_dir, set_of_sequences, set_of_markers, limit, list_of_constraints, deamidation, error, neighbour, all, new_table, config_digestion, config_nb_of_peaks, web):
    # TO DO: display constraints
    print ("PAMPA CLASSIFY\n")
    print_title("MASS SPECTRA")
    print_spectra(spectra_dir, list_of_spectra, web)
    print_title("PEPTIDE MARKERS")
    if peptide_table:
        print_peptide_tables(peptide_table, web)
        print_set_of_markers(set_of_markers)
        markers.check_set_of_markers(set_of_markers)
    else:
        print("  Markers automatically infered from sequences in")
        print_fasta(fasta, fasta_dir, web)
        print_set_of_sequences(set_of_sequences)
        print_digestion(config_digestion)
        print("")
        print("  The corresponding peptide table is in "+new_table, end="\n\n")
    if limit:
        print_title("LIMITS")
        print_file(limit, 'Limit', web)
        print_limit(list_of_constraints)
    print_title("PARAMETERS")
    print("  Minimum number of peaks : " + str(config_nb_of_peaks))
    print("  Near-optimal solutions  : ",end="")
    if neighbour==100:
        print("only solution with the highest number of matching peaks")
    else:
        print("up to "+str(neighbour)+"% matching peaks")
    print("  Selection of solutions  : ", end="")
    if not all:
        print ("peak intensity and inclusion selection")
    else:
        print ("all possible solutions")
    print_error(error)
    if deamidation:
        print ("  Deamidation            : Yes")
    elif peptide_table:
        print ("  Deamidation            : Only those present in the peptide table")
    else:
        print ("  Deamidation            : None")
    #markers.colinearity(set_of_markers)
    #markers.check_set_of_markers(set_of_markers)
    if taxonomy :
        print_title("TAXONOMY")
        print_file(taxonomy, 'Taxonomy', web)
        #ta.table_print(taxonomy_tree)

def create_report_homology(peptide_table, set_of_markers, fasta, fasta_dir, set_of_sequences, taxonomy,  config_digestion, limit, list_of_constraints, web):
    print("PAMPA CRAFT, mode HOMOLOGY\n")
    print_title("INPUT FILES")
    print_fasta(fasta, fasta_dir, web)
    print_peptide_tables(peptide_table, web)
    print_file(limit, 'Limit', web)
    print_file(taxonomy, 'Taxonomy', web)
    print("")
    print_title("FASTA SEQUENCES")
    print_set_of_sequences(set_of_sequences)
    print_title("INPUT PEPTIDE TABLE")
    print_set_of_markers(set_of_markers)
    if limit:
        print_title("LIMITS")
        print_limit(list_of_constraints)
    print_title("PARAMETERS")
    print_digestion(config_digestion)
    # markers.check_set_of_markers(set_of_markers)
    print("")
    
def create_report_selection(spectra_dir, list_of_spectra, peptide_table, set_of_markers, config_selection, error, web):
    print("PAMPA CRAFT, mode SELECTION")
    print_title("MASS SPECTRA")
    print_spectra(spectra_dir, list_of_spectra, web)
    print_title("PEPTIDE MARKERS")
    print_peptide_tables(peptide_table, web)
    print_set_of_markers(set_of_markers)
    print_title("PARAMETERS")
    print("  Minimum proportion of spectra : " + str(config_selection))
    print_error(error)

    
def create_report_deamidation(peptide_table, set_of_markers, set_of_codes, web):
    print("PAMPA CRAFT, mode DEAMIDATION\n")
    if len(set_of_codes)==0:
        print("  Modified peptide markers: all\n")
    else:
        print("  Modified peptide markers: "+ str(set_of_codes)+"\n")
    print("PEPTIDE MARKERS")
    print_peptide_tables(peptide_table, web)
    print_set_of_markers(set_of_markers)
   
def create_report_allpeptides(fasta, fasta_dir, set_of_sequences, config_digestion, limit, list_of_constraints, web, spectra_dir=None, list_of_spectra=None, error=None, config_selection=None):
    print("PAMPA CRAFT, mode ALL PEPTIDES")
    print_title("INPUT SEQUENCES")
    print_fasta(fasta, fasta_dir, web)
    print_set_of_sequences(set_of_sequences)
    if limit:
        print_title("LIMITS")
        print_file(limit, 'Limit', web)
        print_limit(list_of_constraints)
    if spectra_dir:
        print_title("MASS SPECTRA")
        print_spectra(spectra_dir, list_of_spectra, web)
    print_title("PARAMETERS")
    print_digestion(config_digestion)
    if spectra_dir:
        print("  Minimum proportion of spectra  :" + str(config_selection))
    print_error(error)
        
    
def create_report_supplement(peptide_table, set_of_markers, web, taxonomy, list_of_markers=None, set_of_sequences=None):
    print("PAMPA CRAFT, mode SUPPLEMENT \n")
    print_title("INPUT PEPTIDE TABLE")
    print_peptide_tables(peptide_table, web)
    if taxonomy:
        print_title("TAXONOMY")
        print_file(taxonomy, 'Taxonomy', web)
    if set_of_sequences:
        print_title("INPUT SEQUENCES")
        print_set_of_sequences(set_of_sequences)
        print_title("NEW PEPTIDE TABLE")
        print_set_of_markers(set_of_markers)
        #markers.check_set_of_markers(set(list_of_markers))

def create_report_header(command_line, report):
    sys.stdout=open(report, 'w')
    print("=====================================================================\n")
    print("                              P A M P A                              \n")
    print("=====================================================================\n")
    print (time.ctime())
    print("")
    print(command_line)
    print("")
    

def create_report_footer(output_dir, output, report):
    if os.path.getsize(os.path.join(output_dir,'warning.log')) > 0 and os.path.getsize(os.path.join(output_dir,'error.log'))==0:
        print_title("WARNINGS")
        with open(os.path.join(output_dir,'warning.log'), 'r') as file:
            for line in file:
                print("  "+line, end="")
        print("")
    if os.path.getsize(os.path.join(output_dir,'error.log')) > 0:
        print("\n* * * * *    FATAL ERROR    * * * * *")
        with open(os.path.join(output_dir,'warning.log'), 'r') as file:
            for line in file:
                print("  "+line, end="")
        print("\n* * * * *   NO OUTPUT FILE  * * * * *")
    else:
        print_title("OUTPUT FILES")
        print("  Main result file (TSV) : "+output)
        print("  Report (this file)     : "+report)
    print("")
    sys.stdout = sys.__stdout__
    
