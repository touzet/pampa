import os
import time
import sys

from src import taxonomy
# local import
from src import markers
from src import utils

def print_title(title):
    print("-"*30)
    print("   "+title)
    print("-"*30)
    print("")
    
def print_set_of_sequences(set_of_sequences):
    print("  "+str(len(set_of_sequences))+" FASTA sequences selected\n")
    seqid_length=max({len(utils.pretty_print(seq.seqid())) for seq in set_of_sequences})
    for seq in set_of_sequences:
        if seq.taxid() is None:
            print("  " + utils.pretty_print(seq.seqid()).ljust(seqid_length) + "\t" + utils.pretty_print(
                seq.protein()) + "\t " + utils.pretty_print(seq.taxon_name()))
        else:
            print ("  "+utils.pretty_print(seq.seqid()).ljust(seqid_length)+"\t"+utils.pretty_print(seq.protein())+ "\t "+utils.pretty_print(seq.taxon_name())+" (TaxID:"+utils.pretty_print(seq.taxid())+")")
    print("")

def print_set_of_markers_light(set_of_markers):
    list_of_species=sorted(list({m.taxon_name() for m in set_of_markers}))
    set_of_codes={m.code() for m in set_of_markers}
    print("")
    print("Number of clades: ", len(list_of_species))
    print(*list_of_species, sep=", ")

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
        if m.mass() is not None:
            matrix[list_of_species.index((m.taxon_name(), m.taxid()))+1][list_of_codes.index(str(m.code())+"-"+str(m.PTM()))+2]=str(round(float(m.mass()),1)).rjust(6)
        else:
            matrix[list_of_species.index((m.taxon_name(), m.taxid()))+1][list_of_codes.index(str(m.code())+"-"+str(m.PTM()))+2]=""
    for i in range(len(list_of_species)):
        species=utils.pretty_print(list_of_species[i][0])
        s= (species[:taxon_length] if len(species) > taxon_length else species)
        if list_of_species[i][1] is not None:
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
    if len(peptide_table)==1:
        print ("  Peptide table file : ", end="")
    else:
        print ("  Peptide table files : ", end="")
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
    print("")

def print_fasta(fasta, fasta_dir, web):
    if fasta:
        print_file(fasta, "Fasta file", web)
    else:
        if not web:
            print("  Fasta directory: "+fasta_dir)

def print_mature(set_of_sequences):
    processed=""
    for seq in set_of_sequences:
        if "Mature" in seq.field and seq.field["Mature"] is not None:
            min,max=seq.field["Mature"]
            processed += " "+seq.seqid()+"["+str(min)+","+str(max)+"]"
    print("  Sequences edited to extract the helical region:")
    if len(processed) > 0:
        print(processed)
    else:
        print("  None")
    print("")

def print_spectra(spectra_dir, list_of_spectra, web):
    if not web:
        print("  Directory: "+spectra_dir+"\n")
    print("  "+str(len(list_of_spectra))+" spectral files found\n")
    for f in list_of_spectra:
        print("  " + f.name + " (" + str(len(f)) + " peaks)")
    print("")

def  print_spectra_taxid(spectra_dir, list_of_spectra, list_of_target_taxid, web):
    if not web:
        print("  Directory: " + spectra_dir + "\n")
    for taxid in list_of_target_taxid:
        print("  " + taxid+" -> ", end="")
        for f in list_of_target_taxid[taxid]:
            print("  " + f.name + " (" + str(len(f)) + " peaks) ")
print("")

def print_error(error):
    if error is None:
        return
    print("  Error margin tolerance  : "+str(error), end=" ")
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

def print_markers_with_hidden_numbers(set_of_markers, peptide_table, web):
    print_title("INITIAL PEPTIDE MARKERS")
    print_peptide_tables(peptide_table, web)
    print("")
    for m in sorted(set_of_markers, key=lambda m: (len(m.field["magic_number"]), m.field["magic_number"])):
        print("  "+m.field["magic_number"].ljust(6)+str(m))
    print("")

def print_matrices(config_substitution, config_gamma, config_conserved):
    print("  Substitution matrices :", config_substitution)
    print("  Gamma matrices :", config_gamma)
    print("  Conservation matrices :", config_conserved)
    print("")
    print("  These matrices were generated from your FASTA sequences. \n  They are available in the main directory as CSV or JSON files.")
    print("  matrix_custom: one matrix per protein family and per X and Y positions (refer to the G-X-Y repeats in collagens).")
    print("  gamma_custom: one matrix per protein family.")
    print("  conserved_custom: one matrix for all sequences.")
    print("")

def print_digestion(config_digestion):
    print("  In silico digestion :")
    print("     - Enzyme: "+ config_digestion["enzyme"])
    print("     - Maximal number of missed cleavages : "+str(config_digestion["number_of_missed_cleavages"]))
    print("     - Minimal peptide length : "+str(config_digestion["min_peptide_length"]))
    print("     - Maximal peptide length : "+str(config_digestion["max_peptide_length"]))
    print("\n")

def print_deamidation(deamidation, pep_table=False):
    if deamidation:
        print("  Deamidation : Yes\n")
    elif pep_table:
        print("  Deamidation : Only those present in the input peptide table\n")
    else:
        print("  Deamidation : None \n")

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
        print("  Markers automatically inferred from sequences in")
        print_fasta(fasta, fasta_dir, web)
        print("")
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
        print("Only solutions with the highest number of matching peaks")
    else:
        print("up to "+str(neighbour)+"% matching peaks")
    print("  Selection of solutions  : ", end="")
    if not all:
        print ("Peak intensity and inclusion selection")
    else:
        print ("No selection on peak intensity")
    print_error(error)
    print_deamidation(deamidation, peptide_table)
    print("")
    if taxonomy :
        print_title("TAXONOMY")
        print_file(taxonomy, 'Taxonomy', web)
        #ta.table_print(taxonomy_tree)
        print("")

def create_report_homology(peptide_table, set_of_markers, list_of_new_markers, fasta, fasta_dir, set_of_sequences, taxonomy,  config_digestion, limit, list_of_constraints, deamidation, web):
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
    print_deamidation(deamidation, True)
    print("")
    # markers.check_set_of_markers(set_of_markers)
    print_title("OUTPUT PEPTIDE TABLE")
    if len(list_of_new_markers)>0:
        print("   Number of found peptides:" , len(list_of_new_markers))
    else:
        print("   No peptide marker found.")
    set_of_homo_species = {m.taxon_name() for m in list_of_new_markers if m.taxon_name() is not None}
    set_of_initial_species = {m.taxon_name() for m in set_of_markers if m.taxon_name() is not None}
    list_of_new_species= sorted(set_of_homo_species - set_of_initial_species)
    if len(list_of_new_species)>0:
        print("   - Total number of species: ", len(set_of_homo_species))
        print("   - New species: ", end="")
        print(", ".join(list_of_new_species))
    print("\n")
    
def create_report_selection(spectra_dir, list_of_spectra, peptide_table, set_of_markers, config_selection, error, web):
    print("PAMPA CRAFT, mode SELECTION")
    print_title("MASS SPECTRA")
    print_spectra(spectra_dir, list_of_spectra, web)
    print_title("INPUT PEPTIDE TABLE")
    print_peptide_tables(peptide_table, web)
    print_set_of_markers(set_of_markers)
    print_title("PARAMETERS")
    print("  Minimum proportion of spectra : " + str(config_selection))
    print_error(error)

    
def create_report_deamidation(peptide_table, set_of_markers, set_of_codes, web):
    print("PAMPA CRAFT, mode DEAMIDATION\n")
    print_title("INPUT PEPTIDE TABLE")
    print_peptide_tables(peptide_table, web)
    print_set_of_markers(set_of_markers)
    print_title("DEAMIDATION")
    if set_of_codes is None:
        print("  Modified peptide markers: all\n")
    else:
        print("  Modified peptide markers: "+ str(set_of_codes)+"\n")
   
def create_report_allpeptides(fasta, fasta_dir, set_of_sequences, config_digestion, limit, list_of_constraints, deamidation, web,  spectra_dir=None, list_of_spectra=None, error=None, config_selection=None):
    print("PAMPA CRAFT, mode ALL PEPTIDES\n")
    if limit:
        print_title("LIMITS")
        print_file(limit, 'Limit', web)
        print_limit(list_of_constraints)
    print_title("INPUT SEQUENCES")
    print_fasta(fasta, fasta_dir, web)
    print_set_of_sequences(set_of_sequences)
    if spectra_dir:
        print_title("MASS SPECTRA")
        print_spectra(spectra_dir, list_of_spectra, web)
    print_title("PARAMETERS")
    print_digestion(config_digestion)
    print_deamidation(deamidation)
    if spectra_dir:
        print("  Minimum proportion of spectra  :" + str(config_selection))
    print_error(error)
        
    
def create_report_supplement(peptide_table, fasta, fasta_dir,set_of_markers, web, taxonomy, set_of_sequences=None):
    #print_markers_with_hidden_numbers(set_of_markers, peptide_table, web)
    if taxonomy:
        print_title("TAXONOMY")
        print_file(taxonomy, 'Taxonomy', web)
        print("")
    if set_of_sequences:
        print_title("FASTA SEQUENCES")
        print_fasta(fasta, fasta_dir, web)
        print_set_of_sequences(set_of_sequences)
        #print_title("NEW PEPTIDE TABLE")
        #print_set_of_markers(set_of_markers)
        #markers.check_set_of_markers(set(list_of_markers))

def create_report_reconstruction(peptide_table, set_of_markers, target,  targetfile, spectra_dir, list_of_target_taxid, list_of_spectra, taxonomy_file, taxonomy_tree, web, set_of_close_species, config_substitution, config_gamma, config_conserved, set_of_sequences=None ):
    #target = ta.find_all_taxonomic_information_from_taxid(target_taxid, final_taxonomy)
    print("PAMPA CRAFT, mode RECONSTRUCTION \n")
    print_title("INPUT FILES")
    print_title("TARGET TAXONS and MASS SPECTRA")
    print_spectra_taxid(spectra_dir, list_of_spectra, list_of_target_taxid, web)
    print_title("TAXONOMIC MODEL")
    print_peptide_tables(peptide_table, web)
    print_file(taxonomy_file, 'Taxonomy', web)
    taxonomy.table_print(taxonomy_tree)
    print_set_of_markers_light(set_of_markers)
    print("")

    print("")
    if set_of_sequences:
        print_title("FASTA SEQUENCES")
        print_set_of_sequences(set_of_sequences)
        print_mature(set_of_sequences)
        print_title("PARAMETERS for PEPTIDE GENERATION")
        print_matrices("matrix_custom", "gamma_custom", "conserved_custom")
    if config_substitution and config_gamma and config_conserved:
        print_title("PARAMETERS for PEPTIDE GENERATION")
        print_matrices(config_substitution, config_gamma, config_conserved)
        print("  These matrices are available in the 'Matrices' subdirectory as CSV files.\n")
    print_title("PROCESSING")
    print("  Set of close species: \n")
    for taxid in set_of_close_species:
        print ("  ", taxonomy_tree.name[taxid], " ["+taxid+"]")
    print("")

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
        print("\n* * * * *    FATAL ERROR    * * * * *\n")
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
    
