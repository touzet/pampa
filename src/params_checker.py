import os

from src import message
from src import report as rep

def logger_and_outputdir_configuration(output, command_line):
    if output is None:
        rep.create_report_header(command_line, "report.txt")
        message.configure("")
        message.escape("Missing parameter: output (-o).")
    output_dir, output_file = os.path.split(output)
    if len(output_dir)>0 :
        # Ensure the output directory exists. If not, create it.
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    message.configure(output_dir)
    extension=output_file[-4:].lower()
    if extension!=".tsv":
        output_file=output_file+".tsv"
    else:
        output_file=output_file[:-4]+".tsv"
    report_file="report_"+output_file.replace("tsv", "txt")
    output_detail="detail_"+output_file[:-4]+".tsv"
    output_json=output_file[:-4]+".json"
    return output_dir, output_file, report_file, output_detail, output_json

def check_config(config):
    if config is None:
        if not os.path.isfile("config.json"):
            message.escape("File config.json not found. Stopping execution.")
        else:
            return "config.json"
    if not os.path.isfile(config):
        if not os.path.isfile("config.json"):
            message.warning("User file "+config+" not found (-c).")
            message.escape("Default file config.json not found. Stopping execution.")
        else:
            message.warning("File "+config+" not found (-c). Using config.json instead.")
            return "config.json"
    return config
    
def  check_peptide_table(peptide_table):
    if peptide_table is None:
        message.escape("Missing parameter: -p (peptide table)")
    for pep in peptide_table:
        if not os.path.isfile(pep):
            message.escape("File "+pep+" not found (-p).")
            
def check_limit(limit):
    if limit:
        if not os.path.isfile(limit):
            message.warning("File "+limit+" not found. No limit applied.")
        elif os.path.getsize(limit) == 0:
            message.warning("File "+limit+" is empty. No limit applied.")
            
def check_taxonomy(taxonomy):
    if taxonomy:
        if not os.path.isfile(taxonomy):
            message.warning("File "+taxonomy+" not found. Ignored.")
            taxonomy=None
        elif  os.path.getsize(taxonomy) == 0:
            message.warning("File "+taxonomy+" is empty. Ignored.")
            taxonomy=None

def check_sequences(fasta, fasta_dir, mandatory=True):
    if fasta and fasta_dir:
        message.escape("Options -f (fasta file) and -d (directory of fasta files) are mutually exclusive. Stopping execution")
    if not (fasta or fasta_dir) and mandatory:
        message.escape("Missing target sequences (-f or -d). Stopping execution")
    if fasta:
        if not os.path.isfile(fasta):
            message.escape("File "+fasta+" not found (-f). Stopping execution")
        if os.path.getsize(fasta) == 0:
            message.escape("File "+fasta+" is empty. Stopping execution")
    if fasta_dir:
        if not os.path.isdir(fasta_dir):
            message.escape("Directory "+fasta_dir+" not found (-d). Stopping execution")
            
def check_spectra(spectra, mandatory=True):
    if spectra is None :
        if mandatory:
            message.escape("Missing parameter: spectra (-s). Stopping execution")
        else:
            return
    if not os.path.isdir(spectra):
        message.escape("Directory "+spectra+" not found. Stopping execution")
    
def check_error(error, mandatory=True):
    if error is None:
        if mandatory:
            message.escape("Missing parameter: error (-e). Stopping execution")
        else:
            return
    if error<0:
        message.escape("Parameter error (-e) should be a positive value. Stopping execution")
        
def check_spectra_and_error(spectra, error, mandatory=True):
    mandatory= mandatory or spectra or error
    check_spectra(spectra, mandatory)
    check_error(error, mandatory)
        
def useless_parameters(list_of_parameters):
    for p in list_of_parameters:
        if p[0] is not None:
            message.warning("Useless parameter: "+p[1]+" "+str(p[0])+". Ignored.")
            

def check_and_update_parameters_classify(spectra, taxonomy, peptide_table, fasta, fasta_dir, limit, deamidation, error, neighbour, all, mammals, config):
    """
    Parameters checking and fixing. Configuration of loggers
    """
    config=check_config(config)
    check_limit(limit)
    check_spectra_and_error(spectra, error)

    if mammals and (peptide_table or fasta or fasta_dir):
        message.warning("Parameters -p, -f and -d are not compatible with --mammals. Applying --mammals parameter.")
        peptide_table=None
        fasta=None
        fasta_dir=None

    if mammals :
        if not os.path.isfile("Taxonomy/taxonomy_mammals.tsv"):
            message.escape("The file Taxonomy/taxonomy_mammals.tsv is not found.")
        if not os.path.isfile("Peptide_tables/table_mammals.tsv"):
            message.escape("Peptide_tables/table_mammals.tsv is missing.")
        taxonomy="Taxonomy/taxonomy_mammals.tsv"
        peptide_table=["Peptide_tables/table_mammals.tsv"]

    check_taxonomy(taxonomy)
    
    if neighbour not in range(101):
        neighbour=100
        message.warning("Parameter -n (neighbouring): value is 100")
    
    if peptide_table:
        if  fasta or fasta_dir:
            message.escape("Options -p (peptide_table), -f (fasta) and -d (directory of fasta files) are mutually incompatible. Stopping execution")
        else:
            check_peptide_table(peptide_table)
            new_table=None
    elif fasta or fasta_dir:
        check_sequences(fasta, fasta_dir)
    else:
        message.escape("Missing information for marker peptides (-p, -f or -d). Stopping execution")
    
    return (spectra, taxonomy, peptide_table, fasta, fasta_dir, limit, deamidation, error, neighbour, all, config)


def check_and_update_parameters_craft(homology, deamidation, allpeptides, fillin, selection, peptide_table, fasta, fasta_dir, spectra, resolution, limit, taxonomy, config):
    """
    Parameters checking and fixing for PAMPA CRAFT.
    Configuration of loggers
    """
    
    param=sum([homology, deamidation, allpeptides, fillin, selection])
    if param==0:
         message.escape("Missing parameter: --homology, --allpeptides, --fillin, --deamidation, or --selection. Stopping execution")
    if param>1:
        message.escape("Parameters --homology, --allpeptides, --deamidation, --selection and --fillin are mutually exclusive.")
    
    config=check_config(config)
    check_limit(limit)

    if homology :
        check_peptide_table(peptide_table)
        check_sequences(fasta, fasta_dir)
        check_taxonomy(taxonomy)
        useless_parameters([(spectra, '-s'), (resolution,'-e')])
        
    if deamidation:
        check_peptide_table(peptide_table)
        useless_parameters([(fasta, '-f'), (fasta_dir,'-d'), (taxonomy, '-t'), (spectra,'-s'), (resolution,'-e')])

    if allpeptides:
        check_sequences(fasta, fasta_dir)
        check_spectra_and_error(spectra, resolution, False)
        useless_parameters([(peptide_table,'-p')])
    
    if fillin:
        check_peptide_table(peptide_table)
        check_sequences(fasta, fasta_dir, False)
        check_taxonomy(taxonomy)
        check_error(resolution, False)
        useless_parameters([(spectra,'-s')])
         
    if selection:
        check_peptide_table(peptide_table)
        check_spectra_and_error(spectra, resolution)
        useless_parameters([(fasta,'-f'), (fasta_dir,'-d'), (taxonomy,'-t')])
            
    return (homology, deamidation, allpeptides, fillin, selection, peptide_table, fasta, fasta_dir, spectra, resolution, limit, taxonomy, config)

