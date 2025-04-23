import os

from src import message
from src import report as rep

def logger_and_outputdir_configuration(output):
    if output is None:
        rep.create_report_header("report.txt")
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
    output2=os.path.join(output_dir, "detail_"+output_file+".tsv")  ## TO DO
    return output_dir, output_file, report_file

def check_config(config):
    if config is None:
        return "config.json"
    if not os.path.isfile(config):
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
        if not os.path.isdir(directory):
            message.escape("Directory "+directory+" not found (-d). Stopping execution")
            
def check_spectra(spectra, error, mandatory=True):
    if spectra is None :
        if mandatory:
            message.escape("Missing parameter: spectra (-s). Stopping execution")
        else:
            return
    if not os.path.isdir(spectra):
        message.escape("Directory "+spectra+" not found. Stopping execution")
    if error is None:
        message.escape("Missing parameter: error (-e). Stopping execution")
    if error<0:
        message.escape("Parameter error (-e) should be positive. Stopping execution")

def check_and_update_parameters_classify(spectra, taxonomy, peptide_table, fasta, directory, limit, deamidation, error, neighbour, all, output, mammals):
    """
    Parameters checking and fixing. Configuration of loggers
    """

    output,report, report_file, output2=check_output(output_dir, output_file)
    config=check_config(config)
    check_limit(limit)
    check_spectra(spectra, error)

    if mammals and (peptide_table or fasta or directory):
        message.warning("Parameters -p, -f and -d are not compatible with --mammals. Applying --mammals parameter.")
        peptide_table=None
        fasta=None
        directory=None

    if mammals :
        if not os.path.isfile("Taxonomy/taxonomy_mammals.tsv"):
            message.escape("The taxonomy file has been deleted.")
        if not os.path.isfile("Peptide_tables/table_mammals.tsv"):
            message.escape("The peptide table file has been deleted.")
        taxonomy="Taxonomy/taxonomy_mammals.tsv"
        peptide_table=["Peptide_tables/table_mammals.tsv"]

    check_taxonomy(taxonomy)
    
    if neighbour not in range(101):
        neighbour=100
        message.warning("Parameter -n (neighbouring): value is 100")
    
    if peptide_table:
        if  fasta or directory:
            message.escape("Options -p (peptide_table), -f (fasta) and -d (directory of fasta files) are mutually incompatible. Stopping execution")
        else:
            check_peptide_table(peptide_table)
            new_table=None
    elif fasta or directory:
        check_sequences(fasta, directory)
        new_table=os.path.join(output_dir, "table_"+output_file)
    else:
        message.escape("Missing information for marker peptides (-p, -f or -d). Stopping execution")
    
    return (spectra, taxonomy, peptide_table, fasta, directory, limit, deamidation, error, neighbour, all, output, output2, output_dir, report_file, new_table)



def check_and_update_parameters_craft(homology, deamidation, allpeptides, fillin, selection, peptide_table, fasta, directory, spectra, resolution, limit, taxonomy, config):
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
        check_sequences(fasta, directory)
        check_taxonomy(taxonomy)
        
    if deamidation:
        check_peptide_table(peptide_table)
        
    if allpeptides:
        check_sequences(fasta, directory)
        check_spectra(spectra, resolution, False)
    
    if fillin:
        check_peptide_table(peptide_table)
        check_sequences(fasta, fasta_dir, False)
        check_taxonomy(taxonomy)
            
    if selection:
        check_peptide_table(peptide_table)
        check_spectra(spectra, error)
            
    return (homology, deamidation, allpeptides, fillin, selection, peptide_table, fasta, directory, spectra, resolution, limit, taxonomy, config)

