import os

from src import message
from src import report as rep
from src import config as conf

def logger_and_outputdir_configuration(output, command_line):
    if output is None:
        rep.create_report_header(command_line, "report.txt",True)
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

# output: name of the config_file associatd to the taxonomic model
def taxonomic_model(mammals, placentals, birds):
    param = sum([placentals, birds, mammals])
    if param == 0:
        return None
    if param > 1:
        message.escape(
            "Taxonomic models placentals, birds and mammals are mutually exclusive. Stopping execution.")
    if placentals:
        model = "placentals"
    elif birds:
        model = "birds"
    else : #mammals:
        model = "mammals"
    config = "config_" + model + ".json"
    return  check_file_warning(config)

def check_config(config, config_taxo=None):
    if config is not None:
        if not os.path.isfile(config):
            message.escape(f"File {config} not found. Stopping execution.")
        else:
            return config
    if config_taxo is None:
        if not os.path.isfile("config.json"):
            message.escape("File config.json not found. Stopping execution.")
        else:
            return "config.json"
    else:
        return config_taxo

def check_taxonomic_model(config, peptide_table, taxo):
    new_peptide_table = None
    new_taxo = None
    if config:
        new_peptide_table=conf.config_peptide_table(config)
        new_taxo=conf.config_taxonomy(config)
        if new_peptide_table is None or new_taxo is None:
            message.warning("File "+config+" not found. Stopping execution.")
    if new_taxo and taxo:
        message.warning(f"Taxonomic model: parameter -t ({taxo}) ignored.")
    if new_taxo:
        taxo = new_taxo
    if new_peptide_table and peptide_table:
        message.warning(f"Taxonomic model: parameter -p ({peptide_table}) ignored.")
    if new_peptide_table:
        peptide_table = new_peptide_table
    return peptide_table, taxo




def check_model_classify(config, peptide_table, taxo, fasta, fasta_dir):
    if fasta and fasta_dir:
        message.escape(
            f"Parameters -f {fasta} and -d {fasta_dir} are mutually exclusive. Stopping execution.")
    peptide_table, taxo = check_taxonomic_model(config, peptide_table, taxo)
    if not  (fasta or fasta_dir or peptide_table):
        message.escape("No taxonomic model provided. Stopping execution.")
    if not (fasta or fasta_dir):
        return peptide_table, taxo, fasta, fasta_dir
    if (fasta or fasta_dir) and peptide_table:
        message.escape(
            f"Parameters -f {fasta}, -d {fasta_dir} are incompatible with any other taxonomic model (-p, --mammals, --placentals or --birds). Ignored." )
        fasta = None
        fasta_dir = None
    return peptide_table, taxo, fasta, fasta_dir



def  check_peptide_table(peptide_table):
    if peptide_table is None:
        message.escape("Missing parameter: -p (peptide table)")
    for pep in peptide_table:
        if not os.path.isfile(pep):
            message.escape(f"File {pep} not found (-p).")

def check_matrices(config_matrices):
    # TO DO
    """"   
    for gene in list_of_genes:
        path="Matrices/matrix_"+model+"_"+gene+"_X.csv"
        if not os.path.isfile(path):
            message.escape("File "+path+ " is missing.")
        path = "Matrices/matrix_" + model + "_" + gene + "_Y.csv"
        if not os.path.isfile(path):
            message.escape("File " + path + " is missing.")
    """

def check_gamma(config_gamma):
    for gamma_file in config_gamma:
        if not os.path.isfile(gamma_file):
            message.escape("File "+gamma_file+" not found. Stopping execution.")
    """"
    for gene in list_of_genes:
        path="Gamma/gamma_"+model+"_"+gene+".csv"
        if not os.path.isfile(path):
            message.escape("File "+path+ " is missing.")
    """

def check_conserved(config_conserved):
    if not os.path.isfile(config_conserved):
        message.escape("File "+config_conserved+ " not found. Stopping execution.")

def check_file_warning(filename):
    if not os.path.isfile(filename):
        message.warning(f"File {filename} not found. Ignored.")
        return None
    return filename

def check_file_escape(filename, parameter):
    if not os.path.isfile(filename):
        message.escape("File "+filename+" not found (parameter "+parameter+". Stopping execution.")
    elif os.path.getsize(filename) == 0:
        message.escape("File "+filename+" is empty (parameter "+parameter+". Stopping execution.")

def check_limit(limit):
    if limit:
        if not os.path.isfile(limit):
            message.warning("File "+limit+" not found. No limit applied.")
        elif os.path.getsize(limit) == 0:
            message.warning("File "+limit+" is empty. No limit applied.")
            
def check_taxonomy(taxonomy, mandatory=False):
    if taxonomy:
        if not os.path.isfile(taxonomy):
            if mandatory:
                message.escape("File "+taxonomy+" not found (-t). Stopping execution")
            else:
                message.warning("File "+taxonomy+" not found (-t). Ignored.")
                taxonomy=None
        elif  os.path.getsize(taxonomy) == 0:
            if mandatory:
                message.escape("File "+taxonomy+" is empty (-t). Stopping execution")
            else:
                message.warning("File "+taxonomy+" is empty (-t). Ignored.")
                taxonomy=None
    return taxonomy

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
            

def check_and_update_parameters_classify(spectra, taxo, peptide_table, fasta, fasta_dir, limit, deamidation, error, neighbour, allpeptides, mammals, placentals, birds, config):
    """
    Parameters checking and fixing. Configuration of loggers
    """
    config_taxo=taxonomic_model(mammals, placentals, birds)
    peptide_table, taxo, fasta, fasta_dir = check_model_classify(config_taxo, peptide_table, taxo, fasta, fasta_dir)
    if peptide_table:
        check_peptide_table(peptide_table)
    elif fasta or fasta_dir:
        check_sequences(fasta, fasta_dir)
    config=check_config(config, config_taxo)
    check_limit(limit)
    check_spectra_and_error(spectra, error)
    taxo=check_taxonomy(taxo)
    if neighbour not in range(101):
        neighbour=100
        message.warning("Parameter -n (neighbouring): value is 100")
    return spectra, taxo, peptide_table, fasta, fasta_dir, limit, deamidation, error, neighbour, allpeptides, config


def check_and_update_parameters_craft(homology, deamidation, allpeptides, fillin, selection, reconstruction, peptide_table, fasta, fasta_dir, spectra, resolution, limit, taxo, config, placentals, birds, mammals, custom, target, targetfile):
    """
    Parameters checking and fixing for PAMPA CRAFT.
    Configuration of loggers
    """

    param=sum([homology, allpeptides, fillin, selection, reconstruction])
    if param==0 and not deamidation:
         message.escape("Missing parameter: --homology, --allpeptides, --fillin, --deamidation, --selection or reconstruction. Stopping execution")
    if param>1:
        message.escape("Parameters --homology, --allpeptides, --selection, --fillin and --reconstruction are mutually exclusive. Stopping execution")
    check_limit(limit)

    if homology :
        config_taxo=taxonomic_model(mammals, placentals, birds)
        peptide_table, taxo = check_taxonomic_model(config_taxo, peptide_table, taxo)
        if not peptide_table :
            message.escape("No taxonomic model provided. Stopping execution")
        config = check_config(config,config_taxo)
        check_peptide_table(peptide_table)
        check_sequences(fasta, fasta_dir)
        taxo=check_taxonomy(taxo)
        useless_parameters([(spectra, '-s'), (resolution,'-e')])
        
    if deamidation and param==0:
        config = check_config(config)
        check_peptide_table(peptide_table)
        useless_parameters([(fasta, '-f'), (fasta_dir,'-d'), (taxo, '-t'), (spectra,'-s'), (resolution,'-e')])

    if allpeptides:
        config=check_config(config)
        check_sequences(fasta, fasta_dir)
        check_spectra_and_error(spectra, resolution, False)
        useless_parameters([(peptide_table,'-p')])
    
    if fillin:
        config = check_config(config)
        check_peptide_table(peptide_table)
        check_sequences(fasta, fasta_dir, False)
        taxo=check_taxonomy(taxo)
        check_error(resolution, False)
        useless_parameters([(spectra,'-s')])
         
    if selection:
        config = check_config(config)
        check_peptide_table(peptide_table)
        check_spectra_and_error(spectra, resolution)
        useless_parameters([(fasta,'-f'), (fasta_dir,'-d'), (taxo,'-t')])

    if reconstruction:
        config_taxo=taxonomic_model(mammals, placentals, birds)
        peptide_table, taxo = check_taxonomic_model(config_taxo, peptide_table, taxo)
        config = check_config(config,config_taxo)
        if not target and not targetfile:
            message.escape("Missing parameter: -x or -X (target species). Stopping execution")
        if targetfile:
            check_file_escape(targetfile, "-X")
        #if target is None:
        #    message.escape("Missing parameter: -x")
         # TO DO: check config file
        config_substitution, config_gamma, config_conserved = conf.config_matrices_and_co(config)
        check_matrices(config_substitution)
        check_gamma(config_gamma)
        check_conserved(config_conserved)
        check_peptide_table(peptide_table)
        taxo=check_taxonomy(taxo,True)
        check_spectra_and_error(spectra, resolution)

    return homology, deamidation, allpeptides, fillin, selection, reconstruction, peptide_table, fasta, fasta_dir, spectra, resolution, limit, taxo, config, placentals, birds, custom, target, targetfile

