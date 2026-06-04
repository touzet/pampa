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

def check_model(birds, placentals, mammals, peptide_table, taxonomy, config):
    param = sum([placentals, birds, mammals])
    if param > 1:
        message.escape(
            "Taxonomic models placentals, birds and mammals are mutually exclusive. Stopping execution.")
    if param == 0 and not config:
        message.escape(
            "No taxonomic model provided. Stopping execution.")
    if placentals:
        model = "placentals"
    elif birds:
        model = "birds"
    elif mammals:
        model = "mammals"
    if config:
        message.warning(f"Parameter -c overwrites the --{model} mode. Using your config file {config} instead.")
    else:
        config = "config_" + model + ".json"
    if peptide_table:
        message.warning(f"Parameter -p overwrites the --{model} mode. Using your peptide table.")
    else:
        peptide_table = conf.config_peptide_table(config)
    if taxonomy:
        message.warning(f"Parameter -t overwrites the --{model} mode. Using your taxonomy.")
    else:
        taxonomy = conf.config_taxonomy(config)
    return config, peptide_table, taxonomy

def check_model_classify(birds, placentals, mammals, peptide_table, taxonomy, fasta, fasta_dir,config):
    param = sum([placentals, birds, mammals])
    if param > 1:
        message.escape(
            "Taxonomic models placentals, birds and mammals are mutually exclusive. Stopping execution.")
    if param==1:
        if placentals:
            model = "placentals"
        elif birds:
            model = "birds"
        elif mammals:
            model = "mammals"
        if config:
            message.warning("Parameter -c overwrites the --" + model + " mode. Using your config file "+ config+" instead.")
        else:
            config = "config_" + model + ".json"
        if peptide_table or fasta or fasta_dir:
            message.warning("Parameter -p, -f and -d overwrite the --" + model + " mode. Using your peptide table.")
        else:
            peptide_table = conf.config_peptide_table(config)
        if taxonomy:
            message.warning("Parameter -t overwrites the --" + model + " mode. Using your taxonomy.")
        else:
            taxonomy = conf.config_taxonomy(config)
    return config, peptide_table, taxonomy



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
            

def check_and_update_parameters_classify(spectra, taxonomy, peptide_table, fasta, fasta_dir, limit, deamidation, error, neighbour, allpeptides, mammals, placentals, birds, config):
    """
    Parameters checking and fixing. Configuration of loggers
    """
    param = sum([placentals, birds, mammals])
    if param>0:
        config, peptide_table, taxonomy = check_model_classify(birds, placentals, mammals, config, peptide_table, taxonomy, fasta, fasta_dir)
    else:
        if peptide_table:
            if fasta or fasta_dir:
                message.escape(
                    "Options -p (peptide_table), -f (fasta) and -d (directory of fasta files) are mutually incompatible. Stopping execution")
            else:
                check_peptide_table(peptide_table)
        elif fasta or fasta_dir:
            check_sequences(fasta, fasta_dir)
        else:
            message.escape("Missing parameter for marker peptides (-p, -f or -d). Stopping execution")

    config=check_config(config)
    check_limit(limit)
    check_spectra_and_error(spectra, error)
    taxonomy=check_taxonomy(taxonomy)
    
    if neighbour not in range(101):
        neighbour=100
        message.warning("Parameter -n (neighbouring): value is 100")

    return spectra, taxonomy, peptide_table, fasta, fasta_dir, limit, deamidation, error, neighbour, allpeptides, config


def check_and_update_parameters_craft(homology, deamidation, allpeptides, fillin, selection, reconstruction, peptide_table, fasta, fasta_dir, spectra, resolution, limit, taxonomy, config, placentals, birds, mammals, custom, target, targetfile):
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
        config, peptide_table, taxonomy = check_model(birds, placentals, mammals, config, peptide_table, taxonomy)
        config = check_config(config)
        check_peptide_table(peptide_table)
        check_sequences(fasta, fasta_dir)
        taxonomy=check_taxonomy(taxonomy)
        useless_parameters([(spectra, '-s'), (resolution,'-e')])
        
    if deamidation and param==0:
        config = check_config(config)
        check_peptide_table(peptide_table)
        useless_parameters([(fasta, '-f'), (fasta_dir,'-d'), (taxonomy, '-t'), (spectra,'-s'), (resolution,'-e')])

    if allpeptides:
        config=check_config(config)
        check_sequences(fasta, fasta_dir)
        check_spectra_and_error(spectra, resolution, False)
        useless_parameters([(peptide_table,'-p')])
    
    if fillin:
        config = check_config(config)
        check_peptide_table(peptide_table)
        check_sequences(fasta, fasta_dir, False)
        taxonomy=check_taxonomy(taxonomy)
        check_error(resolution, False)
        useless_parameters([(spectra,'-s')])
         
    if selection:
        config = check_config(config)
        check_peptide_table(peptide_table)
        check_spectra_and_error(spectra, resolution)
        useless_parameters([(fasta,'-f'), (fasta_dir,'-d'), (taxonomy,'-t')])

    if reconstruction:
        config, peptide_table, taxonomy = check_model(birds, placentals, mammals, config, peptide_table, taxonomy)
        config = check_config(config)
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
        taxonomy=check_taxonomy(taxonomy,True)
        check_spectra_and_error(spectra, resolution)

    return homology, deamidation, allpeptides, fillin, selection, reconstruction, peptide_table, fasta, fasta_dir, spectra, resolution, limit, taxonomy, config, placentals, birds, custom, target, targetfile

