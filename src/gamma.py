import csv
import subprocess
import re
import os
from src import message, sequences

def read_gamma_file(filename, filedir):
    gamma_dict={}
    file_path = os.path.join(filedir, filename)
    with open(file_path, newline='') as f:
        reader = csv.reader(f)
        next(reader)  # Skip the header row
        for row in reader:
            if len(row) >= 2:
                gamma_dict[int(row[0])] = float(row[1].strip())
    return gamma_dict

def compute_gamma_scores_from_set_of_sequences(set_of_sequences, prot, filename):
    file_path=os.path.join("Gamma", "tmp.fasta")
    min_length=min(len(seq.sequence()) for seq in set_of_sequences)
    set_of_tmp_sequences=set()
    for seq in set_of_sequences:
        seq_tmp=sequences.Sequence(field={})
        seq_tmp.field["Sequence"]=seq.sequence()[0:min_length]
        seq_tmp.field["SeqID"]=seq.seqid()
        set_of_tmp_sequences.add(seq_tmp)
    with open(file_path, "w") as tmp_fasta_file:
        for seq in set_of_tmp_sequences:
            tmp_fasta_file.write(f">{seq.seqid()}\n {seq.sequence()}\n")
    r_script_file = "src/Compute_gamma.r"
    try:
        with open("gamma_r_output.log", "w") as outfile:
            subprocess.run(["Rscript", r_script_file, file_path, filename], stdout=outfile, stderr=subprocess.STDOUT, check=True, text=True)
    except subprocess.CalledProcessError:
        message.warning("Computation of gamma function failed for "+prot+". Ignoring this criterion.")
        return None
    finally:
        os.remove(file_path)
    return read_gamma_file(filename, ".")

def convert_gamma_position(gamma_dict, start, length):
    return {k-start:gamma_dict[k] for k in gamma_dict if k in range(start, start+length)}

def most_common(list_of_values):
    counts = {}
    for item in list_of_values:
        counts[item] = counts.get(item, 0) + 1
    return max(counts, key=counts.get)

# attention à bien traiter toutes les erreurs. Notamment si m.protein() n'a pas de fonction gamma
def map_gamma_to_markers(gamma_prot, set_of_markers):
    gamma_marker = {}
    set_of_codes = {m.code() for m in set_of_markers }
    for code in set_of_codes:
        protein = most_common([m.protein() for m in set_of_markers if m.code() == code])
        start_position = most_common([m.helical() for m in set_of_markers if m.code() == code])
        length = most_common([len(m.sequence()) for m in set_of_markers if m.code() == code])
        gamma_marker[code] = convert_gamma_position(gamma_prot[protein], start_position, length)
    return gamma_marker


def initialize_gamma_from_csv_files(list_of_csv_files, set_of_markers):
    gamma_prot = {}
    for filename in list_of_csv_files:
        pattern = r'_([^_.]+)\.'
        m = re.search(pattern, filename )
        GN = m.group(1) if m else None
        gamma_prot[GN]=read_gamma_file(filename,".")
    return map_gamma_to_markers(gamma_prot, set_of_markers)

def initialize_gamma_from_set_of_sequences(set_of_sequences, seed_name):
    list_of_files = []
    gamma_prot = {}
    set_of_prot={seq.protein() for seq in set_of_sequences}
    for prot in set_of_prot:
        file_name="gamma_"+seed_name+"_"+prot+".csv"
        set_of_prot_sequences={seq for seq in set_of_sequences if seq.protein()==prot}
        gamma_prot[prot]=compute_gamma_scores_from_set_of_sequences(set_of_prot_sequences, prot, file_name)
        if gamma_prot[prot] is not None:
            list_of_files.append(file_name)
    return gamma_prot, list_of_files

def initialize_gamma_markers_from_set_of_sequences(set_of_sequences, set_of_markers):
    gamma_prot = initialize_gamma_from_set_of_sequences(set_of_sequences)
    return map_gamma_to_markers(gamma_prot, set_of_markers)
