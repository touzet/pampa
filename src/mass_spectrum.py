"""
   mass_spectrum.py                 

   All utilies to parse and build mass spectra

"""

from pyteomics import mzml, mgf, auxiliary
import csv
import os
import re

from src import message

class Peak(object):
    """
    Representation of a peak in a mass spectrum.
    """

    def __init__(self, mass, intensity):
        self.mass = mass
        self.intensity = intensity

    def __str__(self):
        return f"mass : {self.mass}; intensity : {self.intensity} "


class Spectrum(object):

    def __init__(self, name="", peaks=[], taxid=""):
        self.name = name
        self.peaks = peaks
        self.taxid = taxid

    def __str__(self):
        if len(self.taxid)==0:
            return f"{self.name} : {self.peaks}"
        else:
            return f"{self.name}[{self.taxid}]:{self.peaks}"

    def __len__(self):
        return len(self.peaks)

    def __getitem__(self,i):
        return ((self.peaks)[i])

    def sort(self):
        self.peaks.sort(key=lambda x: x.mass)

    def __max__(self):
        return max({peak.mass for peak in self})

    def __min__(self):
        return min({peak.mass for peak in self})


def parser_binarymatrix(peak_file_name):
    """
    parser for a binary matrix containing a set of mass spectra
    first row: masses
    other rows: spectra (one spectrum per row)
    columns are separated by spaces
    
    Args:
       peak_file_name (str):  path to the peak file  

    Raises:
       NameError, if the input file is not a csv file

    Returns:
       a list of objects Spectrum

    """
    
    list_of_spectra=[]
    f= open (peak_file_name,"r")
    csv_file = csv.reader(f, delimiter=" ")
    firstline=next(csv_file)
    list_of_masses=[float(m.replace(",",".")) for m in firstline[1:]] 
    for row in csv_file:
        new_spectrum=Spectrum()
        new_spectrum.name=row[0]
        new_spectrum.peaks=[(list_of_masses[i-1],0.0) for (i,b) in enumerate(row) if i>0 and int(b)==1] 
        list_of_spectra.append(new_spectrum)
    f.close()
    return list_of_spectra
    

def peak_parser_csv(peak_file_name, name):
    """
    Parser for a mass spectrum in csv format
    first row: heading
    other rows: one peak (mass, intensity) per row
    columns are separated by comma or semi-columns

    Args:
        peak_file (str): path to the csv file
    Returns:
         list of spectra 
    Raises:
        NameError: if the file is not a csv file
    """

    peak_list = []
    f= open (peak_file_name)
    next(f)
    line=1
    try:
        for row in f:
            peak = re.split(',|;', row)
            if len(peak)>1:
                new_peak=Peak( float(peak[0]), float(peak[1]))
                peak_list.append(new_peak)
            elif len(peak)==1:
                new_peak=Peak( float(peak[0]), 0.0)
                peak_list.append(new_peak)
            line += 1
    except  : 
        message.warning("File "+peak_file_name+", line "+str(line)+":  wrong format. File ignored." )
        return Spectrum(name,peak_list,"")
    
    return [Spectrum(name,peak_list,"")]



def peak_parser_mgf(peak_file_name,name):
    """
    Parser for mass spectra in mgf format

    Args:
        peak_file_name (str): path to the mgf file

    Returns:
        list of spectra 

    """
    peak_list = []
    try:
        f = mgf.read(peak_file_name)
        mass_array = f[0]["m/z array"] #pour avoir les m/z des pics
        intensity_array = f[0]["intensity array"] #pour avoir l'intensité des pics
    except  : 
        message.warning("File "+peak_file_name+" is incorrect (wrong format). File ignored." )
        return Spectrum(name,peak_list,"")

    for mass, intensity in zip(mass_array.tolist(), intensity_array.tolist()):
        peak = Peak(float(mass), float(intensity))
        if peak not in peak_list:
            peak_list.append(peak)
    
    return [Spectrum(name,peak_list,"")]


def peak_parser_mzml(peak_file_name,name):
    """
    Parser for spectra in mzML format, based on the module mzml of pyteomics

    Args:
        peak_file_name (str): path to the peak file

    Returns:
        list of Spectrum  

    """
    list_of_spectra=[]
    file_name, _= os.path.splitext(peak_file_name)
    with mzml.MzML(peak_file_name) as reader:
        for spectrum in reader:
            peak_list = []
            name= file_name+" "+(
            spectrum.get('scanList', {})
            .get('scan', [{}])[0]  # Get first scan if available
            .get('name', spectrum.get('id', 'Unknown'))  # Default to 'id' if 'name' is missing
            )
            #name = spectrum.get('name', spectrum.get('id', 'Unknown'))
            #name= spectrum.get('id', 'Unknown')  # Get spectrum ID
            mz_values = spectrum['m/z array']  # Extract m/z values
            intensity_values = spectrum['intensity array']  # Extract intensity values
            for mass, intensity in zip(mz_values, intensity_values):
                peak = Peak(float(mass), float(intensity))
                peak_list.append(peak)
            list_of_spectra.append(Spectrum(name,peak_list,""))
        
       

    '''
    try:
        #f = mzml.MzML(peak_file_name, use_index=True)
        f = mzml.MzML(peak_file_name)
    except  :
        message.warning("File "+peak_file_name+" is incorrect (wrong format). File ignored." )
        return []
    list_of_spectra=[]   
    for s in f:
        mass_array = s["m/z array"] # m/z 
        intensity_array = s["intensity array"] # intensity
        name=s.get('id', 'Unknown')
        peak_list = []
        for mass, intensity in zip(mass_array.tolist(), intensity_array.tolist()):
            peak = Peak(float(mass), float(intensity))
            if peak not in peak_list:
                peak_list.append(peak)
        list_of_spectra.append(Spectrum(name,peak_list,""))
        '''
        
    return  list_of_spectra


def parser(file_path,name):
    """
    Multi-format parser: csv, mgf or mzML. 
    The type of the file is deduced from the filename extension

    Args:
        file_path(str):  path to the file
        name: name of the spectrum

    Returns:
        list of Spectrum

    """
    if os.path.getsize(file_path) == 0:
        message.warning("File "+file_path+" is empty.")
        return []
        
    _, ext = os.path.splitext(file_path)
    if ext == ".csv" or ext == ".CSV" or ext == ".txt" or ext == ".TXT": # csv file 
        list_of_spectra = peak_parser_csv(file_path,name)
    elif ext == ".mgf" or ext == ".MGF" : # mgf file 
        list_of_spectra = peak_parser_mgf(file_path,name)
    elif ext == ".mzML": # mzML file
        list_of_spectra = peak_parser_mzml(file_path,name)
    else:
        message.warning("File "+file_path+" is not recognized (unknown format). Accepted formats are csv, mgf and mzML. File ignored.")
        return []

    return list_of_spectra



