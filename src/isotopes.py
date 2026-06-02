"""
    assignment.py
"""

from pyteomics import mass
from scipy.special import gammaln
import numpy as np
from src import utils
from src import compute_masses, assignment

ISOTOPE_DELTA = 1.003355

# Carbon counts per amino acid
C_TABLE = {
    'A':3,'R':6,'N':4,'D':4,'C':3,'E':5,'Q':5,'G':2,
    'H':6,'I':6,'L':6,'K':6,'M':5,'F':9,'P':5,'S':3,
    'T':4,'W':11,'Y':9,'V':5
}

def count_carbons(peptide):
    comp = mass.Composition(sequence=peptide)
    nC = comp.get('C', 0)
    nC2= sum(C_TABLE[aa] for aa in peptide)
    nC3=0.044 * compute_masses.peptide_mass(peptide)
    if not nC==nC2:
        print("something wrong with the carbons", peptide)
    return nC3


def count_nitrogens(peptide):
    comp = mass.Composition(sequence=peptide)
    nN = comp.get('N', 0)
    return nN

def one_N15_abundance(nN, p=0.00364):
    return nN * p * ((1 - p) ** (nN - 1))

def carbon_envelope(nC, k_max=4, p=0.0107):
    k = np.arange(k_max + 1)

    log_binom = (
            gammaln(nC + 1)
            - gammaln(k + 1)
            - gammaln(nC - k + 1)
    )

    log_prob = (
            log_binom
            + k * np.log(p)
            + (nC - k) * np.log(1 - p)
    )

    return np.exp(log_prob)

def fast_isotope_abundance(peptide):
    nC = count_carbons(peptide)
    envelope = carbon_envelope(nC)
    return envelope

def fast_isotope_abundance_mass(mass):
    nC = mass * 0.044 # number of carbons (approximation)
    envelope = carbon_envelope(nC)
    return envelope

def isotope_abundance(peptide):
    dict_C={}
    print("isotope abundance for", peptide)
    dist = mass.isotopologues(peptide, report_abundance=True,overall_threshold=0.001)
    for isotope in dist:
        print(isotope[0], isotope[1])
        nb_C12=isotope[0].get('C[12]', 0)
        nb_C13 = isotope[0].get('C[13]', 0)
        abundance=dict_C.get((nb_C12, nb_C13), 0)
        if isotope[1]>abundance:
            dict_C[(nb_C12, nb_C13)]=isotope[1]
    print("------")
    max_C12 = max(x for (x,y) in dict_C.keys())
    isotopes=[]
    k=max_C12
    fini=False
    while not fini:
        if (k, max_C12-k) in dict_C:
            isotopes.append(dict_C[(k, max_C12-k)])
        else:
            fini=True
        k+=-1
    for i in dist:
        print(i[0], "mass=" , i.mass, i[1])
    return isotopes


def range_intensity_isotope(peak_list, intensity, isotope0, isotope1):
    expected_intensity=intensity*isotope1/isotope0
    for p in peak_list:
        print(0.5 * expected_intensity," < ", p.intensity , " < ", 4*expected_intensity)
    return any(0.5 * expected_intensity < p.intensity and p.intensity < 4*expected_intensity for p in peak_list)


def check_isotope_pattern2(spectrum, i, mass_markers_list, j, resolution):
    theoretical_mass, markers = mass_markers_list[j]
    print("- - - - Attempt ", theoretical_mass)
    peak_mass, peak_intensity = spectrum[i].mass, spectrum[i].intensity
    j2 = j + 1
    while j2 < len(mass_markers_list) and utils.matching_masses(peak_mass, mass_markers_list[j2][0], resolution):
        markers.update(mass_markers_list[j2][1])
        j2 += 1
    min_intensity = 1.1 * min(peak.intensity for peak in spectrum)
    first_isotopes = assignment.matching_peaks(theoretical_mass + ISOTOPE_DELTA, resolution, spectrum)
    if len(first_isotopes) == 0:
        print("No two isotopes")
        return False, None
    print("from mass")
    isotopes = fast_isotope_abundance_mass(theoretical_mass)
    print(isotopes)
    print("from sequences")
    for peptide in {m.sequence() for m in markers}:
        print("fast: ", fast_isotope_abundance(peptide))
        #print("slow: ", isotope_abundance(peptide))
    isotopes = [x for x in isotopes if peak_intensity * x / isotopes[0] > min_intensity]
    if len(isotopes) < 2:
        print("No two isotopes, bis")
        return False, None
    print("Trying ", spectrum[i].mass)
    print("isotopes", isotopes)
    print("k=0 ", peak_intensity)
    for k in range(1, len(isotopes)):
        print("k=", k, " ", peak_intensity * isotopes[k] / isotopes[0])
    for k in range(2, min(len(isotopes), 4)):
        potential_peaks = assignment.matching_peaks(theoretical_mass + k * ISOTOPE_DELTA, resolution, spectrum)
        if len(potential_peaks) == 0:
            print("fail, no potential peaks")
            return False, None
        if not range_intensity_isotope(potential_peaks, peak_intensity, isotopes[0], isotopes[k]):
            print("fail: intensity range")
            return False, None
    return True, markers
    pre_isotope_peaks = matching_peaks(theoretical_mass - ISOTOPE_DELTA, resolution, spectrum)
    if len(pre_isotope_peaks) == 0:
        return True, markers
    isotopes = fast_isotope_abundance_mass(theoretical_mass -  ISOTOPE_DELTA)
    for peak in pre_isotope_peaks:
        if range_intensity_isotope([spectrum[i]], peak.intensity, isotopes[0], isotopes[1]):
            print("fail previous peak")
            return False, None
    return True, markers

# i : peak index in spectrum
# base_mass: first peak
def check_isotope_pattern(spectrum, i, mass_markers_list, j, resolution):
    theoretical_mass, markers=mass_markers_list[j]
    peak_mass, peak_intensity=spectrum[i].mass, spectrum[i].intensity
    j2=j+1
    while j2<len(mass_markers_list) and utils.matching_masses(peak_mass,mass_markers_list[j2][0], resolution):
        markers.update(mass_markers_list[j2][1])
        j2+=1
    min_intensity=1.1*min(peak.intensity for peak in spectrum)
    first_isotopes=assignment.matching_peaks(theoretical_mass + ISOTOPE_DELTA, resolution, spectrum)
    if len(first_isotopes)==0:
        return False,None
    set_of_peptides = {m.sequence() for m in markers}
    for peptide in set_of_peptides:
        # _=isotope_abundance(peptide)
        isotopes=fast_isotope_abundance(peptide)
    # add N[15] for first isotope
    nN=count_nitrogens(peptide)
    #isotopes[1]+=one_N15_abundance(nN)
    isotopes=[x for x in isotopes if peak_intensity*x/isotopes[0]>min_intensity]
    if len(isotopes)<2:
        return False, None
    print("Trying ", spectrum[i].mass)
    print("isotopes", isotopes)
    print("k=0 ", peak_intensity)
    for k in range(1, len(isotopes)):
        print("k=", k," ", peak_intensity*isotopes[k]/isotopes[0])
    for k in range(2, min(len(isotopes),4)):
        potential_peaks= assignment.matching_peaks(theoretical_mass + k * ISOTOPE_DELTA, resolution, spectrum)
        if len(potential_peaks)==0:
            return False, None
        if not range_intensity_isotope(potential_peaks, peak_intensity, isotopes[0], isotopes[k]):
            print("fail ")
            return False, None
    pre_isotope_peaks = assignment.matching_peaks(theoretical_mass - ISOTOPE_DELTA, resolution, spectrum)
    if len(pre_isotope_peaks)==0:
        print("success", end="\n")
        return True, markers

