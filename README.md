# PAMPA 


Pampa (Protein Analysis by Mass Spectrometry for Ancient Species) is a versatile software program designed for species identification in ZooMS (Zooarcheological by Mass Spectrometry) through the usage of marker peptides. It offers the following features :

- It can handle any number of mass spectra in a single run.
- It allows to automatically generate marker peptides, and users can also define their own markers by creating custom peptide tables. 
- It enables in-depth exploration of various assignment possibilities within the taxonomic space.

## How to install the program ? 

Dependencies: Biopython, pyteomics.

## How to run the program ?

pampa has two modules: 
- _assign_, for species identifications 
- _build_, for construction of custom marker peptides.

Type `pampa assign -h` or `pampa build -h` respectively to print the help.

## Assign module

```
usage: 
 
 pampa assign [-h] (-s SPECTRA PATH) (-e ERROR MARGIN) (-o OUTPUT FILE) (-p PEPTIDE TABLE | -f FASTA file | -d FASTA dir) [-l LIMIT] [-t TAXONOMY] [-n NEIGHBOURING] [-a]

This module is for species identification.

options:
  -h, --help        show this help message and exit

general options:
  -s SPECTRA        Path to the spectra files (one spectrum per file). Authorized formats: cvs, mgd, mzML.
  -e ERROR          Maximal error margin for the observation (in Dalton or ppm). 
  -o OUTPUT         Output path (should include the output file name)

options for organism selection:
  -p PEPTIDE_TABLE  Peptide table (TSV file)
  -f FASTA          Fasta sequences
  -d DIRECTORY      Directory where to find Fasta files.
  -l LIMIT          Limit the set of peptides or fasta sequences to organisms, molecules or sequence ID specified in this file (TSV file), optional.
  -t TAXONOMY       Taxonomy (TSV file), optional.

options for suboptimal solutions:
  -n NEIGHBOUR      Provide near-optimal solutions within a specified percentage margin, ranging between 0 and 100. Default is 100. With this value, only optimal solutions are provided.
  -a                Provide all solutions, and not only suboptimal solutions, within the percentage margin specified with option -n. 

```

### Mass spectra (-s)

The program processes a batch of mass spectra simultaneously. All mass spectra files are contained within the same folder, with one file dedicated to each mass spectrum. These files should have one of the following extensions: .csv (in CSV format), .mgf (in MGF format), or .mzML (in mzML format). Any other files present will be disregarded. You can specify the path to the folder using the '-s' option.

_CSV format_: It consists of two columns. The first column is designated for mass (m/z), and the second column records intensity (I). Columns are separated by either a comma (',') or a semicolon (';'). The initial row serves as the header.

_MGF format_: Mascot Generic Format

_mzML format_: see https://www.psidev.info/mzML

### Error margin (-e)

The error margin is related to the resolution of the mass spectrometer, that is its ability to distinguish closely spaced peaks. We employ it to set an upper bound on the deviation between a peak and the theoretical mass of the marker peptide. This option is mandatory, and can be expressed in Daltons or in ppm.
 -  If the value is smaller than 1, it is assumed to be in Da. In this case, recommended values are  0.1 for maldi TOF, and 0.01 for maldi FT.
 -  If the value is larger than 1, it is assumend to be in ppm. In this case, recommended values are  XXfor maldi TOF, and 5 for maldi FT.

### Output files (-o)

Name of the main output file, in TSV format. This file contains the list of species found for each mass spectrum.
Two other accompanying files are automatically created, in the same directory.

- detail_<outputfile> (TSV file): this file contains the detail of the assignment (which markers are found for which species)
- report_<outputfile> (TXT file): this file contains a report on the run's inputs (number of mass spectra, number of species tested,  parameters...)

### Peptide table (-p)

The first way to use pampa for species identification is to provide a list of marker peptides. This list should be organized within a Tab-Separated Values (TSV) file, featuring the following columns:

- Rank: Taxonomic rank
- Taxid: Taxonomic identifier
- Taxon name: Scientific name
- Sequence: Marker peptide sequence
- PTM: Description of post-translational modifications applied to the marker peptide
- Name: Marker name
- Masses: Peptide mass
- Gene: Gene name, e.g., COL1A1
- SeqId: Sequence identifier(s) of the protein sequence from which the marker peptide is derived
- Begin: Start position of the peptide marker within the protein sequence
- End: End position of the peptide marker within the protein sequence
- Comment: Additional comments about the marker

The first row of the file should contain column headings. 

Most of these fields are optional and are here for reference. The following information is mandatory:

- You must provide a taxid for the peptide marker. Rank and taxon names are included primarily to enhance the clarity of results.
- You should furnish either a sequence or a mass for your marker peptide. If the sequence is provided without a mass, the program will automatically compute the mass from it. To do so, it will utilize either the PTM description (when available) or infer potential PTMs from the sequence.

_PTM description_: This is a concise representation of the number of proline oxidations (P) and deamidations (D) necessary to compute the mass of a peptide sequence. For instance, '2P1D' signifies two oxyprolines and one deamidation, '1D4P' represents one deamidation and four oxyprolines, '2P' corresponds to two oxyprolines without any deamidation, and '1D' indicates one deamidation without oxyprolines. To specify that no PTM should be applied, use '0P' or '0D'.
In cases where the 'PTM' field is omitted, the program will estimate the potential oxyprolines from the peptide sequence, but deamidations will not be considered.

_Examples_: Sample peptide tables can be found in the 'Peptide_tables' folder.

### Running the program without peptide tables (-f, -d and -l options)



### Taxonomy (-t)

The program offers the optional possibility to add taxonomic information to the species identification. In this case, you should supply a taxonomy file.

The taxonomy must be in the form of a Tab-Separated Values (TSV) file comprising five columns: Taxid, Common name, Scientific name, Parent (taxid), and Rank (species, genus, etc.). 
You can obtain this type of file directly from UniProt (https://www.uniprot.org/taxonomy) by following these steps:

  1. Use the search bar to find your desired clade, entering its common name, scientific name, or taxid.
  2. Select the clade of interest and click on 'Browse all descendants.'
  3. Locate the 'download' link.
  4. Choose the TSV format and customize the columns in the following order: Common name, Scientific name, Parent, and Rank.
  5. Proceed to download the taxonomy file.

_Examples_: Some examples are available in the folder 'Taxonomy'.

### Neighbouring (-n and -a)

The option -n allows to obtain also near-optimal solutions for species identification.
The suboptimality range is specified as a percentage, ranging between 0 and 100. 
Default value is 100, meaning that only optimal solutions (with the maximal number of marker petides) are provided. 

The option -a comes with -n. By default, the -n option will generate near-optimal solutions that are not included in any other solution. When the option -a is activated, the program computes all solutions, even those that are included in other solutions.  

## Build module

This module allows to build a new peptide table by homology.

```
usage: 
 
 pampa build [-h] (-p PEPTIDE TABLE) (-o OUTPUT FILE) (-f FASTA file | -d FASTA dir) [-l LIMIT]

This module is for the construction of custom peptide tables.

options:
  -h, --help        show this help message and exit
  -p PEPTIDE_TABLE  TSV file that contains model peptide markers, with sequences.
  -o OUTPUT         Output path (should include the output file name)
  -f FASTA          FASTA file that contains new sequences with header with species.
  -d DIRECTORY      Directory that contains FASTA files
  -l LIMIT          TSV file
  ```

### Peptide table (-p)

This table contains the list of marker peptides that will be used as models to find new markers in new sequences by homology.

### Output files (-o)


### Target sequences (-f, -d and -l)

The target sequences are the amino-acids sequences in which the new markers are searched. Those sequences can be available either in a (multi-)FASTA file (option -f), or in a directory containing FASTA files (option -d). In both cases, the set of sequences can be 
_limited_ to a subset of organisms, molecules or sequence identifiers with option -l.

_Option -f : _ The specified file can contain an arbitrary number of FASTA sequences, coming from various organisms. Two types of FASTA heading are recognized. 
 - Uniprot compliant, with the sequence identifier at the beginning of the heading SeqID and mandatory fields OS, OX and GN:  
   `>P02453 CO1A1_BOVIN Collagen alpha-1(I) chain OS=Bos taurus OX=9913 GN=COL1A1 `

 - NCBI compliant: 

_ Option -d:_



New petides are found by sequence alignment, allowing up to 10% mismatches. Masses are automatically computed.
