# PAMPA


PAMPA (Protein Analysis by Mass Spectrometry for Ancient Species) is a versatile software program designed to efficiently manage a wide range of tasks related to handling ZooMS (Zooarcheological by Mass Spectrometry) data through the usage of marker peptides. 

It consists of three primary scripts: 
  - pampa_light: enables fast and easy species identification from mass spectra,
  - pampa_classify:  also dedicated to taxonomic assignment, with a comprehensive range of advanced options,
  - pampa_craft: to build custom sets of marker peptides.

Moreover, PAMPA offers the following features :

- It can handle any number of mass spectra in a single run.
- It allows to automatically generate marker peptides, and users can also define their own markers by creating custom peptide tables. 
- It enables in-depth exploration of various assignment possibilities within the taxonomic space.

## How to install the program ?

PAMPA is written in Python 3.7, and can be installed either by downloading the source code or cloning this repository.  

 - downloading, as a zip file: button _code<>_ on the right-hand side of the screen
 - cloning: `git clone https://github.com/touzet/anc_prot.git`

Additionaly, the full version of PAMPA necessitates the Biopython and pyteomics libraries.  

 - biopython (https://biopython.org/): `pip install biopython`
 - pyteomics (https://pypi.org/project/pyteomics/): `pip install pyteomics`

PAMPA LIGHT requires no external dependencies.

## How to use the program ?

We recommend starting by reviewing the General formats and definitions section, as it covers fundamental concepts shared across all scripts, such as peptide tables, format of mass spectra, error margin, handling PTMs, format of output files, etc. Afterwards, you may refer to the comprehensive documentation of [PAMPA LIGHT](#PAMPA-LIGHT), [PAMPA CLASSIFY](#PAMPA-CLASSIFY) and [PAMPA CRAFT](#PAMPA-CRAFT)) for more in-depth details.


## PAMPA LIGHT

PAMPA LIGHT is designed for rapid and straightforward species identification. It takes a set of mass spectra as input and attempts to determine the best taxonomic assignment for each of them. The assignment utilizes marker peptides from representative species, which are compiled into peptide tables and taxonomies that are included with the program code. These resources can be invoked as preset options for choosing clades. 

```
usage: python3 pampa_light.py 
	[-h]
	(-s SPECTRA PATH)
	(-e ERROR MARGIN)
	(-o OUTPUT FILE)
	(--mammals)
	[-l LIMIT]
	[-n NEIGHBOURING]
	[-a]
```

### Command-Line Example

```
python3 pampa_light.py -s SpectraFolder -e 0.1 -o resultFile.tsv --mammals
```

This command executes the program on all mass spectra within the 'SpectraFolder' directory, compiling the primary results into the TSV file 'resultFile.tsv'. The taxonomic assignment model utilized is 'mammals', and the error margin for masses is set to 0.1.

### Mass spectra (-s) 


__This option is mandatory.__

The program processes a batch of mass spectra simultaneously. All mass spectra files are contained within the same folder, with one file dedicated to each mass spectrum. These files should have one of the following extensions: .csv or .txt (in CSV format), .mgf (in MGF format), or .mzML (in mzML format). Any other files present will be disregarded. You can specify the path to the folder using the '-s' option.


### Error margin (-e)

__This option is mandatory.__

The error margin is related to the resolution of the mass spectrometer, that is its ability to distinguish closely spaced peaks. We employ it to set an upper bound on the deviation between a peak and the theoretical mass of the marker peptide. This option is mandatory, and can be expressed in Daltons or in ppm.
 -  If the value is smaller than 1, it is assumed to be in Da (Daltons). In this case, recommended values are  0.1 for maldi TOF, and 0.01 for maldi FT.
 -  If the value is larger than 1, it is assumed to be in ppm (parts per million). In this case, recommended values are 50 for maldi TOF, and 5 for maldi FTICR.

### Output files (-o)

__This option is mandatory.__

This option enables you to specify the name of the main output file, in TSV format. You can include the desired path to the directory where the file should be created. If the specified directory does not exist, it will be created automatically.

For each spectrum, the output file will give the best assignment, based on the highest number of marker peptides. It contains the following information:
- Peaks from the spectrum that match the marker petides,  
- Score: the total number of marker peptides, 
- Assignment: largest subtree of the taxonomy that is compatible with the marker peptides found,
- Rank: Taxonomic rank of the assignment (e.g. species, genus, family),
- Species: the list of species supporting the assignment.

Two other accompanying files are automatically created, in the same directory.

- detail_&lt;outputfile&gt; (TSV file): this file contains the detail of the assignment (which markers are found for which species). It also provides the intensity of the peaks used in the assignment.
- report_&lt;outputfile&lt; (TXT file): this file contains a report on the run's inputs (number of mass spectra, number of species tested,  parameters...)

### <a id="peptide"></a>Organism selection (-- mammals)

__This option is mandatory.__

PAMPA LIGHT utilizes a predefined set of marker peptides in conjunction with the NCBI taxonomy for species identification. The markers are accessible through _peptide tables_, which are stored in TSV files distributed with the code (see the [Peptide tables](#Peptide-tables) section), together with the taxonomy (see the [Taxonomy](#Taxonomy) section).

### <a id="limit"></a> Limiting search (-l)

It is possible to filter the peptide table to limit the search according to various criteria such as organism, gene name, sequence identifier, or PTMs. For that, you can use the '-l' option together with a _limit file_ (see the [Limiting searches](#Limiting-searches) section). 

### Neighbouring (-n and -a)

By default, PAMPA identifies the species with the highest number of marker peptides. 
The -n option allows to obtain also near-optimal solutions. For that, you can set the suboptimality 
range as a percentage from 0 to 100, with the default being 100 
(corresponding to solutions with the highest number of marker peptides). 
For example, if the optimal solutions has 11 marker peptides, '-n 80' will provide solutions with 9 markers or more.

By default, the '-n' option will generate only near-optimal solutions that are not included in any other solution.
When used together with '-n,'  the '-a' option allows to change this, so that the program computes all solutions, even those that are included in other solutions.  

## General formats and definitions

### Mass spectra 

PAMPA recognizes three formats for mass spectra:

  - _CSV format_: It consists of two columns. The first column is designated for mass (m/z), and the second column records intensity (I). Columns are separated by either a comma (',') or a semicolon (';'). The initial row serves as the header.

  - _MGF format_: Mascot Generic Format

  - _mzML format_: see https://www.psidev.info/mzML

In all cases, we recommend deisotoping the mass spectra before processing them.

In practice, the program can process a batch of mass spectra simultaneously. All mass spectra files are contained within the same folder, with one file dedicated to each mass spectrum. These files should have one of the following extensions: .csv or .txt (in CSV format), .mgf (in MGF format), or .mzML (in mzML format). Any other files present will be disregarded. 

### Peptide tables

PAMPA handles set of markers that serves to species identification.
Those markers are accessible through _peptide tables_, which are stored in TSV files distributed with the code. Peptide tables should include 12 columns corresponding to the 12 fields below: 

- Rank: Taxonomic rank
- Taxid: Taxonomic identifier
- Taxon name: Scientific name
- Sequence: Marker peptide sequence
- PTM: Description of post-translational modifications (PTMs) applied to the marker peptide (see below)
- Name: Marker name
- Masses: Peptide mass
- Gene: Gene name, e.g., COL1A1
- SeqId: Sequence identifier(s) of the protein sequence from which the marker peptide is derived
- Begin: Start position of the peptide marker within the protein sequence
- End: End position of the peptide marker within the protein sequence
- Comment: Additional comments about the marker

The first row of the file should contain column headings. 

Most of these fields are optional and are here for reference and traceability. Only the following information is mandatory:
 - You must provide a taxid for the peptide marker. Rank and taxon names are included primarily to enhance the clarity of results.
 - You should furnish either a sequence, possibly with a PTM description,  or a mass for your marker peptide. If the sequence is provided without a mass, the program will automatically compute the mass from it. To do so, it will utilize either the PTM description (when available) or infer potential PTMs from the sequence.

Several pre-defined peptide tables  are distributed with the code of PAMPA. You van find them in the Peptide_table directrory. Do not move them or modify them, since they are used by PAMA LIGHT. Alternatively, peptide table files can be created manually with any spreadsheet software by opting for the TSV export format. Lastly,  [PAMPA CRAFT](#PAMPA-CRAFT) provides automated methods to generate  peptide tables.  

### PTM description 

PAMPA recognizes three types of PTMs: 
 - oxylation of prolines (indicated by the single-letter code 'O'), 
 - deamidation of asparagine and glutamine (indicated by the single-letter code 'D'), 
 - phosphorylation of serine, threonine, and tyrosine (indicated by the single-letter code 'P'). 

The _PTM description_ is  a concise representation of the number of oxylations, deamidations and phosphorylations necessary to compute the mass of a peptide sequence. For instance, '2O1D' signifies two oxyprolines and one deamidation, '1P4O' represents one phosphorylation and four oxyprolines, '2O' corresponds to two oxyprolines without any deamidation and phosphorylation. When no PTM applies, the description should be "0O", or "0D", etc. 

When the PTM description this is left empty in the peptide table,  it is assumed that PTMs are not known. In this case, they are directly infered by PAMPA following two rules:
 
  - No deamidation and phosphorylation are added.
  - The number of oxyprolines is determined empirically using the following formula: Let 'p' represent the total number of prolines in the peptide, and 'pp' represent the number of prolines involved in the pattern 'GxP'. If the difference 'p-pp' is less than 3, then 'pp' oxyprolines are applied. If 'p-pp' is 3 or greater, 'pp' oxyprolines and 'pp+1' oxyprolines are applied.


### Taxonomy

PAMPA offers the optional possibility to add taxonomic information to the species identification. In this case, it uses a taxonomy file.
Such files are available in the Taxonomy folder. Alternatively, users can create their own taxonomy file. 
The taxonomy must be in the form of a Tab-Separated Values (TSV) file comprising five columns: Taxid, Common name, Scientific name, Parent (taxid), and Rank (species, genus, etc.). 
You can obtain this type of file directly from UniProt (https://www.uniprot.org/taxonomy) by following these steps:

  1. Use the search bar to find your desired clade, entering its common name, scientific name, or taxid.
  2. Select the clade of interest and click on 'Browse all descendants.'
  3. Locate the 'download' link.
  4. Choose the TSV format and customize the columns in the following order: Common name, Scientific name, Parent, and Rank.
  5. Proceed to download the taxonomy file.

### FASTA sequences

Two types of FASTA heading are recognized. 

 - UniprotKB-like, with some sequence identifier at the beginning of the heading, and mandatory fields OS (scientific name of the organism) , OX (taxonomomic identifier of the organism, such as assigned by the NCBI) and GN (Gene Name):
     
   `>P02453 CO1A1_BOVIN Collagen alpha-1(I) chain OS=Bos taurus OX=9913 GN=COL1A1 `

 - NCBI-like : to do

### Limiting searches

It is possible to filter the peptide table or the peptide sequences to limit the search according to various criteria such as organism, gene name, sequence identifier, or PTMs. For that, you can that follow  these guidelines:

- Create a text file that outlines your filtering criteria.
- Each line in the file should comprise one or several of these fields:
	- "OS=" for authorized taxon names  
 	- "OX=" for authorized taxids	
	- "GN=" for authorized gene names
	- "PTM=" for authorized PTMs
	- "SeqID=" for authorized sequence identifiers
- Separate multiple elements for a field with commas.

For example, If you want to limit your search to a specific set of organisms, your file might look like this: 
```
OS= Diceros bicornis, Cervus elaphus, Bos taurus, Equus caballus   
```   
Of course, you can combine constraints to narrow down your search. For instance, limiting the search to markers coming from COL1A2 gives:
```
OS= Diceros bicornis, Cervus elaphus, Bos taurus, Equus caballus GN=COL1A2
```  
This means that the search will focus on markers from COL1A2 within the specified organisms.

Assume now that you want to further refine this selection and  exclude certain PTMs, such as deamidation and  phosporylation. Then you have to add one contraint to authorize only proline oxylation. This gives:
```
OS= Diceros bicornis, Cervus elaphus, Bos taurus, Equus caballus GN=COL1A1 PTM=O
```  
The limit file can comprise an arbitrary number of lines, with each line representing a distinct constraint. The resulting selection is determined by the union of all these constraints.

Finally, in the presence of a taxonomy (as is the case with PAMPA light), the OS and OX fields become applicable to clades at any taxonomic rank (e.g. genus, family, order). In such instances, the constraint will choose all descendants accordingly. 

```
OS=Pecora GN=COL1A1
```
This will select all COL1A1 markers for species from the Pecora infraorder. Equivalently, you might have used the taxid of Pecora, with the OX field.
```
OX=35500 GN=COL1A1
```

### Error margin 

The error margin is related to the resolution of the mass spectrometer, that is its ability to distinguish closely spaced peaks. We employ it to set an upper bound on the deviation between a peak and the theoretical mass of the marker peptide. This option is mandatory, and can be expressed in Daltons or in ppm.
 -  If the value is smaller than 1, it is assumed to be in Da (Daltons). In this case, recommended values are  0.1 for maldi TOF, and 0.01 for maldi FT.
 -  If the value is larger than 1, it is assumed to be in ppm (parts per million). In this case, recommended values are 50 for maldi TOF, and 5 for maldi FTICR.



## PAMPA CLASSIFY

PAMPA CLASSIFY represents an advanced version of PAMPA light, offering users the capability to utilize their personalized set of marker peptides for taxonomic assignment. 
This functionality is enabled by custom peptide tables, which are TSV files containing list of marker peptides. Users have the option to create these tables manually following the description given in the [PAMPA light](#PAMPA-light) section, or generate them using [PAMPA CRAFT](#PAMPA-CRAFT). Additionally, in situations where no marker peptides are available, it is possible to supply FASTA sequences for the automatic inference of peptides through in silico digestion.

```
usage: 
 
 python3 pampa_classify.py 
	[-h]
	(-s SPECTRA PATH)
	(-e ERROR MARGIN)
	(-o OUTPUT FILE)
	(-p PEPTIDE TABLE... | -f FASTA file | -d FASTA dir)
	[-l LIMIT]
	[-t TAXONOMY]
	[-n NEIGHBOURING]
	[-a]

options:
  -h, --help        show this help message and exit

general options:
  -s SPECTRA        Path to the spectra files (one spectrum per file). Authorized formats: cvs, mgd, mzML.
  -e ERROR          Maximal error margin for the observation (in Dalton or ppm). 
  -o OUTPUT         Output path (should include the output file name)

options for organism selection:
  -p PEPTIDE_TABLE  Peptide table (TSV file)
  -f FASTA          Fasta sequences, for in silico digestion
  -d DIRECTORY      Directory where to find Fasta files
  -l LIMIT          Limit the set of peptides or fasta sequences to organisms, molecules or sequence ID specified
                    in this file (TSV file), optional.
  -t TAXONOMY       Taxonomy (TSV file), optional.

options for suboptimal solutions:
  -n NEIGHBOUR      Provide near-optimal solutions within a specified percentage margin, ranging between 0 and 100.
                    Default is 100. With -n 100, only optimal solutions are provided.
  -a                Provide all solutions, and not only suboptimal solutions, within the percentage margin specified
                    with option -n. 

```

The  **-s** (mass spectra), **-e** (error margin), **-o** (output), **-l** (limit), **-n** (neighbouring) and **-a**  options are the same as with PAMPA light and the documentation can be found in the [dedicated section](#PAMPA-light).

Hereafter, we provide description of **-p**, **-f**, **-d** and **-t** options,  each of which is specific to this module.

### Peptide table (-p)

This option allows you to employ your own set of marker peptides. This set should be structured within a _peptide table_, formatted as a TSV (Tab-Separated Values) file. The specific format details for peptide tables are descrided  in the [PAMPA light](#PAMPA-light) section.  Such file can be created manually with any spreadsheet software by opting for the TSV export format. Alternatively, the module [PAMPA CRAFT](#PAMPA-CRAFT) provides automated methods to generate  peptide tables.  

Peptide tables should contain 12 columns corresponding to  Rank, Taxid, Taxon name, Sequence, PTM, Name, Masses, Gene, SeqID, Begin, End, Comment. Most of these fields are optional and are here for reference and traceability. Only the following information is mandatory:

- You must provide a taxid for the peptide marker. Rank and taxon names are included primarily to enhance the clarity of results.
- You should furnish either a sequence, possibly with a PTM description,  or a mass for your marker peptide. If the sequence is provided without a mass, the program will automatically compute the mass from it. To do so, it will utilize either the PTM description (when available) or infer potential PTMs from the sequence.

Note that you can specify multiple peptide tables by using the -p option repeatedly.  

_Inference of PTMs_:  TO DO 

_Examples_: Sample peptide tables can be found in the 'Peptide_tables' folder.

### Running the program without peptide tables (-f and -d)

When no marker peptides are available, it is possible to provide FASTA sequences for the representative species instead. These sequences will undergo in silico digestion to identify all tryptic peptides, allowing for up to one missed cleavage. Masses are then automatically computed. 

The provided sequences can be available either in a (multi-)FASTA file ('-f'option), or in a directory containing FASTA files ('-d' option). In both cases, the set of sequences can optionnally be _limited_ to a subset of organisms, molecules or sequence identifiers with '-l' option.

_Option -f :_ The specified file can contain an arbitrary number of FASTA sequences, coming from various organisms. Two types of FASTA heading are recognized. 
 - UniprotKB-like, with some sequence identifier at the beginning of the heading, and mandatory fields OS (scientific name of the organism) , OX (taxonomomic identifier of the sorganism, such as assigned by the NCBI) and GN (Gene Name):
     
   `>P02453 CO1A1_BOVIN Collagen alpha-1(I) chain OS=Bos taurus OX=9913 GN=COL1A1 `

 - NCBI-like 

_Option -d:_ The directory can contain an arbitrary number of FASTA files, following the same requirements as with '-f' option.
Only files with extension _.fa_ or _.fasta_ will be examined. 


### Taxonomy (-t)

The program offers the optional possibility to add taxonomic information to the species identification. In this case, you should supply a taxonomy file.

The taxonomy must be in the form of a Tab-Separated Values (TSV) file comprising five columns: Taxid, Common name, Scientific name, Parent (taxid), and Rank (species, genus, etc.). 
You can obtain this type of file directly from UniProt (https://www.uniprot.org/taxonomy) by following these steps:

  1. Use the search bar to find your desired clade, entering its common name, scientific name, or taxid.
  2. Select the clade of interest and click on 'Browse all descendants.'
  3. Locate the 'download' link.
  4. Choose the TSV format and customize the columns in the following order: Common name, Scientific name, Parent, and Rank.
  5. Proceed to download the taxonomy file.

When a taxonomy is provided, the software will indicate, for each spectrum, the taxonomic resolution of the assignment. This is computed as the largest clade of the taxonomy that is compatible with the prediction.

_Examples_: example files are available in the folder 'Taxonomy'.


## PAMPA CRAFT 

PAMPA CRAFT  is for the design of custom peptide tables, that can then be utilized by [PAMPA CLASSIFY](#PAMPA-CLASSIFY).

```
usage: 
pampa_craft  [-h] 
   --homology | --denovo | --fillin 
   [-f FASTA | -d DIRECTORY] [-p PEPTIDE_TABLE] [-l LIMIT] -o OUTPUT 
```

The three options —-homology, --denovo, and --fillin represent three distinct approaches for constructing a new peptide table: 
  - by homology from a set of existing peptide markers,
  - de novo, using in silico tryptic digestion for a collection of protein sequences,
  - by filling in an existing peptide table, for which masses are missing.

The full usage description is given below.

### --homology 

With this option, the input consists of a set of well-defined marker peptides, and the goal is to search a set of target protein sequences for similar peptides. New peptides are discovered through sequence alignment, allowing up to 10% mismatches between peptides. The algorithm also ensures that the new peptides can undergo tryptic digestion and infers new cleavage sites when necessary. Masses are automatically computed. 

'''
 usage:
 python3 pampa_craft --homology
   -p PEPTIDE_TABLE [PEPTIDE_TABLE]
   Peptide table(s) that contain model peptide markers
   -f FASTA      Fasta file for new species
   -d DIRECTORY  Directory containing Fasta files for new species
   -l LIMIT      Limit file that 
   -o OUTPUT     Path to the output file (new peptide table)
'''

#### Peptide table (-p)

This table contains the list of marker peptides that will be used as models to find new markers in new sequences by homology.
The format of this table is described in section [Organism selection](#peptide). This option is required.


#### Target sequences (-f, -d and -l)

__The use of either -f or -d is mandatory.__

The target sequences are the amino-acids sequences in which the new markers are searched. Those sequences can be available either in a (multi-)FASTA file (-f option), or in a directory containing FASTA files (-d option). In both cases, the set of sequences can optionnally be _limited_ to a subset of organisms, molecules or sequence identifiers with -l option.

_Option -f :_ The specified file can contain an arbitrary number of FASTA sequences, coming from various organisms. Two types of FASTA heading are recognized. 
 - UniprotKB-like, with some sequence identifier at the beginning of the heading, and mandatory fields OS (scientific name of the organism) , OX (taxonomomic identifier of the sorganism, such as assigned by the NCBI) and GN (Gene Name):
     
   `>P02453 CO1A1_BOVIN Collagen alpha-1(I) chain OS=Bos taurus OX=9913 GN=COL1A1 `

 - NCBI-like 

_Option -d:_ The directory can contain an arbitrary number of FASTA files, following the same requirements as with '-f' option.
Only files with extension _.fa_ or _.fasta_ will be examined. 

_Option -l:_ This option allows to filter the set of FASTA sequences to limit the selection according to the organism (OS=), the taxid (OX=), the gene name (GN=), the sequence identifier (SeqID=). The full description of the syntax is given in section [Limiting search](#limit).

#### Output file (-o)

This is the name of the new peptide table created by the program. This option is required.

### De novo

This option allows to infer all tryptic peptides from a set of FASTA sequences through in silico digestion, and to compute their masses. 


```
usage:
pampa_craft --denovo 
   -p PEPTIDE_TABLE [PEPTIDE_TABLE]
   Peptide table(s) that contain model peptide markers
   -f FASTA      Fasta file for new species
   -d DIRECTORY  Directory containing Fasta files for new species
   -l LIMIT      Limit file that 
   -o OUTPUT     Path to the output file (new peptide table)
```

#### Target sequences (-f, -d and -l)

The target sequences consist of amino acid sequences in FASTA format that will undergo tryptic digestion. These sequences can be provided either as individual sequences within a (multi-)FASTA file (using the -f option) or as multiple FASTA files within a directory (using the -d option). The use of either -f or -d is required. Additionally, users have the option to selectively limit the set of sequences to a specific subset of organisms, molecules, or sequence identifiers using the -l option.

_Option -f :_ The specified file can contain an arbitrary number of FASTA sequences, coming from various organisms. Two types of FASTA heading are recognized. 
 - UniprotKB-like, with some sequence identifier at the beginning of the heading, and mandatory fields OS (scientific name of the organism) , OX (taxonomomic identifier of the sorganism, such as assigned by the NCBI) and GN (Gene Name):
     
   `>P02453 CO1A1_BOVIN Collagen alpha-1(I) chain OS=Bos taurus OX=9913 GN=COL1A1 `

 - NCBI-like 

_Option -d:_ The directory can contain an arbitrary number of FASTA files, following the same requirements as with '-f' option.
Only files with extension _.fa_ or _.fasta_ will be examined. 

_Option -l:_ This option allows to filter the set of FASTA sequences to limit the selection according to the organism (OS=), the taxid (OX=), the gene name (GN=), the sequence identifier (SeqID=). The full description of the syntax is given in section [Limiting search](#limit).

#### Output file (-o)

This is the name of the new table containing tryptic peptides created by the program. This option is required. 

### Fill in 



```
 usage:
 python3 pampa_craft --fillin 
   -p PEPTIDE_TABLE Peptide table to complete
   -o OUTPUT        Path to the output file (new peptide table)
```

#### Peptide table (-p)

Name of the peptide table to complete.

#### Output file (-o)

This is the name of the new table obtained by completion of the input table.  This option is required. 


## Bug Report

PAMPA is still under development. If you come across any bugs or unexpected behavior, please take a moment to report it using our GitHub Issues page. Just click on the "Issues" tab above and create a new issue. You may also contact the author (helene.touzet@univ-lille.fr).

We value your all feedback and contributions.  Thank you for helping us make the PAMPA project better!

