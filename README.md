# PAMPA 


Pampa (Protein Analysis by Mass Spectrometry for Ancient Species) is a versatile software program designed for species identification in ZooMS (Zooarcheological by Mass Spectrometry) through the usage of marker peptides. It offers the following features :

- It can handle any number of mass spectra in a single run.
- Users can define their own markers by creating custom peptide tables. 
- The software facilitates the incorporation of taxonomic data.
- It enables in-depth exploration of various assignment possibilities within the taxonomic space.

## How to run the program ?

```
usage: pampa [-h] -s SPECTRA [-p PEPTIDE_TABLE] [-f FASTA] [-E ENTRY] [-d DIRECTORY] -r RESOLUTION [-t TAXONOMY] [-a AMPLITUDE] -o OUTFILE

options:
  -h, --help            show this help message and exit
  -s SPECTRA, --spectra SPECTRA
                        Path to the spectra files (one spectrum per file).
                        Authorized formats: CSV, MGF, mzML.
  -p PEPTIDE_TABLE, --peptide_table PEPTIDE_TABLE
                        List of peptide markers (TSV file)
  -f FASTA, --fasta FASTA
                        Fasta sequences
  -E ENTRY, --Entry ENTRY
                        List of Fasta files (TSV files)
  -d DIRECTORY, --directory DIRECTORY
                        Directory where to find Fasta files.
  -r RESOLUTION, --resolution RESOLUTION
                        Resolution. Recommended values: 0.01 for maldi FT and 0.1 for maldi TOF.
  -t TAXONOMY, --taxonomy TAXONOMY
                        Taxonomy (TSV file)
  -a AMPLITUDE, --amplitude AMPLITUDE
                        Minimal number of markers (percentage between 0 and 100)
  -o OUTFILE, --outfile OUTFILE
                        Output file name
```

## Mass spectra (-s)

The program processes a batch of mass spectra simultaneously. All mass spectra files are contained within the same folder, with one file dedicated to each mass spectrum. These files should have one of the following extensions: .csv (in CSV format), .mgf (in MGF format), or .mzML (in mzML format). Any other files present will be disregarded. You can specify the path to the folder using the '-s' option.

**CSV format**: It consists of two columns. The first column is designated for mass (m/z), and the second column records intensity (I). Columns are separated by either a comma (',') or a semicolon (';'). The initial row serves as the header.

**MGF format**: Mascot Generic Format

**mzML format**: see https://www.psidev.info/mzML

## Peptide tables (-p)

The list of peptide markers used for species identification should be organized within a Tab-Separated Values (TSV) file, featuring the following columns:

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

**PTM description**: This is a concise representation of the number of proline oxidations (P) and deamidations (D) necessary to compute the mass of a peptide sequence. For instance, '2P1D' signifies two oxyprolines and one deamidation, '1D4P' represents one deamidation and four oxyprolines, '2P' corresponds to two oxyprolines without any deamidation, and '1D' indicates one deamidation without oxyprolines. To specify that no PTM should be applied, use '0P' or '0D'.
In cases where the 'PTM' field is omitted, the program will estimate the potential oxyprolines from the peptide sequence, but deamidations will not be considered.

**Examples**: Sample peptide tables can be found in the 'Peptide_tables' folder.

## Running the program without peptide tables (-f, -E and -d options)

## Resolution (-r)

Resolution refers to the mass spectrometer's ability to distinguish closely spaced peaks. In the program, we employ it to set an upper limit on the deviation between a peak and the theoretical mass of the marker peptide. This option is mandatory. Recommended values are  0.1 for maldi TOF, and 0.01 for maldi FT.

## Taxonomy (-t)

The program offers the optional possibility to add taxonomic information to the species identification. In this case, you should supply a taxonomy file.

The taxonomy must be in the form of a Tab-Separated Values (TSV) file comprising five columns: Taxid, Common name, Scientific name, Parent (taxid), and Rank (species, genus, etc.). 
You can obtain this type of file directly from UniProt (https://www.uniprot.org/taxonomy) by following these steps:

  1. Use the search bar to find your desired clade, entering its common name, scientific name, or taxid.
  2. Select the clade of interest and click on 'Browse all descendants.'
  3. Locate the 'download' link.
  4. Choose the TSV format and customize the columns in the following order: Common name, Scientific name, Parent, and Rank.
  5. Proceed to download the taxonomy file.

**Examples**: Some examples are available in the folder 'Taxonomy'.

## Amplitude (-a)

Default value is 100. 

## Output files (-o)

Name of the main output file, in TSV format. This file contains the list of species found for each mass spectrum.
Two other accompanying files are automatically created:

- //detail_//<outputfile> (TSV file): this file contains the detail of the assignment (which markers are found for which species)
- //report_//<outputfile> (TXT file): this file contains a report on the inputs of the run (number of mass spectra, number of species tested,  parameters...)
