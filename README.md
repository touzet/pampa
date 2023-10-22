# PAMPA 


Pampa (Protein Analysis by Mass Spectrometry for Ancient Species) is a versatile software program designed for species identification in ZooMS (Zooarcheological by Mass Spectrometry) through the usage of marker peptides. It offers the following features :

- Users can define their own markers by creating custom peptide tables. 
- The software facilitates the incorporation of taxonomic data.
- It enables in-depth exploration of various assignment possibilities within the taxonomic space.

## Installation

## How to run the program ?

## Mass spectra



## Peptide tables (-p)

The list of peptide markers used for species identification should be organized within a Tab-Separated Values (TSV) file, featuring the following columns

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

Most of these fields are optional and are here for reference. The following information is mandatory:

- You must provide a taxid for the peptide marker. Rank and taxon names are included primarily to enhance the clarity of results.
- You should furnish either a sequence or a mass for your marker peptide. If the sequence is provided without a mass, the program will automatically compute the mass from it. To do so, it will utilize either the PTM description (when available) or infer potential PTMs from the sequence.

A 'PTM description' provides a concise representation of the number of proline oxidations (P) and deamidations (D) necessary to compute the mass of a peptide sequence. For instance, '2P1D' signifies two oxyprolines and one deamidation, '1D4P' represents one deamidation and four oxyprolines, '2P' corresponds to two oxyprolines without any deamidation, and '1D' indicates one deamidation without oxyprolines. To specify that no PTM should be applied, use '0P' or '0D'.
In cases where the 'PTM' field is omitted, the program will estimate the potential oxyprolines from the peptide sequence, but deamidations will not be considered.


The first row of the file should contain column headings. 


For reference and guidance, sample peptide tables can be found in the 'Peptide_tables' directory.

## Running the program without peptide tables

## Taxonomy

Taxonomy information should be provided in the form of a Tab-Separated Values (TSV) file comprising five columns: Taxid, Common name, Scientific name, Parent (taxid), and Rank (species, genus, etc.). 
You can obtain this type of file directly from UniProt (https://www.uniprot.org/taxonomy) by following these steps:

  1. Use the search bar to find your desired clade, entering its common name, scientific name, or taxid.
  2. Select the clade of interest and click on 'Browse all descendants.'
  3. Locate the 'download' link.
  4. Choose the TSV format and customize the columns in the following order: Common name, Scientific name, Parent, and Rank.
  5. Proceed to download the taxonomy file.

Some examples are available in the directory 'Taxonomy'.
