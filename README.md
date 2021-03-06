# LRGparser: Usage Instructions

**Parses Locus Reference Genomic (LRG) XML file and outputs CSV file containing information about each exon in each transcript.**

For LRG file details, see Locus Reference Genomic website: http://www.lrg-sequence.org/

**NOTE:** This parser only works for LRG files that had public status as of **November 2016**.

**Authors**

Paul Acklam and Andrew Bond, November 2016


**Requirements:**

* `LRGparser.py`
* Python 2.7
* Internet connection with access to GitHub and LRG XML files

LRGparser.py should be downloaded from the following GitHub repository:

https://github.com/andyb3/LRG_Parse


**Python version**

LRGparser was written and tested in Python 2.7 . It is not currently compatible with Python 3.X.

All modules used are standard Python 2.7 modules, no third party installations should be required.


**Execution**

Navigate to the directory containing LRGparser.py and execute the following command:

&gt;&gt; python `LRGparser.py` [option] ... [gene | LRG] ...

-g [gene symbol]
User provides HGNC approved gene symbol for the gene of interest

-l [LRG code]
User provides the LRG code corresponding to the gene of interest

For example:

    >> python LRGparser.py -g BRCA1

    >> python LRGparser.py -l LRG_292

**Output**

If a valid gene or LRG code is provided by the user, LRGparser downloads the relevant XML file from the LRG website, extracts the required information about each exon in the transcript, and outputs this in a CSV file with the following columns:

LRG_Number | Build | HGNCID | Gene | Transcript_ID | Exon_no | Chrom | Start | End | Sequence

The parser outputs the genomic sequence and coordinates, using the genome build specified as the 'main assembly' in the LRG file. The parser accounts for all differences between the LRG and genomic sequence when converting to the latter.

The file is output into the working directory. The file will have the following name format:
LRG_[number]_output.csv

An example output file (for the NF1 gene) is included in this repository (LRG_214_output_example.csv)

**Testing**

A test suite is available for the program. The easiest way to run the test suite is to download the complete GitHub repository (https://github.com/andyb3/LRG_Parse) and run the tests.py script. **NOTE:** On Windows the repository should be downloaded as a .zip file and extracted, rather than using git clone. (Depending on the way git is configured, the pickle files can become corrupted when using git clone on Windows.)

The program has been successfully tested on Ubuntu 16.10 and Windows 10 (64 bit).

The program output has been validated against the genomic coordinates and sequence found in the Ensembl genome browser for the following genes:
* BRCA1
* NF1
* COL1A1
* MTM1

See the 'EnsemblVsLRGParser...' csv files in the Test_Files directory to view these comparisons.

All currently available LRG files (as of Nov 2016) have been run through the parser (using the allLRGsTest.sh shell script and allLRGs.txt files in Test_Files folder) without the program generating any error messages.
