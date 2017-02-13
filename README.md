# iHS_calc

# Version log

v1.4 (used in Johnson et al):
 
- updated code used in Voight et al (2006) with the following:
 a. Outputs all SNP tests, with flags to denote if the iHS calculation was kosher 
    (e.g. if close to the edge, provide a warning)
 b. Outputs 
 c. Slight change in iHS calculation resulting in a dramatic speedup with no loss in accuracy.

# Compilation instructions

Use the attached makefile, and run make on the command line

%> make

to make sure the program works, run

%> ./iHS_calc compilecheck

if you read "OK!" then the compilation works and the calculation is ready to run.

# Required Files.

The calculator requires two files.

"INFO" file: see sample_info. This is a map file. Each row is a marker. The columns contain:

 Column 1: the SNP identifier
 Column 2: the physical position (not used)
 Column 3: The recombination position, in units of rho (=4Nr) the population recombination rate parameter
 Column 4, 5: Alleles for the marker (not used)

"DATA" file: see sample_data. This is a flat file (no headers); each row is a chromosome, each column is an allele at the given marker. 
The markers are polarized by ancestral state. "0" is ancestral, "1" is derived.   

# Command line instructions

This program is also used in conjunction with WHAMM, a script/tool (still in progress) which
can perform some of the processing steps (automation by chromosome, filter markers, etc.

The calculator can be used in two modes.

Mode One. Provide an info and data file

  %> ./iHS_calc sample_info sample_data

this will analyze all markers in the given files.

Mode Two. Focus on a subset of the data

  %> ./iHS_calc sample_info sample_data 0 9

this will analyze only the first 10 markers in the data file provided.
