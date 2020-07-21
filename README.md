# calculate_common_coverage - v0.1

## What does this app do?
This app runs the dockerised [calculate_common_coverage](https://github.com/moka-guys/calculate_common_coverage) script which calculates the % of each gene covered at the required read depth in all samples involved in an analysis.


## What data are required for this app to run?
1. BAM and BAM index files. Up to 3 indexed BAM files can be provided. At least 2 BAM files are required
2. A BED file in the sambamba format.
3. The minimum read depth required (integer).
4. Count overlapping mate reads once? This True/False boolean denotes whether to apply the -m flag. True counts overlapping mate reads once (default) whereas False counts each overlapping mate read seperately.
5. Min base quality to be included (default = 25)
6. Output filename


## What does this app output?
This app produces a single output, which can be used to import gene level coverage into MOKA in the same way the exome coverage data is imported

This file is output to /coverage


## Created by
This app was created within the Viapath Genome Informatics section
