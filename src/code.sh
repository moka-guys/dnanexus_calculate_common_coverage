#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail


#######################################################################
# Make folders
#######################################################################
# output folder
mkdir -p ~/out/trio_wide_coverage/coverage/
# input folder
mkdir input_files

#######################################################################
# Download and move inputs into docker share
#######################################################################
dx-download-all-inputs --parallel

mv $bamfile1_path input_files 
mv $bamfile2_path input_files
mv $bamindex1_path input_files
mv $bamindex2_path input_files
mv $sambamba_bed_path input_files

# test if optional inputs provided before moving
if [[ ! -z $bamfile3_path ]]; then
	mv $bamfile3_path input_files
fi
if [[ ! -z $bamindex3_path ]]; then
	mv $bamindex3_path input_files
fi


#######################################################################
# Set up Docker
#######################################################################
# load docker image
tar -xf calculate_common_coverage_1-0.tar.gz
docker load  --input calculate_common_coverage_1-0.tar

#######################################################################
# build options
# script options:
# '-r' = minimum coverage level
# '-m' (boolean) = merge_overlapping_mate_reads
# '-q' = minimum base quality score
# '-b' = '/path/to/bedfile.bed'
# '-o' = '/path/to/outputfile'
# '-t' = '/path/to/sample1_markdup_sorted.bam'
# '-u' = '/path/to/sample2_markdup_sorted.bam' 
# '-v' = '/path/to/sample3_markdup_sorted.bam'
#######################################################################
# create string of options with required args
opts="-r $coverage_level -q $min_base_qual"

# add optional arguments
# if -m flag is given
if [[ "$merge_overlapping_mate_reads" == true ]]; then
	opts="$opts -m"
fi
# if a third BAM file is provided
if [[ ! -z $bamfile3_path ]]; then
	opts="$opts -v /input_files/$bamfile3_name"
fi

# run docker 
docker run \
-v /home/dnanexus/input_files:/input_files \
-v /home/dnanexus/out/trio_wide_coverage/coverage/:/output_files \
calculate_common_coverage:1.0 \
python3 /calculate_common_coverage/calculate_common_coverage.py \
-b /input_files/$sambamba_bed_name \
-o /output_files/$output_file_name \
-t /input_files/$bamfile1_name \
-u /input_files/$bamfile2_name \
$opts

#######################################################################
# upload output
#######################################################################
dx-upload-all-outputs --parallel