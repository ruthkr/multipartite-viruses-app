#!/bin/bash

# Script to divide data into subsets
# Authors: Ruth Kristianingsih <ruth.kristianingsih30@gmail.com>

directory=data
infile=$1
outfile_name=$2

# Looping over N
for G in 0.65; do
	awk -v gamma=$G -f 6eq_3d_phase_portrait_clean_points.awk $directory/$infile > $directory/$outfile_name$G.csv;
	awk -v gamma=$G -f 6eq_3d_phase_portrait_clean_separatrix_points.awk $directory/$infile > $directory/$outfile_name$G-separatrix.csv;
	awk -v gamma=$G -f 6eq_3d_phase_portrait_clean_paint.awk $directory/$infile > $directory/$outfile_name$G-paint.csv;
	awk -v gamma=$G -f 6eq_3d_phase_portrait_clean_vectors.awk $directory/$infile > $directory/$outfile_name$G-vectors.csv;
	Rscript 6eq_3d_phase_portrait_clean_surface.R $directory/$outfile_name$G
done

