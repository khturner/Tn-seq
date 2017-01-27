#!/usr/bin/env python
import argparse
import sys
import subprocess

parser = argparse.ArgumentParser(description = 'Perform obligate essentiality analysis on TnSeq insertion site count data')
parser.add_argument('-r', '--reference', help = 'Reference genome to map against (FASTA)', required = True)
parser.add_argument('-g', '--gff', help = 'GFF file describing genomic annotations', required = True)
parser.add_argument('-f', '--features', help = 'Features of interest to consider in essentiality analysis (default: %(default)s)', default = 'gene')
parser.add_argument('-3', '--threeprimetrim', help = 'Percent of the 3\' end of genes to ignore insertions in (default: %(default)s)', default = 10)
parser.add_argument('-i', '--ignore', help = 'Number of highest-read sites to ignore (default: %(default)s)', default = 0)
parser.add_argument('-m', '--minreads', help = 'Minimum number of reads required to consider a site (default: %(default)s)', default = 2)
parser.add_argument('-o', '--output', help = 'Output file prefix', required = True)
parser.add_argument('-c', '--control', help = 'Insertion site read count data TSV file (specify additional times for replicates)', nargs = '+', required = True)
#parser.add_argument('-t', '--test', help = 'Test condition site count data file (specify additional times for replicates)', action = 'append', required = True)

args = parser.parse_args()

cmd = ['Rscript', '--vanilla', sys.path[0] + '/../R/obligate_essentiality_analysis.R', args.reference, args.gff, args.features,
       str(args.threeprimetrim), str(args.ignore), str(args.minreads), args.output] + args.control
print(cmd)
p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
result, err = p.communicate()
if p.returncode != 0:
  raise IOError(err)
print(result)
