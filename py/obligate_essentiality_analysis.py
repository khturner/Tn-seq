#!/usr/bin/env python
import argparse
import sys
import subprocess

parser = argparse.ArgumentParser(description = 'Perform obligate essentiality analysis on TnSeq insertion site count data')
parser.add_argument('-r', '--reference', help = 'Reference genome to map against (FASTA)', required = True)
parser.add_argument('-g', '--gff', help = 'GFF file describing genomic annotations', required = True)
parser.add_argument('-f', '--features', help = 'Features of interest to consider in essentiality analysis (default: %(default)s)', default = 'CDS')
parser.add_argument('-5', '--fiveprimetrim', help = 'Percent of the 5\' end of genes to ignore insertions in (default: %(default)s)', default = 0)
parser.add_argument('-3', '--threeprimetrim', help = 'Percent of the 3\' end of genes to ignore insertions in (default: %(default)s)', default = 10)
parser.add_argument('-i', '--ignore', help = 'Number of highest-read sites to ignore (default: %(default)s)', default = 0)
parser.add_argument('-m', '--minreads', help = 'Minimum number of reads required to consider a site (default: %(default)s)', default = 2)
parser.add_argument('--correctgc', help = 'Include local GC content in smoothing model (default: %(default)s)', default = True)
parser.add_argument('-p', '--pseudoreps', help = 'Number of pseudodata reps to calculate (default: %(default)s, >2,000 suggested for final data)', default = 100)
parser.add_argument('-a', '--attributetag', help = 'Tag in GFF attributes containing gene name (default: %(default)s', default = 'locus_tag')
parser.add_argument('-o', '--output', help = 'Output file prefix', required = True)
parser.add_argument('-c', '--control', help = 'Insertion site read count data TSV file (specify additional times for replicates)', nargs = '+', required = True)

args = parser.parse_args()

if args.correctgc:
  correct_gc_bias = 1
else:
  correct_gc_bias = 0
cmd = ['Rscript', '--vanilla', sys.path[0] + '/../R/obligate_essentiality_analysis.R', args.reference, args.gff, args.features,
       str(args.fiveprimetrim), str(args.threeprimetrim), str(args.ignore), str(args.minreads), str(correct_gc_bias), str(args.pseudoreps),
       args.attributetag, args.output] + args.control
print(cmd)
p = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
result, err = p.communicate()
if p.returncode != 0:
  raise IOError(err)
print(result)
