#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description = 'Perform conditional essentiality analysis on TnSeq insertion site count data')
parser.add_argument('-g', '--gff', help = 'Features to be considered in TnSeq analysis (GFF) - Note: prefilter for features of interest', required = True)
parser.add_argument('-5', '--fiveprimetrim', help = 'Percent of the 5\' end of genes to ignore insertions in', default = 10)
parser.add_argument('-3', '--threeprimetrim', help = 'Percent of the 3\' end of genes to ignore insertions in', default = 10)
parser.add_argument('-c', '--control', help = 'Control condition site count data file (specify additional times for replicates)', action = 'append', required = True)
parser.add_argument('-t', '--test', help = 'Test condition site count data file (specify additional times for replicates)', action = 'append', required = True)
parser.add_argument('-o', '--output', help = 'Output file prefix', required = True)

args = parser.parse_args()
print(args)
