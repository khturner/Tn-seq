#!/usr/bin/env python
import argparse
import time
import subprocess

parser = argparse.ArgumentParser(description = 'Process a TnSeq sequencing run to prepare a count of reads per insertion site, suitable for conditional or absolute essentiality analysis')
parser.add_argument('-r', '--reference', help = 'Reference genome to map against (FASTA)', required = True)
parser.add_argument('-p', '--primer', help = 'Sequence of the primer used to amplify your TnSeq library', required = True)
parser.add_argument('-i', '--invertedrepeat', help = 'Sequence of the Tn end after your primer sequence but before the insertion site', required = True)
parser.add_argument('-m', '--mismatches', help = 'Number of mismatches to allow when searching reads for the primer/repeat signature', default = 1)
parser.add_argument('-1', '--read1', help = 'Read 1 file (FASTQ). This is the read that the Tn end is presumed to be on', required = True)
parser.add_argument('-2', '--read2', help = 'Read 2 file (FASTQ) if paired end sequencing performed. This will be used in mapping')

args = parser.parse_args()
# print(args) # debug

# Report number of reads
print('### ' + time.strftime("%c") + ' ###')
print('Beginning TnSeq read parsing')
p = subprocess.Popen(['wc', '-l', args.read1], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
result, err = p.communicate()
if p.returncode != 0:
  raise IOError(err)
print('Reads in file: ' + str(int(result.strip().split()[0]) / 4))
print('Searching for reads with primer...')

# Report number of reads with primer
print('### ' + time.strftime("%c") + ' ###')
