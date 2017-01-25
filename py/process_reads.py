#!/usr/bin/env python
import argparse
import time
import subprocess
import os.path

parser = argparse.ArgumentParser(description = 'Process a TnSeq sequencing run to prepare a count of reads per insertion site, suitable for conditional or absolute essentiality analysis')
parser.add_argument('-r', '--reference', help = 'Reference genome to map against (FASTA)', required = True)
parser.add_argument('-p', '--primer', help = 'Sequence of the primer used to amplify your TnSeq library', required = True)
parser.add_argument('-i', '--invertedrepeat', help = 'Sequence of the Tn end after your primer sequence but before the insertion site', required = True)
parser.add_argument('-m', '--mismatches', help = 'Number of mismatches to allow when searching reads for the primer/repeat signature', default = 1)
parser.add_argument('-1', '--read1', help = 'Read 1 file (FASTQ). This is the read that the Tn end is presumed to be on', required = True)
parser.add_argument('-2', '--read2', help = 'Read 2 file (FASTQ) if paired end sequencing performed. This will be used in mapping')
parser.add_argument('-t', '--threads', help = 'Number of compute threads', default = 1)
parser.add_argument('-o', '--output', help = 'Output file prefix', required = True)

args = parser.parse_args()
# print(args) # debug

# Report number of reads
print('### ' + time.strftime("%c") + ' ###')
print('Beginning TnSeq read parsing...')
p = subprocess.Popen(['wc', '-l', args.read1], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
result, err = p.communicate()
if p.returncode != 0:
  raise IOError(err)
print('Reads in file: ' + str(int(result.strip().split()[0]) / 4))
print('')

# Report number of reads with primer
print('### ' + time.strftime("%c") + ' ###')
print('Searching for reads with primer...')
#DEBUG: skip computation
#p = subprocess.Popen(['fqgrep', '-m', str(args.mismatches), '-C', '-p', args.primer, args.read1],
#                     stdout = subprocess.PIPE, stderr = subprocess.PIPE)
#result, err = p.communicate()
#if p.returncode != 0:
#  raise IOError(err)
#reads_with_primer = int(result.strip().split()[2])
# Was read 2 specified? Count it too if so
if args.read2 is not None:
  p = subprocess.Popen(['fqgrep', '-m', str(args.mismatches), '-C', '-p', args.primer, args.read1],
                       stdout = subprocess.PIPE, stderr = subprocess.PIPE)
  result, err = p.communicate()
  if p.returncode != 0:
    raise IOError(err)
  reads_with_primer += int(result.strip().split()[2])
#print('Reads with primer: ' + str(reads_with_primer))
print('')

# Filter reads containing the expected Tn end
# Done for the night 11/24 - something's screwy here with this damn pipe
print('### ' + time.strftime("%c") + ' ###')
print('Searching for reads with an inverted repeat in the proper position on the read...')
min_ir_position = len(args.primer) + len(args.invertedrepeat) - 2
max_ir_position = len(args.primer) + len(args.invertedrepeat) + 8
cmd = 'fqgrep -m ' + str(args.mismatches) + ' -r -p ' + args.primer + args.invertedrepeat + ' ' + args.read1 + \
      ' | awk -v min=' + str(min_ir_position) + ' -v max=' + str(max_ir_position) + \
      ' -F "\\t" \'(($8 >= min && $8 <= max) || $1 == "read name")\' | trimmer --5-prime > ' + \
      args.output + '.Tnreads.trimmed.fastq'
#DEBUG: skip computation
#p = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
#result, err = p.communicate()
#if p.returncode != 0:
#  raise IOError(err)
# Was read 2 specified? Include it too if so
if args.read2 is not None:
  cmd = 'fqgrep -m ' + str(args.mismatches) + ' -r -p ' + args.primer + args.invertedrepeat + ' ' + args.read2 + \
        ' | awk -v min=' + str(min_ir_position) + ' -v max=' + str(max_ir_position) + \
        ' -F "\\t" \'(($8 >= min && $8 <= max) || $1 == "read name")\' | trimmer --5-prime >> ' + \
        args.output + '.Tnreads.trimmed.fastq'
  p = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
  result, err = p.communicate()
  if p.returncode != 0:
    raise IOError(err)
p = subprocess.Popen(['wc', '-l', args.output + '.Tnreads.trimmed.fastq'], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
result, err = p.communicate()
if p.returncode != 0:
  raise IOError(err)
print('Tn-adjacent reads to be mapped onto reference : ' + str(int(result.strip().split()[0]) / 4))
print('')

# Map Tn-adjacent reads onto reference
print('### ' + time.strftime("%c") + ' ###')
print('Mapping Tn-adjacent reads onto reference...')
if not os.path.isfile(args.reference + '.1.bt2'):
  p = subprocess.Popen(['bowtie2-build', args.reference, args.reference], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
  result, err = p.communicate()
  if p.returncode != 0:
    raise IOError(err)
#DEBUG: skip computation
#p = subprocess.Popen(['bowtie2', '--end-to-end', '-p', str(args.threads), '-a', '-x', args.reference, '-U',
#                      args.output + '.Tnreads.trimmed.fastq', '-S', args.output + '.sam'],
#                     stdout = subprocess.PIPE, stderr = subprocess.PIPE)
#result, err = p.communicate()
#if p.returncode != 0:
#  raise IOError(err)
#print(err)


## 1/25/17 - at here, have some test lines in a SAM file, next prototype the sorting/tallying lines



# Map and convert - feel free to change bowtie2 parameters yourself
#echo "$PREFIX: Mapping with Bowtie2..."
#echo "Bowtie2 report:" >> $PREFIX-TnSeq.txt
#bowtie2 --end-to-end -p 16 -a -x $BOWTIEREF -U $PREFIX-IR-clip.fastq -S $PREFIX.sam 2>> $PREFIX-TnSeq.txt
#grep '^@' $PREFIX.sam > $PREFIX-mapped.sam
#cat $PREFIX.sam | grep -v '^@' | awk -F "\t" '(and($2, 0x4) != 0x4)' | sort -u -k1,1 >> $PREFIX-mapped.sam
#echo "Number of reads mapping at high enough score:" >> $PREFIX-TnSeq.txt
#cat $PREFIX-mapped.sam | wc -l >> $PREFIX-TnSeq.txt
#echo "$PREFIX: Tallying mapping results..."
#grep -v '^@' $PREFIX-mapped.sam | awk -F "\t" 'and($2, 0x100) != 0x100 {if (and($2, 0x10) != 0x10) print $4; else print $4+length($10)}' | grep '[0-9]' | sort | uniq -c | sort -#n -r > $PREFIX-sites.txt
#echo "Number of insertion sites identified:" >> $PREFIX-TnSeq.txt
#wc -l $PREFIX-sites.txt >> $PREFIX-TnSeq.txt
#echo "Most frequent sites:" >> $PREFIX-TnSeq.txt
#head $PREFIX-sites.txt >> $PREFIX-TnSeq.txt
