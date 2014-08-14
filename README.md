Tn-seq
===========
Custom scripts for analyzing (parsing, mapping, and tallying) Tn-seq reads and
determining differentially abundant transposon insertion mutants.

Copyright (c) 2014 Keith H. Turner, Jake Everett, Urvish Trivedi, Kendra P.
Rumbaugh, and Marvin Whiteley

The scripts contained herein can be used to automatically analyze high-throughput
(Illumina) sequencing reads derived from transposon-genome junctions. First, each
individual dataset is analyzed with TnSeq.sh or TnSeq2.sh, and then a control and
test condition and their specified data files are compared with TnSeqAnalysis.sh.
See below for specific usage details and software dependencies. Please direct any
questions on their use or construction to Keith H. Turner (khturner at utexas.edu).


TnSeq.sh
===========
This script takes two FASTQ files specifying read 1 and read 2 of a paired-end
sequencing run done on a Tn-seq library. This sequencing library should have been
prepared with an end-blind method (i.e. the transposon end could be on either read
1 or on read 2), and should contain a "tag" or "IR" sequence derived from the end
of the transposon to identify which reads are transposon-derived. These reads are
found, trimmed of sequence that will not map to the genome (both transposon- and
sequencing adapter-derived), and mapped to your genome with bowtie2. Finally,
insertion site locations and read counts are tallied. All results and run
information is put in a directory named for your files, and these results can be
fed directly into TnSeqAnalysis.sh (see below).

Usage: ./TnSeq.sh [-p \<primer seq\>] [-i \<IR seq\>] [-a \<assembly\>] [-m \<#\>] \<pfx\>

Arguments:

\<primer seq\> - The sequence of your Tn-seq primer specific to your transposon

\<IR seq\>     - The sequence of the transposon end sequence remaining (for junction
   authentication)

\<assembly\>   - The name of the assembly you're using (e.g. "PAO1")

-m \<#\>     - The number of mismatches/indels you want to tolerate during search

\<pfx\>        - the file prefix for your sequence files (If your sequence files are
   named condition1_R1.fastq and condition1_R2.fastq, the prefix is "condition1")

Dependencies:

-fqgrep (https://github.com/indraniel/fqgrep)

-bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

-$REFGENOME defined in your environment. This should specify a directory that
   contains bowtie2 references and your genome annotations in GFF format in a
   subdirectory with the same name as the genome.
   (e.g. if $REFGENOME is "/home/user/ref_genome", then files such as
   "/home/user/ref_genome/PAO1/PAO1.gff", "/home/user/ref_genome/PAO1/PAO1.1.bt2",
   etc. should be present)

-trimmer (this package) should be available through your $PATH

-flexbar (http://sourceforge.net/projects/flexbar/)

-the file "~/adapters/3_adapter_seq.fasta" should specify the Illumina adapters:

\>index_sp

AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC

\>3_adapter_seq

TCGTATGCCGTCTTCTGCTTG


TnSeq2.sh
===========
This script takes one FASTQ file specifying read 1 single-end sequencing run done
on a Tn-seq library. This sequencing library should have been prepared with an
end-specific method (i.e. the transposon end should be on read 1), and should
contain a "tag" or "IR" sequence derived from the end of the transposon to identify
which reads are transposon-derived. These reads are found, trimmed of sequence that
will not map to the genome (both transposon- and sequencing adapter-derived), and
mapped to your genome with bowtie2. Finally, insertion site locations and read
counts are tallied. All results and run information is put in a directory named for
your files, and these results can be fed directly into TnSeqAnalysis.sh (see below).

Usage: ./TnSeq2.sh [-p \<primer seq\>] [-i \<IR seq\>] [-a \<assembly\>] [-m \<#\>] \<pfx\>

Arguments:

\<primer seq\> - The sequence of your Tn-seq primer specific to your transposon

\<IR seq\>     - The sequence of the transposon end sequence remaining (for junction
   authentication)

\<assembly\>   - The name of the assembly you're using (e.g. "PAO1")

-m \<#\>     - The number of mismatches/indels you want to tolerate during search

\<pfx\>        - the file prefix for your sequence files (If your sequence file is
   named condition1_R1.fastq, the prefix is "condition1")

Dependencies:

-fqgrep (https://github.com/indraniel/fqgrep)

-bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

-$REFGENOME defined in your environment. This should specify a directory that
   contains bowtie2 references and your genome annotations in GFF format in a
   subdirectory with the same name as the genome.
   (e.g. if $REFGENOME is "/home/user/ref_genome", then files such as
   "/home/user/ref_genome/PAO1/PAO1.gff", "/home/user/ref_genome/PAO1/PAO1.1.bt2",
   etc. should be present)

-trimmer (this package) should be available through your $PATH

-flexbar (http://sourceforge.net/projects/flexbar/)

-the file "~/adapters/3_adapter_seq.fasta" should specify the Illumina adapters:

\>index_sp

AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC

\>3_adapter_seq

TCGTATGCCGTCTTCTGCTTG


GFFTrim.pl
===========
This program takes a GFF file on STDIN, trims the specified percentage from the 3'
and 5' ends of every feature in that GFF file, and returns a GFF file suffixed
".trunc.gff" with these new starts/ends. This works best with this software package
if only features marked as "gene" are present in the GFF file.

Usage: ./GFFTrim.pl (5p) (3p) (genome).gff

Arguments:

(5p)         - An integer [0-50] specifying the percentage of the 5' end to trim
   off of each feature

(3p)         - An integer [0-50] specifying the percentage of the 3' end to trim
   off of each feature

(genome).gff - The GFF file to be trimmed


PullKegg.pl
===========
This program takes a GFF file on STDIN, searches http://www.kegg.jp for the
locus_tag and returns a new GFF file with KO and Kegg pathway information on
STDOUT. This works best with this software package if only features marked as
"gene" are present in the GFF file.

Usage: ./PullKegg.pl < (in).gff > (out.gff)

Dependencies:

-wget (Unix)


TnSeqAnalysis.sh
===========
This script takes the results of TnSeq.sh or TnSeq2.sh for sequence files derived
from one or more replicates of two conditions that you wish to compare for
differential mutant abundance. The insertion site locations and read counts data is
smoothed with LOESS smoothing (to correct for genomic position-dependent effects on
apparent insertion abundance), normalized with DESeq, the number of reads per gene
and the number of independent insertions identified per gene is tallied, and
differential mutant abundance is calculated using a negative binomial test.

Usage: ./TnSeqAnalysis.sh [-i \<#\>] [-a \<assembly\>] [-o \<output\>] [-c \<name\>]
   [-x \<#\>] [-t \<name\>] [-y \<#\>] \<pfx1\> \<pfx2\> \<pfx3\> ... \<pfxn\>

Arguments:

-i \<#\>     - The number of the most represented sites to ignore

\<assembly\> - The name of the assembly you're using (e.g. "PAO1")

\<output\>   - The name for the output file

-c \<name\>  - The name for the control condition

-x \<#\>     - The number of replicates for the control condition

-t \<name\>  - The name for the test condition

-y \<#\>     - The number of replicates for the test condition

\<pfx#\>     - The file prefixes to be considered, listed with the control conditions
   followed by the test conditions (e.g. "./TnSeqAnalysis.sh -i 50 -a PAO1
   -o Example -c control -x 2 -t test -y 2 C1 C2 T1 T2")

Dependencies:

-$REFGENOME defined in your environment. This should specify a directory that
   contains your genome annotations in GFF format in a subdirectory with the same
   name as the genome. (e.g. if $REFGENOME is "/home/user/ref_genome", then the
   files "/home/user/ref_genome/PAO1/PAO1.trunc.gff" and 
   "/home/user/ref_genome/PAO1/PAO1.gene.products.kegg.txt" should be present)
   
-R (http://www.r-project.org/)

-TnSeqDESeq.R (this package) placed in ~/local/bin/

-DESeq (CRAN)

-(assembly).gene.products.kegg.txt, a tab-separated file containing any annotation
   information you want automatically appended to your DESeq results file. This is
   so named because our file includes (locus)-(gene name)-(product)-(KEGG
   Orthology number)-(KEGG Pathway) information, but you can include whatever you'd
   like. If you want to change the name of this file, also do so on TnSeqDESeq.R
   line 69. This file should be in the location described above.
