#!/usr/bin/perl
# GFFTrim - Keith H. Turner - Whiteley Lab
# Takes a GFF file and returns a GFF file with the specified amount trimmed from the 3p or 5p ends
# 11/26/12

sub usage();

my $FivepToTrim = 0;
my $ThreepToTrim = 0;

unless (@ARGV == 3 && 
	($FivepToTrim = $ARGV[0]) =~ /^\d+$/ && 
	($ThreepToTrim = $ARGV[1]) =~ /^\d+$/ && 
	($filename = $ARGV[2]) =~ /^(.+)\.gff$/ &&
	$FivepToTrim >= 0 && $FivepToTrim <= 50 &&
	$ThreepToTrim >= 0 && $ThreepToTrim <= 50) {
	usage();
	}
$prefix = $1;

open GENES, "${prefix}.gff" or die "Unable to open file ${prefix}.gff for reading data: $!\n";
open OUT, ">${prefix}.trunc.gff" or die "Unable to open file ${prefix}.trunc.gff for outputting data: $!\n";
while (<GENES>) {
	chomp;
	@line = split(/\t/);
	unless ($line[0] =~ /^#/) {
		$start = $line[3];
		$end = $line[4];
		$length = $end - $start + 1;
		$FivepOff = int($length * ($FivepToTrim/100));
		$ThreepOff = int($length * ($ThreepToTrim/100));
		$line[3] = $start + $FivepOff;
		$FivepEnd = $line[3] - 1;
		$line[4] = $end - $ThreepOff;
		$ThreepStart = $line[4] + 1;
		print OUT "${line[0]}\t${line[1]}\t${line[2]}\t${line[3]}\t${line[4]}\t${line[5]}\t${line[6]}\t${line[7]}\t${line[8]}\n";
	}
	if ($line[0] =~ /^#/) {
		print OUT $_;
		print OUT "\n";
	}
}
close GENES;
close OUT;

sub usage() {
print<<EOF;
GFFTrim, by Keith H. Turner (khturner@utexas.edu), Nov 26, 2012

This program takes a GFF file corresponding to all features in a genome, creates a file
derived from it containing only those features with type "gene" (suffixed ".genes.gff"),
and creates a second file derived from that first one that truncates the "gene" features
by a specified percentage off of the 5' and 3' ends (suffixed ".trunc.gff").

usage: $0 (5p) (3p) (file).gff

Arguments:

(5p)		An integer [0-50] specifying the percentage of the 5' end to trim off of each gene
(3p)		An integer [0-50] specifying the percentage of the 3' end to trim off of each gene
(file).gff	Your GFF file containing all features of your genome

EOF

exit 1;
}
