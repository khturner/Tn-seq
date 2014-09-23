#!/usr/bin/perl
# TnSeqPileup - Keith H. Turner - Whiteley Lab
# Takes a GFF file and a tab-separated (position-count-count-count...) file and piles up how many reads map to each gene
# 11/26/12

sub usage();

unless (@ARGV == 2 && 
	($GFFFile = $ARGV[0]) =~ /\.gff$/ && 
	($SitesFile = $ARGV[1]) =~ /^sites.smoothed.txt$/) {
		usage();
	}

# Store count data in an AoA
open SITES, "$SitesFile" or die "Unable to open file $SitesFile for reading data: $!\n";
my @counts;
$header = <SITES>;
chomp $header;
@headers = split(/\t/, $header);
foreach $prefix (@headers[1 .. $#headers]) {
	push @prefixes, $prefix;
}
while (<SITES>) {
	chomp;
	@line = split(/\t/);
	push @counts, [ @line ];
}
close SITES;
$lastsiteindex = $#counts;

# Go through each gene and sum the number of reads contained between start and stop
open GFF, "$GFFFile" or die "Unable to open file $GFFFile for reading data: $!\n";
my @outdata;
while (<GFF>) {
	chomp;
	@line = split(/\t/);
	unless ($line[0] =~ /^#/) {
		$start = $line[3];
		$end = $line[4];
		$name = "N/A";
		$locus_tag = "N/A";
		@field8 = split(/;/, $line[8]); # attributes
		foreach $tag (@field8) {
			if ($tag =~ /^Name=(.+)$/) {
				$name = $1;
			}
			if ($tag =~ /^locus_tag=(.+)$/) {
				$locus_tag = $1;
			}
		}
		my @outline = ($locus_tag, $name, ($end-$start+1));
		
		# Hack to speed up the insertion site scanning. Could potentially be faster
		if ($counts[$lastsiteindex*.75][0] < $start) {
			$i = int($lastsiteindex*.75);
		}
		elsif ($counts[$lastsiteindex*.5][0] < $start) {
			$i = int($lastsiteindex*.5);
		}
		elsif ($counts[$lastsiteindex*.25][0] < $start) {
			$i = int($lastsiteindex*.25);
		}
		else {
			$i=0;
		}
		
		# Scan and sum
		push @outline, (0) x ($#{$counts[0]});
		while ($counts[$i][0] <= $end && $i <= $lastsiteindex) {
			if ($counts[$i][0] >= $start) {
				for ($j=1; $j<=$#{$counts[$i]}; $j++) {
					$outline[$j+2] += $counts[$i][$j];
				}
			}
			$i++;
		}
		push @outdata, [ @outline ];
	}
}
close GFF;

# Output data
open OUT, ">genes.pileup.txt" or die "Unable to open file $genes.pileup.txt for outputting data: $!\n";
print OUT "Locus_Tag\tName\tLength\t";
print OUT join ("\t", @prefixes);
print OUT "\n";
for ($i=0; $i<=$#outdata; $i++) {
	print OUT join ("\t", @{$outdata[$i]});
	print OUT "\n";
}
close OUT;

sub usage() {
print<<EOF;
TnSeqPileup, by Keith H. Turner (khturner@utexas.edu), Nov 26, 2012

This program takes a GFF file containing only features from a genome with type "gene" and
a tab-separated value file containing Tn insertion sites and the number of reads at each
site, and outputs a larger tab-separated value file, with extra columns reflecting data
extracted from the GFF file and the total number of reads mapping to each gene.

usage: $0 (genes).gff (sites).txt

Arguments:

(genes).gff\tYour GFF file containing only the genes of your genome
(sites).txt\tA tab-separated value file containing location-reads pairs

EOF

exit 1;
}