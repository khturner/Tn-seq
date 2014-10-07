#!/usr/bin/env perl

sub usage();

unless (@ARGV == 2 &&
	($boundariesFile = $ARGV[0]) =~ /\.boundaries\.tsv$/ &&
	($sitesFile = $ARGV[1]) =~ /\.sitecounts\.tsv$/) { usage(); }

# Fill an AoA with the boundaries of the bins, pre-filled with 1 read each by the R script.
# This is where read counts per bin will be kept.
open BOUNDARIES, "$boundariesFile" or die "Unable to open file $boundariesFile for reading data: $!\n";
my @boundaries;
$bheader = <BOUNDARIES>;
while (<BOUNDARIES>) {
	chomp;
	@line = split(/\t/);
	push @boundaries, [ @line ];
}

# Create a second AoA of the same size as the boundaries AoA, but pre-filled with 0 sites each.
# This is where independent insertions per bin will be kept.
my @numsites;
for ($i = 0; $i <= $#boundaries; $i++) {
	push @numsites, [ @{$boundaries[$i]} ]
}
for ($i = 0; $i <= $#numsites; $i++) {
	for ($j = 2; $j <= $#{ $numsites[0] }; $j++) {
		$numsites[$i][$j] = 0;
	}
}

# Read through the counts per site file line-by-line and tally read counts and site counts.
open SITES, "$sitesFile" or die "Unable to open file $sitesFile for reading data: $!\n";
my $line_counter = 0;
$sheader = <SITES>;
while (<SITES>) {
	chomp;
	@line = split(/\t/);
	$pos = $line[0]; # consider the next insert position
	if ($line_counter > 0) { $line_counter--; }
	while ($line_counter <= $#boundaries && $boundaries[$line_counter][0] <= $pos) { # is the start of the current boundary less than the position?
		if ($boundaries[$line_counter][1] >= $pos) { # is the end of the current boundary more than the position? If so, winner!
			for($i = 1; $i <= $#line; $i++) {
				$boundaries[$line_counter][$i+1] += $line[$i];
				if ($line[$i] > 0) {
					$numsites[$line_counter][$i+1]++; 
				}
			}
		}
		$line_counter++;
	}
}

# Output read counts per bin.
open BOUT, ">".$boundariesFile.".out" or die "Unable to open file $boundariesFile.out for outputting data: $!\n";
print BOUT $bheader;
open SOUT, ">".$boundariesFile.".numsites.out" or die "Unable to open file $boundariesFile.numsites.out for outputting data: $!\n";
print SOUT $bheader;
for ($i=0; $i<=$#boundaries; $i++) {
	print BOUT join ("\t", @{$boundaries[$i]});
	print BOUT "\n";
	print SOUT join ("\t", @{$numsites[$i]});
	print SOUT "\n";
}
close BOUT;

# Output site counts per bin.

for ($i=0; $i<=$#numsites; $i++) {

}
close SOUT;

sub usage() {
print<<EOF;
TnGeneBin.pl, by Keith H. Turner (khturner@utexas.edu), Sep 17, 2014

This is an accessory script for the TnSeq suite of R scripts, used to handle binning
individual sites and the associated counts of reads per site into predefined bins.
These bins are specified by a start and a stop position, and often denote 3'-trimmed
gene boundaries. However, any set of boundaries can be specified.

usage: $0 (out_pfx).boundaries.tsv (out_pfx).sitecounts.tsv

Arguments:

(out_pfx).boundaries.tsv - starts, stops, and 1's for every bin
(out_pfx).sitecounts.tsv - position and read count info

EOF

exit 1;
}