#!/usr/bin/env perl
$|=1;
use strict;

# A very simple utility to find pairs of genomes with the same identifier but not similar enough to each other, and pairs of genomes with different identifiers but similar enough to be the same individual.
# Simply pipe into this script the output of comparing fingerprints in one or more sets of genomes.

my $cutoff = 0.75;
while (<>) {
	chomp;
	my($q, $t, $c) = split /\t/;
	(undef, $t) = split /\./, $t if $t =~ /\./;
	if ($q eq $t) {
		print "$q\t$c\n" if $c<$cutoff;
	} else {
		print "$q\t$t\t$c\n" if $c>=$cutoff;
	}
}

